#!/usr/bin/env python3
"""
Abundance of every gene of the non-redundant gene catalog in every sample.

Each gene inherits the mean depth of its contig, already computed for binning
(jgi_summarize_bam_contig_depths, one table per assembly). No reads are mapped
again.

How the depth of a catalog gene (a cluster) is built, for each sample:

  1. Within one assembly, the depths of the cluster members are summed. Reads
     of a gene with several variants in the same assembly are split between
     them, so the sum recovers the depth of the gene.
  2. Across assemblies, those sums are averaged. With --map_all_samples (or a
     single co-assembly), every sample is mapped against every assembly, so the
     same reads are counted once per assembly that contains the gene. Summing
     across assemblies would inflate genes assembled many times (common,
     abundant genes) over genes assembled only once.

     Only assemblies with a depth column for that sample count: with
     --assembly_mode group or per_sample without --map_all_samples, a sample is
     mapped only against its own assembly.

TPM: mean depth is already reads x read length / contig length, the analogue of
reads per kilobase. TPM is the depth of each gene divided by the sum of the
depths of all genes in that sample, times 10^6. No further length correction.

Gene names are <assembly>|<contig>_<n> (Prodigal names, prefixed by the
assembly in PRODIGAL_CONTIGS). Depth tables are named <assembly>.depth.txt.

Usage:
  gene_abundance.py --clusters gene_catalog.clusters.tsv.gz \
                    --depths T0.depth.txt T1.depth.txt ... \
                    --out-depth gene_depth.tsv.gz --out-tpm gene_tpm.tsv.gz \
                    --out-stats gene_abundance_stats.tsv
"""

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_clusters(path):
    """Representative index of every member, grouped by assembly."""
    rep_index = {}
    rep_names = []
    by_assembly = {}  # assembly -> (list of contigs, list of rep indices)

    with open_text(path) as fh:
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")
            ri = rep_index.get(rep)
            if ri is None:
                ri = len(rep_names)
                rep_index[rep] = ri
                rep_names.append(rep)

            assembly, gene = member.split("|", 1)
            contig = gene.rsplit("_", 1)[0]
            contigs, reps = by_assembly.setdefault(assembly, ([], []))
            contigs.append(contig)
            reps.append(ri)

    return rep_names, by_assembly


def read_depth(path):
    """Mean depth per contig and sample: contigName, then <sample>.sorted.bam columns."""
    header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    mean_cols = [c for c in header[3:] if not c.endswith("-var")]
    table = pd.read_csv(path, sep="\t", usecols=["contigName"] + mean_cols,
                        index_col="contigName")
    samples = [c[: -len(".sorted.bam")] if c.endswith(".sorted.bam") else c
               for c in mean_cols]
    return table, samples


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clusters", required=True,
                   help="MMseqs2 cluster table: representative<TAB>member")
    p.add_argument("--depths", required=True, nargs="+",
                   help="<assembly>.depth.txt tables from MERGE_DEPTHS")
    p.add_argument("--out-depth", required=True)
    p.add_argument("--out-tpm", required=True)
    p.add_argument("--out-stats", required=True)
    args = p.parse_args()

    rep_names, by_assembly = read_clusters(args.clusters)
    n_reps = len(rep_names)
    print(f"{n_reps} catalog genes, members in {len(by_assembly)} assemblies",
          file=sys.stderr)

    # Sample columns of every assembly, and their union in a fixed order
    depth_files = {}
    for path in args.depths:
        name = os.path.basename(path)
        if not name.endswith(".depth.txt"):
            sys.exit(f"Depth tables must be named <assembly>.depth.txt: {name}")
        depth_files[name[: -len(".depth.txt")]] = path

    missing = sorted(set(by_assembly) - set(depth_files))
    if missing:
        sys.exit(f"No depth table for assemblies: {', '.join(missing)}")

    tables = {}
    all_samples = set()
    for assembly in sorted(by_assembly):
        table, samples = read_depth(depth_files[assembly])
        tables[assembly] = (table, samples)
        all_samples.update(samples)
    samples = sorted(all_samples)
    col = {s: i for i, s in enumerate(samples)}

    sums = np.zeros((n_reps, len(samples)), dtype=np.float64)
    n_asm = np.zeros((n_reps, len(samples)), dtype=np.uint16)

    for assembly in sorted(by_assembly):
        contigs, reps = by_assembly[assembly]
        table, asm_samples = tables[assembly]
        cols = np.array([col[s] for s in asm_samples])

        rows = table.index.get_indexer(contigs)
        if (rows < 0).any():
            n_bad = int((rows < 0).sum())
            sys.exit(f"{assembly}: {n_bad} genes sit on contigs absent from its depth table")

        reps = np.asarray(reps)
        values = table.to_numpy(dtype=np.float64)[rows]
        # step 1: sum the members of each cluster within this assembly
        np.add.at(sums, (reps[:, None], cols[None, :]), values)
        # step 2: count this assembly once per cluster, for the samples it has
        present = np.unique(reps)
        n_asm[np.ix_(present, cols)] += 1
        print(f"{assembly}: {len(contigs)} genes, {len(present)} catalog genes, "
              f"{len(asm_samples)} samples", file=sys.stderr)
        del tables[assembly]

    # In place, to keep a single genes x samples matrix in memory: sums
    # becomes the mean depth, then the TPM. Where no assembly counted, the sum
    # is already 0.
    depth = sums
    np.divide(depth, n_asm, out=depth, where=n_asm > 0)
    del n_asm

    gz = {"method": "gzip", "compresslevel": 1}
    index = pd.Index(rep_names, name="gene")
    pd.DataFrame(depth, index=index, columns=samples, copy=False).to_csv(
        args.out_depth, sep="\t", float_format="%.6g", compression=gz)

    totals = depth.sum(axis=0)
    detected = (depth > 0).sum(axis=0)
    scale = np.divide(1e6, totals, out=np.zeros_like(totals), where=totals > 0)
    tpm = depth
    tpm *= scale
    pd.DataFrame(tpm, index=index, columns=samples, copy=False).to_csv(
        args.out_tpm, sep="\t", float_format="%.6g", compression=gz)

    pd.DataFrame({
        "sample": samples,
        "n_genes_detected": detected,
        "total_depth": totals,
    }).to_csv(args.out_stats, sep="\t", index=False, float_format="%.6g")


if __name__ == "__main__":
    main()
