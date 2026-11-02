#!/usr/bin/env python3
"""
Community profiles from the non-redundant gene catalog.

Joins the TPM of every catalog gene (gene_abundance.py) with its annotations
and sums the TPM of the genes of each feature in each sample:

  ko_per_sample.tsv        KEGG orthologue (eggNOG-mapper). A gene with several
                           KOs counts towards each of them.
  canthyd_per_sample.tsv   CANT-HYD marker (one per gene, best model)
  hadeg_per_sample.tsv     HADEG gene, with its compound, pathway and subpathway
  taxonomy_per_sample.tsv  taxon at each rank, domain to species. A gene takes
                           the CAT classification of its contig; genes on
                           contigs without a classification at that rank are
                           counted as "Unclassified".

Long tables: one row per feature and sample, with n_genes (catalog genes with
that feature) and tpm. Also written:

  catalog_annotation.tsv.gz  one row per catalog gene: KOs, CANT-HYD and HADEG
                             markers, and the lineage of its contig
  contig_taxonomy.tsv        CAT classification of every contig of every assembly
  profile_stats.tsv          fraction of TPM annotated by each source, per sample

Files named NO_FILE are ignored, and a table whose source is missing is not
written. Gene names are <assembly>|<contig>_<n>; CAT tables are named
<assembly>.CAT.named.txt.
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
UNCLASSIFIED = "Unclassified"


def real_files(paths):
    return [p for p in (paths or []) if os.path.basename(p) != "NO_FILE" and os.path.getsize(p) > 0]


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_eggnog(paths):
    """gene -> list of KOs"""
    kos = {}
    for path in paths:
        c_query = c_ko = None
        with open_text(path) as fh:
            for line in fh:
                if line.startswith("#query"):
                    header = line.rstrip("\n").split("\t")
                    c_query, c_ko = header.index("#query"), header.index("KEGG_ko")
                    continue
                if line.startswith("#") or c_query is None:
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) <= c_ko or f[c_ko] in ("-", ""):
                    continue
                ko = [k.replace("ko:", "") for k in f[c_ko].split(",") if k]
                if ko:
                    kos[f[c_query]] = ko
    return kos


def read_table(paths, columns):
    """protein -> tuple of columns, from the CANTHYD / HADEG tables"""
    out = {}
    for path in paths:
        t = pd.read_csv(path, sep="\t", dtype=str)
        for row in t[["protein"] + columns].itertuples(index=False):
            out[row[0]] = tuple(row[1:])
    return out


def clean_taxon(value):
    """'g__Alcanivorax: 0.95' -> 'g__Alcanivorax'; 'no support' or empty -> None"""
    value = value.strip()
    if not value or value in ("no support", "NA"):
        return None
    name = value.rsplit(": ", 1)[0] if ": " in value else value
    return name.rstrip("*").strip() or None


def read_cat(paths):
    """(assembly, contig) -> tuple of taxa (domain..species, None if unknown), and a long table"""
    taxonomy = {}
    rows = []
    for path in paths:
        name = os.path.basename(path)
        assembly = name.split(".CAT.")[0]
        with open(path) as fh:
            header = fh.readline().rstrip("\n").lstrip("# ").split("\t")
            rank_cols = []
            for rank in RANKS:
                key = "superkingdom" if rank == "domain" else rank
                rank_cols.append(header.index(key) if key in header else None)
            c_class = header.index("classification")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                taxa = tuple(clean_taxon(f[c]) if c is not None and c < len(f) else None
                             for c in rank_cols)
                taxonomy[(assembly, f[0])] = taxa
                rows.append((assembly, f[0], f[c_class]) + tuple(t or "" for t in taxa))
    table = pd.DataFrame(rows, columns=["assembly", "contig", "classification"] + RANKS)
    return taxonomy, table


def gene_contig(gene):
    assembly, rest = gene.split("|", 1)
    return assembly, rest.rsplit("_", 1)[0]


class Profile:
    """Sum of TPM per feature (a tuple) and sample, plus the genes of each feature."""

    def __init__(self, n_samples):
        self.n_samples = n_samples
        self.tpm = defaultdict(lambda: np.zeros(n_samples))
        self.n_genes = defaultdict(int)

    def add(self, keys, values):
        """keys: one feature per row of values (genes x samples)"""
        if not keys:
            return
        codes, uniques = pd.factorize(pd.Series(keys, dtype=object))
        sums = np.zeros((len(uniques), self.n_samples))
        np.add.at(sums, codes, values)
        counts = np.bincount(codes, minlength=len(uniques))
        for i, key in enumerate(uniques):
            self.tpm[key] += sums[i]
            self.n_genes[key] += int(counts[i])

    def write(self, path, columns, samples):
        keys = sorted(self.tpm, key=lambda k: tuple(str(x) for x in k))
        out = []
        for key in keys:
            for j, sample in enumerate(samples):
                out.append(tuple(key) + (sample, self.n_genes[key], self.tpm[key][j]))
        pd.DataFrame(out, columns=columns + ["sample", "n_genes", "tpm"]).to_csv(
            path, sep="\t", index=False, float_format="%.6g")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tpm", required=True, help="gene_tpm.tsv.gz from gene_abundance.py")
    p.add_argument("--eggnog", nargs="*", default=[])
    p.add_argument("--canthyd", nargs="*", default=[])
    p.add_argument("--hadeg", nargs="*", default=[])
    p.add_argument("--cat", nargs="*", default=[])
    p.add_argument("--outdir", default=".")
    p.add_argument("--chunksize", type=int, default=500000)
    args = p.parse_args()

    eggnog_files, canthyd_files = real_files(args.eggnog), real_files(args.canthyd)
    hadeg_files, cat_files = real_files(args.hadeg), real_files(args.cat)

    kos = read_eggnog(eggnog_files) if eggnog_files else None
    canthyd = read_table(canthyd_files, ["marker", "description"]) if canthyd_files else None
    hadeg = read_table(hadeg_files, ["marker", "compound", "pathway", "subpathway"]) if hadeg_files else None
    taxonomy, contig_table = read_cat(cat_files) if cat_files else (None, None)

    def out(name):
        return os.path.join(args.outdir, name)

    if contig_table is not None:
        contig_table.to_csv(out("contig_taxonomy.tsv"), sep="\t", index=False)

    samples = pd.read_csv(args.tpm, sep="\t", nrows=0).columns[1:].tolist()
    n = len(samples)
    prof_ko, prof_canthyd, prof_hadeg = Profile(n), Profile(n), Profile(n)
    prof_tax = Profile(n)
    total = np.zeros(n)
    annotated = {k: np.zeros(n) for k in ("eggnog_ko", "canthyd", "hadeg", "cat_domain")}

    gz = {"method": "gzip", "compresslevel": 1}
    first = True
    for chunk in pd.read_csv(args.tpm, sep="\t", index_col=0, chunksize=args.chunksize):
        genes = chunk.index.tolist()
        values = chunk.to_numpy(dtype=np.float64)
        total += values.sum(axis=0)
        ann = pd.DataFrame(index=chunk.index)

        if kos is not None:
            rows, keys = [], []
            for i, g in enumerate(genes):
                for ko in kos.get(g, ()):
                    rows.append(i)
                    keys.append((ko,))
            prof_ko.add(keys, values[rows])
            has = np.array([g in kos for g in genes], dtype=bool)
            annotated["eggnog_ko"] += values[has].sum(axis=0)
            ann["ko"] = [",".join(kos.get(g, ())) for g in genes]

        for name, table, prof in (("canthyd", canthyd, prof_canthyd), ("hadeg", hadeg, prof_hadeg)):
            if table is None:
                continue
            rows = [i for i, g in enumerate(genes) if g in table]
            prof.add([table[genes[i]] for i in rows], values[rows])
            annotated[name] += values[rows].sum(axis=0)
            ann[f"{name}_marker"] = [table[g][0] if g in table else "" for g in genes]
            if name == "hadeg":
                ann["hadeg_pathway"] = [table[g][2] if g in table else "" for g in genes]

        if taxonomy is not None:
            lineages = [taxonomy.get(gene_contig(g), (None,) * len(RANKS)) for g in genes]
            for r, rank in enumerate(RANKS):
                prof_tax.add([(rank, lin[r] or UNCLASSIFIED) for lin in lineages], values)
            has = np.array([lin[0] is not None for lin in lineages], dtype=bool)
            annotated["cat_domain"] += values[has].sum(axis=0)
            ann["lineage"] = [";".join(t for t in lin if t) for lin in lineages]

        ann.to_csv(out("catalog_annotation.tsv.gz"), sep="\t", mode="w" if first else "a",
                   header=first, compression=gz)
        first = False

    if kos is not None:
        prof_ko.write(out("ko_per_sample.tsv"), ["ko"], samples)
    if canthyd is not None:
        prof_canthyd.write(out("canthyd_per_sample.tsv"), ["marker", "description"], samples)
    if hadeg is not None:
        prof_hadeg.write(out("hadeg_per_sample.tsv"),
                         ["marker", "compound", "pathway", "subpathway"], samples)
    if taxonomy is not None:
        prof_tax.write(out("taxonomy_per_sample.tsv"), ["rank", "taxon"], samples)

    # Share of the TPM of each sample that each source annotates
    stats = pd.DataFrame({"sample": samples, "total_tpm": total})
    used = {"eggnog_ko": kos, "canthyd": canthyd, "hadeg": hadeg, "cat_domain": taxonomy}
    for key, source in used.items():
        if source is not None:
            stats[f"fraction_{key}"] = np.divide(annotated[key], total,
                                                 out=np.zeros(n), where=total > 0)
    stats.to_csv(out("profile_stats.tsv"), sep="\t", index=False, float_format="%.4g")


if __name__ == "__main__":
    main()
