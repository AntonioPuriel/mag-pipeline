#!/usr/bin/env python3
"""Classify genomes by the MIMAG standard (Bowers et al. 2017, doi:10.1038/nbt.3893).

  high    completeness > 90 %, contamination < 5 %, 5S, 16S and 23S rRNA genes
          present, and tRNAs for at least 18 of the 20 standard amino acids
  medium  completeness >= 50 %, contamination < 10 %
  low     anything below medium

A genome that meets the high thresholds of CheckM2 but lacks an rRNA gene or
enough tRNAs is medium. The column `limiting` says what keeps each genome below
high quality.

rRNA and tRNA genes are searched with both the bacterial and the archaeal
models. The domain comes from GTDB-Tk when its classification is given;
otherwise the model that finds more rRNA types, then more tRNA types, is kept
(bacterial on ties).

One row per genome in the rRNA tables, which only list the genomes that were
searched (the dereplicated representatives with --dereplicate).

Usage:
  mimag_quality.py --bin-summary bin_summary.tsv --checkm2 C [C ...]
                   --rrna R [R ...] --trna T [T ...]
                   [--gunc G [G ...]] [--gtdbtk X [X ...]]
Files named NO_FILE are ignored.
"""
import argparse
import csv
import os
import sys

RRNA_TYPES = ("5S", "16S", "23S")
MIN_TRNA_TYPES = 18


def real_files(paths):
    return [p for p in (paths or []) if os.path.basename(p) != "NO_FILE" and os.path.getsize(p) > 0]


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def per_model(paths, columns):
    """genome -> model -> {column: value} from the rRNA or tRNA summaries."""
    table = {}
    for path in real_files(paths):
        for row in read_tsv(path):
            table.setdefault(row["genome"], {})[row["model"]] = {c: row[c] for c in columns}
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-summary", required=True)
    parser.add_argument("--checkm2", nargs="+", required=True)
    parser.add_argument("--rrna", nargs="+", required=True)
    parser.add_argument("--trna", nargs="+", required=True)
    parser.add_argument("--gunc", nargs="*", default=[])
    parser.add_argument("--gtdbtk", nargs="*", default=[])
    args = parser.parse_args()

    # Size and contiguity. A bin that DAS_Tool kept unchanged appears twice
    # (under its binner and under dastool) with the same contigs: either row works
    stats = {}
    for path in real_files([args.bin_summary]):
        for row in read_tsv(path):
            if row["bin"] not in stats or row["binner"] == "dastool":
                stats[row["bin"]] = row

    checkm2 = {}
    for path in real_files(args.checkm2):
        for row in read_tsv(path):
            checkm2[row["Name"]] = (row.get("Completeness"), row.get("Contamination"))

    gunc = {}
    for path in real_files(args.gunc):
        for row in read_tsv(path):
            flag = row.get("pass.GUNC")
            gunc[row["genome"]] = "pass" if flag == "True" else "fail" if flag == "False" else "NA"

    taxonomy = {}
    for path in real_files(args.gtdbtk):
        for row in read_tsv(path):
            taxonomy[row["user_genome"]] = (row.get("classification", ""), row.get("closest_genome_ani", "NA"))

    rrna = per_model(args.rrna, ["rrna_5S", "rrna_16S", "rrna_23S", "rrna_partial"])
    trna = per_model(args.trna, ["trna_types", "trna_total", "trna_missing"])

    header = [
        "genome", "assembly", "binner", "n_contigs", "total_length", "n50", "gc_pct",
        "completeness", "contamination", "gunc",
        "domain", "domain_source",
        "rrna_5S", "rrna_16S", "rrna_23S", "rrna_partial",
        "trna_types", "trna_total", "trna_missing",
        "mimag_quality", "limiting", "classification", "closest_ani",
    ]
    rows = []

    for genome in sorted(rrna):
        classification, closest = taxonomy.get(genome, ("", "NA"))
        domain_tag = classification.split(";")[0] if classification else ""

        def score(model):
            r = rrna[genome].get(model, {})
            t = trna.get(genome, {}).get(model, {})
            rrna_types = sum(1 for k in RRNA_TYPES if int(r.get("rrna_" + k, 0)) > 0)
            return (rrna_types, int(t.get("trna_types", 0)),
                    sum(int(r.get("rrna_" + k, 0)) for k in RRNA_TYPES), int(t.get("trna_total", 0)))

        if domain_tag == "d__Archaea":
            model, source = "arc", "gtdbtk"
        elif domain_tag == "d__Bacteria":
            model, source = "bac", "gtdbtk"
        else:
            # max() keeps the first of equal scores: bacterial on ties
            model, source = max(("bac", "arc"), key=score), "rrna_trna"

        r = rrna[genome].get(model, {})
        t = trna.get(genome, {}).get(model, {"trna_types": "0", "trna_total": "0", "trna_missing": "NA"})

        comp, cont = checkm2.get(genome, (None, None))
        c, k = number(comp), number(cont)

        limiting = []
        if c is None or k is None:
            quality = "NA"
            limiting.append("no CheckM2 result")
        else:
            if c <= 90:
                limiting.append("completeness <= 90")
            if k >= 5:
                limiting.append("contamination >= 5")
            missing_rrna = [x for x in RRNA_TYPES if int(r.get("rrna_" + x, 0)) == 0]
            if missing_rrna:
                limiting.append("no " + "/".join(missing_rrna) + " rRNA")
            if int(t["trna_types"]) < MIN_TRNA_TYPES:
                limiting.append("tRNA types %s < %d" % (t["trna_types"], MIN_TRNA_TYPES))

            if not limiting:
                quality = "high"
            elif c >= 50 and k < 10:
                quality = "medium"
            else:
                quality = "low"

        s = stats.get(genome, {})
        rows.append([
            genome, s.get("assembly", "NA"), s.get("binner", "NA"),
            s.get("n_contigs", "NA"), s.get("total_length", "NA"), s.get("n50", "NA"), s.get("gc_pct", "NA"),
            comp if comp is not None else "NA", cont if cont is not None else "NA",
            gunc.get(genome, "NA"),
            "Archaea" if model == "arc" else "Bacteria", source,
            r.get("rrna_5S", "0"), r.get("rrna_16S", "0"), r.get("rrna_23S", "0"), r.get("rrna_partial", "0"),
            t["trna_types"], t["trna_total"], t["trna_missing"],
            quality, "; ".join(limiting) if limiting else "-",
            classification or "NA", closest or "NA",
        ])

    order = {"high": 0, "medium": 1, "low": 2, "NA": 3}
    rows.sort(key=lambda row: (order[row[19]], -(number(row[7]) or 0), row[0]))

    out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    out.writerow(header)
    out.writerows(rows)

    counts = {q: sum(1 for row in rows if row[19] == q) for q in order}
    sys.stderr.write("MIMAG quality: %d high, %d medium, %d low, %d without CheckM2\n"
                     % (counts["high"], counts["medium"], counts["low"], counts["NA"]))


if __name__ == "__main__":
    main()
