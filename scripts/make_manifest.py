#!/usr/bin/env python3
"""Write manifest.json and report/ for the results of an existing run.

Since 2.1.0, mag-pipeline writes manifest.json and the report/ folder at the
end of every run (workflow.onComplete in nextflow.config). This script builds
the same two things for results that already exist, without running the
pipeline again: runs of 2.0.0 or earlier, or runs whose cache is gone.

The pipeline, run, database and parameter fields come from run_info.txt,
which every run since 1.6 writes. The versions of the tools cannot be
recovered after the run, so `software` is left empty. The manifest records
that it was rebuilt in `source`.

Usage:
    python3 scripts/make_manifest.py RESULTS_DIR [--force] [--no-report]

Only the Python standard library is needed (Python 3.6 or later).
"""

import argparse
import datetime
import gzip
import json
import os
import re
import shutil
import sys

SCHEMA_VERSION = "1.1"

# table name -> (folder under the results, file name pattern). Keep in sync
# with `known` in workflow.onComplete (nextflow.config):
# tests/test_make_manifest.sh checks that both lists are the same.
TABLES = [
    ("assembly_stats",       "02_assembly",              r".+\.contigs\.stats\.tsv"),
    ("mapping_summary",      "03_mapping",               r"mapping_summary\.tsv"),
    ("depth",                "03_mapping",               r".+\.depth\.txt"),
    ("bin_summary",          "04_binning",               r"bin_summary\.tsv"),
    ("contig2bin",           "04_binning",               r"contig2bin\.tsv"),
    ("mag_abundance",        "04_binning",               r"mag_abundance\.tsv(\.gz)?"),
    ("dastool_summary",      "04_binning/dastool",       r".+_DASTool_summary\.tsv"),
    ("mag_quality",          "05_quality",               r"mag_quality\.tsv"),
    ("mimag_quality",        "05_quality",               r"mimag_quality\.tsv"),
    ("checkm2",              "05_quality/checkm2",       r".+\.checkm2\.tsv"),
    ("gunc",                 "05_quality/gunc",          r".+\.gunc\.tsv"),
    ("dereplication",        "05_quality/dereplication", r"dereplication_clusters\.tsv"),
    ("gtdbtk",               "06_taxonomy",              r".+\.gtdbtk\.tsv"),
    ("mag_tree",             "06_taxonomy/mag_tree",     r".+\.tree"),
    ("ko_per_mag",           "07_function",              r"ko_per_mag\.tsv"),
    ("annotation_stats",     "07_function",              r"annotation_stats\.tsv"),
    ("hydrocarbon_per_mag",  "07_function",              r"hydrocarbon_per_mag\.tsv"),
    ("gene_catalog_stats",   "08_contigs/catalog",       r"gene_catalog\.stats\.tsv"),
    ("gene_abundance_stats", "08_contigs/abundance",     r"gene_abundance_stats\.tsv"),
    ("profile_stats",        "08_contigs/profiles",      r"profile_stats\.tsv"),
    ("taxonomy_per_sample",  "08_contigs/profiles",      r"taxonomy_per_sample\.tsv(\.gz)?"),
    ("ko_per_sample",        "08_contigs/profiles",      r"ko_per_sample\.tsv(\.gz)?"),
    ("canthyd_per_sample",   "08_contigs/profiles",      r"canthyd_per_sample\.tsv(\.gz)?"),
    ("hadeg_per_sample",     "08_contigs/profiles",      r"hadeg_per_sample\.tsv(\.gz)?"),
    ("singlem",              "09_singlem",               r"singlem_profiles_with_extras\.tsv(\.gz)?"),
    ("singlem_per_sample",   "09_singlem/profiles",      r".+\.singlem\.tsv(\.gz)?"),
    ("mag_tree_msa",         "06_taxonomy/mag_tree",     r".+\.user_msa\.fasta(\.gz)?"),
]

# table name -> (folder under the results, pattern of the path below it):
# files in subfolders. Keep in sync with `nested` in workflow.onComplete.
NESTED = [
    ("gtdbtk_msa",           "06_taxonomy/gtdbtk",       r"[^/]+/align/[^/]+\.user_msa\.fasta(\.gz)?"),
    ("gtdbtk_trees",         "06_taxonomy/gtdbtk",       r"[^/]+/classify/[^/]+\.tree"),
]

# Left out of report/: one row per contig (GB in large runs), and the GTDB-Tk
# alignments and placement trees (thousands of reference genomes)
EXCLUDED_FROM_REPORT = ["depth", "contig2bin", "gtdbtk_msa", "gtdbtk_trees"]

# Labels of run_info.txt -> keys of `databases` in the manifest
DATABASE_KEYS = {
    "CheckM2": "checkm2", "GTDB-Tk": "gtdbtk", "eggNOG": "eggnog",
    "CANT-HYD": "canthyd", "HADEG": "hadeg", "CAT": "cat", "SingleM": "singlem",
    "GUNC": "gunc",
}


def find_tables(outdir):
    tables = {}
    for name, folder, pattern in TABLES:
        directory = os.path.join(outdir, folder)
        if not os.path.isdir(directory):
            continue
        hits = sorted(
            "{}/{}".format(folder, f) for f in os.listdir(directory)
            if re.fullmatch(pattern, f) and os.path.isfile(os.path.join(directory, f))
        )
        if hits:
            tables[name] = hits
    for name, folder, pattern in NESTED:
        directory = os.path.join(outdir, folder)
        hits = []
        for root, _, files in os.walk(directory):
            for f in files:
                rel = os.path.relpath(os.path.join(root, f), directory).replace(os.sep, "/")
                if re.fullmatch(pattern, rel):
                    hits.append("{}/{}".format(folder, rel))
        if hits:
            tables[name] = sorted(hits)
    return tables


def parse_value(text):
    """Parameter values of run_info.txt, typed as Nextflow would write them."""
    if text in ("null", ""):
        return None
    if text in ("true", "false"):
        return text == "true"
    if re.fullmatch(r"-?\d+", text):
        return int(text)
    if re.fullmatch(r"-?(\d+\.\d*|\.\d+)([eE]-?\d+)?|-?\d+[eE]-?\d+", text):
        return float(text)
    return text


def read_run_info(path):
    header, command, databases, params = {}, [], {}, {}
    section = "header"
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if re.fullmatch(r"[a-z][a-z ]*:", line):
                section = line[:-1]
                continue
            if not line.strip() or line == "mag-pipeline run information":
                continue
            if section == "header" and ":" in line:
                key, value = line.split(":", 1)
                header[key.strip()] = value.strip()
            elif section == "command line":
                command.append(line.strip())
            elif section == "reference databases" and ":" in line:
                label, value = (part.strip() for part in line.split(":", 1))
                if value == "not used":
                    continue
                match = re.fullmatch(r"(.*?)\s*\(last modified (.*)\)", value)
                path_, modified = (match.group(1), match.group(2)) if match else (value, None)
                if modified and not re.fullmatch(r"\d{4}-\d{2}-\d{2}", modified):
                    modified = None
                key = DATABASE_KEYS.get(label, re.sub(r"[^a-z0-9]", "", label.lower()))
                databases[key] = {"path": path_, "modified": modified}
            elif section == "parameters" and "=" in line:
                key, value = line.split("=", 1)
                params[key.strip()] = parse_value(value.strip())

    pipeline = header.get("pipeline", "")
    name, _, version = pipeline.rpartition(" ")
    return {
        "pipeline": {
            "name": name or pipeline or None,
            "version": version or None,
            "revision": None,
            "commit": None,
            "nextflow": header.get("nextflow"),
        },
        "run": {
            "name": None,
            "session_id": None,
            "started": header.get("started"),
            "completed": header.get("completed"),
            "duration": header.get("duration"),
            "success": parse_value(header.get("success", "")),
            "exit_status": None,
            "profile": header.get("profile"),
            "container_engine": header.get("container engine"),
            "command_line": " ".join(command) or None,
        },
        "databases": databases,
        "params": dict(sorted(params.items())),
    }


def write_json(path, data):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=4)
        handle.write("\n")


def build_report(outdir, manifest):
    """report/: the summary tables gzipped, as workflow.onComplete does."""
    report = os.path.join(outdir, "report")
    if os.path.isdir(report):
        shutil.rmtree(report)
    os.makedirs(report)
    bundled = {}
    for name, files in manifest["tables"].items():
        if name in EXCLUDED_FROM_REPORT:
            continue
        copies = []
        for rel in files:
            src = os.path.join(outdir, rel)
            compress = re.fullmatch(r".+\.(tsv|txt)", rel) is not None
            dest_rel = rel + ".gz" if compress else rel
            dest = os.path.join(report, dest_rel)
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            if compress:
                with open(src, "rb") as fin, gzip.open(dest, "wb") as fout:
                    shutil.copyfileobj(fin, fout)
            else:
                shutil.copy2(src, dest)
            copies.append(dest_rel)
        bundled[name] = copies
    info = os.path.join(outdir, "run_info.txt")
    if os.path.isfile(info):
        shutil.copy2(info, os.path.join(report, "run_info.txt"))

    report_manifest = dict(manifest)
    report_manifest["tables"] = bundled
    report_manifest["bundle"] = {
        "excluded": [name for name in EXCLUDED_FROM_REPORT if name in manifest["tables"]]
    }
    write_json(os.path.join(report, "manifest.json"), report_manifest)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("results", help="results directory of the run (--outdir)")
    parser.add_argument("--run-info", help="run_info.txt (default: RESULTS/run_info.txt)")
    parser.add_argument("--force", action="store_true",
                        help="overwrite an existing manifest.json")
    parser.add_argument("--no-report", action="store_true", help="do not build report/")
    args = parser.parse_args(argv)

    outdir = os.path.abspath(args.results)
    if not os.path.isdir(outdir):
        sys.exit("Results directory not found: " + outdir)
    target = os.path.join(outdir, "manifest.json")
    if os.path.exists(target) and not args.force:
        sys.exit(target + " already exists (written by the pipeline?): use --force to replace it")

    run_info = args.run_info or os.path.join(outdir, "run_info.txt")
    if not os.path.isfile(run_info):
        sys.exit("run_info.txt not found: " + run_info)

    info = read_run_info(run_info)
    tables = find_tables(outdir)
    if not tables:
        sys.exit("No mag-pipeline table found in " + outdir)

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "source": {
            "written_by": "scripts/make_manifest.py",
            "date": datetime.date.today().isoformat(),
            "from": os.path.basename(run_info),
        },
        "pipeline": info["pipeline"],
        "run": info["run"],
        "databases": info["databases"],
        "params": info["params"],
        "software": {},
        "tables": tables,
    }
    write_json(target, manifest)
    n_files = sum(len(files) for files in tables.values())
    print("{}: {} tables, {} files".format(target, len(tables), n_files))

    if not args.no_report:
        report = build_report(outdir, manifest)
        size = sum(os.path.getsize(os.path.join(root, f))
                   for root, _, names in os.walk(report) for f in names)
        print("{}: {:.1f} MB".format(report, size / 1e6))


if __name__ == "__main__":
    main()
