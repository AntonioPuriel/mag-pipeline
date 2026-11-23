#!/usr/bin/env bash
# Test of scripts/make_manifest.py: manifest.json and report/ rebuilt from the
# results of a run, and the same list of tables as workflow.onComplete in
# nextflow.config.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root="$here/.."
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
res="$tmp/results"

# The two lists of tables must agree
python3 - "$root/nextflow.config" "$root/scripts/make_manifest.py" <<'EOF'
import re, sys
config = open(sys.argv[1]).read()
block = config[config.index("def known = ["):config.index("def tables = [:]")]
groovy = re.findall(r"(\w+)\s*:\s*\['([^']+)',\s*/([^/]+)/\]", block)
script = open(sys.argv[2]).read()
python = re.findall(r'\("(\w+)",\s*"([^"]+)",\s*r"([^"]+)"\)', script)
assert groovy, "no tables found in nextflow.config"
assert groovy == python, "nextflow.config and make_manifest.py list different tables:\n{}\n{}".format(groovy, python)
print("tables in sync:", len(groovy))
EOF

mkdir -p "$res"/{02_assembly,03_mapping,04_binning,05_quality/gunc,06_taxonomy/mag_tree,09_singlem/profiles}
printf "stage\tn_contigs\nfiltered\t10\n" > "$res/02_assembly/T0.contigs.stats.tsv"
printf "contigName\tcontigLen\n" > "$res/03_mapping/T0.depth.txt"
printf "sample\tread_pairs\ns1\t100\n" > "$res/03_mapping/mapping_summary.tsv"
printf "contig\tbin\n" > "$res/04_binning/contig2bin.tsv"
printf "bin\ts1\n" | gzip > "$res/04_binning/mag_abundance.tsv.gz"
printf "bin\tcompleteness\n" > "$res/05_quality/mag_quality.tsv"
printf "genome\tpass.GUNC\n" > "$res/05_quality/gunc/T0_dastool.gunc.tsv"
printf "genome\tpass.GUNC\n" > "$res/05_quality/gunc/T1_dastool.gunc.tsv"
printf "(a,b);\n" > "$res/06_taxonomy/mag_tree/bac120.tree"
printf "taxonomy\n" > "$res/09_singlem/profiles/s1.singlem.tsv"
printf "not a table\n" > "$res/05_quality/notes.txt"

cat > "$res/run_info.txt" <<'EOF'
mag-pipeline run information

pipeline         : AntonioPuriel/mag-pipeline 2.0.0
nextflow         : 26.04.6
started          : 2026-11-11T09:17:54+01:00
completed        : 2026-11-13T05:58:06+01:00
duration         : 1d 20h 40m 12s
success          : true
profile          : pyrene
container engine : none (conda or local tools)
work directory   : /work/x
launch directory : /work/y
user             : someone

command line:
  nextflow run /path/mag-pipeline -profile pyrene --workflow both

reference databases:
  CheckM2 : /db/uniref100.KO.1.dmnd (last modified 2021-03-23)
  GTDB-Tk : /db/gtdbtk (last modified unknown)
  CAT     : not used

parameters:
  assembly_mode           = group
  binners                 = concoct,metabat2
  cat_db                  = null
  catalog_chunk_size      = 1000000
  catalog_min_id          = 0.95
  refine                  = true
EOF

python3 "$root/scripts/make_manifest.py" "$res"

# A second run refuses to overwrite the manifest, unless --force
if python3 "$root/scripts/make_manifest.py" "$res" 2>/dev/null; then
    echo "FAIL: overwrote manifest.json without --force"; exit 1
fi
python3 "$root/scripts/make_manifest.py" "$res" --force > /dev/null

python3 - "$res" <<'EOF'
import gzip, json, os, sys
res = sys.argv[1]
m = json.load(open(os.path.join(res, "manifest.json")))
assert m["schema_version"] == "1.0"
assert m["source"]["written_by"] == "scripts/make_manifest.py"
assert m["pipeline"] == {"name": "AntonioPuriel/mag-pipeline", "version": "2.0.0",
                         "revision": None, "commit": None, "nextflow": "26.04.6"}, m["pipeline"]
assert m["run"]["success"] is True
assert m["run"]["command_line"] == "nextflow run /path/mag-pipeline -profile pyrene --workflow both"
assert m["databases"] == {"checkm2": {"path": "/db/uniref100.KO.1.dmnd", "modified": "2021-03-23"},
                          "gtdbtk": {"path": "/db/gtdbtk", "modified": None}}, m["databases"]
assert m["params"] == {"assembly_mode": "group", "binners": "concoct,metabat2", "cat_db": None,
                       "catalog_chunk_size": 1000000, "catalog_min_id": 0.95, "refine": True}, m["params"]
assert m["software"] == {}
expected = {
    "assembly_stats": ["02_assembly/T0.contigs.stats.tsv"],
    "mapping_summary": ["03_mapping/mapping_summary.tsv"],
    "depth": ["03_mapping/T0.depth.txt"],
    "contig2bin": ["04_binning/contig2bin.tsv"],
    "mag_abundance": ["04_binning/mag_abundance.tsv.gz"],
    "mag_quality": ["05_quality/mag_quality.tsv"],
    "gunc": ["05_quality/gunc/T0_dastool.gunc.tsv", "05_quality/gunc/T1_dastool.gunc.tsv"],
    "mag_tree": ["06_taxonomy/mag_tree/bac120.tree"],
    "singlem_per_sample": ["09_singlem/profiles/s1.singlem.tsv"],
}
assert m["tables"] == expected, m["tables"]

r = json.load(open(os.path.join(res, "report", "manifest.json")))
assert r["bundle"] == {"excluded": ["depth", "contig2bin"]}
assert "depth" not in r["tables"] and "contig2bin" not in r["tables"]
assert r["tables"]["mag_quality"] == ["05_quality/mag_quality.tsv.gz"]
assert r["tables"]["mag_abundance"] == ["04_binning/mag_abundance.tsv.gz"]  # already gzipped
assert r["tables"]["mag_tree"] == ["06_taxonomy/mag_tree/bac120.tree"]      # trees stay as they are
for files in r["tables"].values():
    for f in files:
        assert os.path.isfile(os.path.join(res, "report", f)), f
with gzip.open(os.path.join(res, "report", "05_quality/mag_quality.tsv.gz"), "rt") as h:
    assert h.read() == "bin\tcompleteness\n"
assert os.path.isfile(os.path.join(res, "report", "run_info.txt"))
print("make_manifest.py: OK")
EOF
