#!/usr/bin/env bash
# Merge the bacterial and archaeal GTDB-Tk summaries into a single table.
# Usage: gtdbtk_summary.sh <gtdbtk_out_dir>
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $(basename "$0") <gtdbtk_out_dir>" >&2
    exit 1
fi

out_dir=$1
shopt -s nullglob

summaries=("$out_dir"/*.summary.tsv "$out_dir"/classify/*.summary.tsv)

if [ "${#summaries[@]}" -eq 0 ]; then
    echo "No GTDB-Tk summary found in $out_dir" >&2
    printf "user_genome\tclassification\tclosest_genome_ani\n"
    exit 0
fi

head -n 1 "${summaries[0]}"
for f in "${summaries[@]}"; do
    tail -n +2 "$f"
done
