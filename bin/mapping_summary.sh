#!/usr/bin/env bash
# Summarise Bowtie2 logs into one table: sample, read pairs and overall alignment rate.
# Usage: mapping_summary.sh <sample1.bowtie2.log> [sample2.bowtie2.log ...]
set -euo pipefail

if [ "$#" -eq 0 ]; then
    echo "Usage: $(basename "$0") <sample.bowtie2.log> [...]" >&2
    exit 1
fi

printf "sample\tread_pairs\toverall_alignment_rate_pct\n"

for log in "$@"; do
    sample=$(basename "$log" .bowtie2.log)
    pairs=$(awk '/reads; of these:/ { print $1; exit }' "$log")
    rate=$(awk '/overall alignment rate/ { sub("%", "", $1); print $1; exit }' "$log")
    printf "%s\t%s\t%s\n" "$sample" "${pairs:-NA}" "${rate:-NA}"
done | sort
