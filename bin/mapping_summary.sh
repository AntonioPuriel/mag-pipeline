#!/usr/bin/env bash
# Summarise Bowtie2 logs into one table: assembly, sample, read pairs and overall
# alignment rate. Logs are named <assembly>__<sample>.bowtie2.log, since one
# sample can be mapped against several assemblies.
# Usage: mapping_summary.sh <assembly>__<sample>.bowtie2.log [...]
set -euo pipefail

if [ "$#" -eq 0 ]; then
    echo "Usage: $(basename "$0") <sample.bowtie2.log> [...]" >&2
    exit 1
fi

printf "assembly\tsample\tread_pairs\toverall_alignment_rate_pct\n"

for log in "$@"; do
    name=$(basename "$log" .bowtie2.log)
    assembly=${name%%__*}
    sample=${name#*__}
    pairs=$(awk '/reads; of these:/ { print $1; exit }' "$log")
    rate=$(awk '/overall alignment rate/ { sub("%", "", $1); print $1; exit }' "$log")
    printf "%s\t%s\t%s\t%s\n" "$assembly" "$sample" "${pairs:-NA}" "${rate:-NA}"
done | sort
