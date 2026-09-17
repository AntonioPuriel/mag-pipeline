#!/usr/bin/env bash
# Combine per-bin statistics with CheckM2 results and assign a quality category.
#   high:   completeness >= 90 and contamination < 5
#   medium: completeness >= 50 and contamination < 10
#   low:    everything else
# CheckM2 reports must be named <assembly>_<binner>.checkm2.tsv
# Usage: mag_quality.sh <bin_summary.tsv> <assembly>_<binner>.checkm2.tsv [...]
set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $(basename "$0") <bin_summary.tsv> <assembly>_<binner>.checkm2.tsv [...]" >&2
    exit 1
fi

bin_summary=$1
shift

printf "assembly\tbinner\tbin\tn_contigs\ttotal_length\tn50\tgc_pct\tcompleteness\tcontamination\tquality\n"

awk -F '\t' -v OFS='\t' -v summary="$bin_summary" '
    # CheckM2 reports: store completeness and contamination per assembly/binner/bin
    FILENAME != summary {
        if (FNR == 1) {
            base = FILENAME
            sub(/.*\//, "", base)
            sub(/\.checkm2\.tsv$/, "", base)
            n = split(base, parts, "_")
            binner = parts[n]
            assembly = substr(base, 1, length(base) - length(binner) - 1)
            for (i = 1; i <= NF; i++) {
                if ($i == "Name")          c_name = i
                if ($i == "Completeness")  c_comp = i
                if ($i == "Contamination") c_cont = i
            }
            next
        }
        key = assembly SUBSEP binner SUBSEP $c_name
        comp[key] = $c_comp
        cont[key] = $c_cont
        next
    }
    # Bin summary: keep bins assessed by CheckM2
    FNR == 1 { next }
    {
        key = $1 SUBSEP $2 SUBSEP $3
        if (!(key in comp)) next
        c = comp[key] + 0
        k = cont[key] + 0
        quality = (c >= 90 && k < 5) ? "high" : (c >= 50 && k < 10) ? "medium" : "low"
        print $1, $2, $3, $4, $5, $6, $7, comp[key], cont[key], quality
    }
' "$@" "$bin_summary" | sort -t $'\t' -k1,1 -k2,2 -k8,8gr
