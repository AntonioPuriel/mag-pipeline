#!/usr/bin/env bash
# Combine per-bin statistics with CheckM2 results and (optionally) GTDB-Tk taxonomy,
# and assign a quality category:
#   high:   completeness >= 90 and contamination < 5
#   medium: completeness >= 50 and contamination < 10
#   low:    everything else
# Input files are recognised by their name:
#   <assembly>_<binner>.checkm2.tsv   CheckM2 report
#   <assembly>_<binner>.gtdbtk.tsv    GTDB-Tk classification
# Usage: mag_quality.sh <bin_summary.tsv> <report.tsv> [...]
set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $(basename "$0") <bin_summary.tsv> <report.tsv> [...]" >&2
    exit 1
fi

bin_summary=$1
shift

printf "assembly\tbinner\tbin\tn_contigs\ttotal_length\tn50\tgc_pct\tcompleteness\tcontamination\tquality\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tclosest_ani\n"

awk -F '\t' -v OFS='\t' -v summary="$bin_summary" '
    function rank(classification, prefix,   i, n, parts, value) {
        n = split(classification, parts, ";")
        for (i = 1; i <= n; i++) {
            if (substr(parts[i], 1, length(prefix)) == prefix) {
                value = substr(parts[i], length(prefix) + 1)
                return (value == "" ? "NA" : value)
            }
        }
        return "NA"
    }

    # Header of every input file: work out assembly, binner and column positions
    FNR == 1 && FILENAME != summary {
        base = FILENAME
        sub(/.*\//, "", base)
        kind = (base ~ /\.gtdbtk\.tsv$/) ? "gtdbtk" : "checkm2"
        sub(/\.(checkm2|gtdbtk)\.tsv$/, "", base)
        n = split(base, parts, "_")
        binner = parts[n]
        assembly = substr(base, 1, length(base) - length(binner) - 1)

        c_name = c_comp = c_cont = c_class = c_ani = 0
        for (i = 1; i <= NF; i++) {
            if ($i == "Name" || $i == "user_genome") c_name  = i
            if ($i == "Completeness")                c_comp  = i
            if ($i == "Contamination")               c_cont  = i
            if ($i == "classification")              c_class = i
            if ($i == "closest_genome_ani" || $i == "fastani_ani") c_ani = i
        }
        next
    }

    FILENAME != summary {
        key = assembly SUBSEP binner SUBSEP $c_name
        if (kind == "checkm2") {
            comp[key] = $c_comp
            cont[key] = $c_cont
        } else {
            tax[key] = (c_class ? $c_class : "")
            ani[key] = (c_ani  ? $c_ani  : "NA")
        }
        next
    }

    # Bin summary: keep the bins assessed by CheckM2
    FNR == 1 { next }
    {
        key = $1 SUBSEP $2 SUBSEP $3
        if (!(key in comp)) next

        c = comp[key] + 0
        k = cont[key] + 0
        quality = (c >= 90 && k < 5) ? "high" : (c >= 50 && k < 10) ? "medium" : "low"

        classification = (key in tax) ? tax[key] : ""
        closest = (key in ani) ? ani[key] : "NA"

        print $1, $2, $3, $4, $5, $6, $7, comp[key], cont[key], quality,
              rank(classification, "d__"), rank(classification, "p__"),
              rank(classification, "c__"), rank(classification, "o__"),
              rank(classification, "f__"), rank(classification, "g__"),
              rank(classification, "s__"), closest
    }
' "$@" "$bin_summary" | sort -t $'\t' -k1,1 -k2,2 -k8,8gr
