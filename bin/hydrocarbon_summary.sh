#!/usr/bin/env bash
# Number of hydrocarbon degradation genes per MAG and marker, from the tables
# of CANTHYD (*.canthyd.tsv) and HADEG (*.hadeg.tsv). One row per genome,
# database and marker:
#   CANT-HYD: category = model description (e.g. "alkane hydrolase")
#   HADEG:    category = compound (A_Alkanes, C_Aromatics, ...), pathway = HADEG pathway
# Files named NO_FILE are ignored.
#
# Usage: hydrocarbon_summary.sh <table.tsv> [...]
set -euo pipefail

printf "genome\tdatabase\tmarker\tcategory\tpathway\tn_genes\n"

files=()
for f in "$@"; do
    [ "$(basename "$f")" = "NO_FILE" ] && continue
    [ -s "$f" ] && files+=("$f")
done
[ "${#files[@]}" -gt 0 ] || exit 0

awk -F '\t' -v OFS='\t' '
    FNR == 1 {
        c_cat = c_path = 0
        for (i = 1; i <= NF; i++) {
            if ($i == "description" || $i == "compound") c_cat = i
            if ($i == "pathway") c_path = i
        }
        next
    }
    {
        # columns 1-4 are always genome, protein, database, marker
        key = $1 SUBSEP $3 SUBSEP $4 SUBSEP (c_cat ? $c_cat : "NA") SUBSEP (c_path ? $c_path : "NA")
        n[key]++
    }
    END {
        for (key in n) {
            split(key, k, SUBSEP)
            print k[1], k[2], k[3], k[4], k[5], n[key]
        }
    }
' "${files[@]}" | sort -t $'\t' -k1,1 -k2,2 -k3,3
