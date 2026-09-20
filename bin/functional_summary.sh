#!/usr/bin/env bash
# Summarise eggNOG-mapper annotations per MAG.
#   ko:    one row per MAG and KEGG orthologue, with the number of genes
#   stats: one row per MAG with genes, annotated genes and distinct KOs
# Annotation files must be named <assembly>_<binner>.emapper.annotations and the
# query names must be <bin>|<gene>, as written by the pipeline.
# eggNOG-mapper only writes the genes it could annotate, so the total number of
# genes comes from the Prodigal tables (<assembly>_<binner>.gene_counts.tsv),
# which must sit next to the annotation files for the stats mode.
# Usage: functional_summary.sh <ko|stats> <file.emapper.annotations> [...]
set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $(basename "$0") <ko|stats> <file.emapper.annotations> [...]" >&2
    exit 1
fi

mode=$1
shift

case "$mode" in
    ko)    printf "assembly\tbinner\tbin\tko\tn_genes\n" ;;
    stats) printf "assembly\tbinner\tbin\tn_genes\tn_annotated\tn_kos\n" ;;
    *)     echo "Unknown mode: $mode" >&2; exit 1 ;;
esac

for file in "$@"; do
    base=$(basename "$file" .emapper.annotations)
    binner=${base##*_}
    assembly=${base%_*}

    counts="$(dirname "$file")/${base}.gene_counts.tsv"
    [ -f "$counts" ] || counts=/dev/null

    awk -F '\t' -v OFS='\t' -v assembly="$assembly" -v binner="$binner" -v mode="$mode" -v counts="$counts" '
        # Genes called by Prodigal, read before the annotations
        FILENAME == counts {
            if (FNR == 1) next
            called[$1] = $2
            next
        }

        # Header line of the annotation table
        /^#query/ {
            for (i = 1; i <= NF; i++) {
                if ($i == "#query")   c_query = i
                if ($i == "KEGG_ko")  c_ko = i
            }
            next
        }
        /^#/ { next }
        NF < 2 { next }

        {
            split($c_query, parts, "|")
            bin = parts[1]
            n_annot[bin]++

            if (c_ko && $c_ko != "-" && $c_ko != "") {
                n = split($c_ko, kos, ",")
                for (i = 1; i <= n; i++) {
                    ko = kos[i]
                    sub(/^ko:/, "", ko)
                    if (ko == "") continue
                    key = bin SUBSEP ko
                    if (!(key in ko_count)) ko_list[bin] = ko_list[bin] " " ko
                    ko_count[key]++
                }
            }
        }

        END {
            # Every bin with called genes, plus any bin seen only in the annotations
            for (bin in called) seen_bin[bin] = 1
            for (bin in n_annot) seen_bin[bin] = 1

            for (bin in seen_bin) {
                if (mode == "stats") {
                    n_kos = split(ko_list[bin], seen, " ")
                    distinct = 0
                    for (i = 1; i <= n_kos; i++) if (seen[i] != "") distinct++
                    n_called = (bin in called) ? called[bin] : n_annot[bin] + 0
                    print assembly, binner, bin, n_called, n_annot[bin] + 0, distinct
                } else {
                    n_kos = split(ko_list[bin], seen, " ")
                    for (i = 1; i <= n_kos; i++) {
                        ko = seen[i]
                        if (ko == "") continue
                        print assembly, binner, bin, ko, ko_count[bin SUBSEP ko]
                    }
                }
            }
        }
    ' "$counts" "$file"
done | sort -t $'\t' -k1,1 -k2,2 -k3,3 -k4,4
