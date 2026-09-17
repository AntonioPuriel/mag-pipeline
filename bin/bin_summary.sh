#!/usr/bin/env bash
# Summarise genome bins: number of contigs, size, N50 and GC content per bin.
# Bin directories must be named <assembly>_<binner>_bins.
# Usage: bin_summary.sh <assembly>_<binner>_bins [...]
set -euo pipefail

if [ "$#" -eq 0 ]; then
    echo "Usage: $(basename "$0") <assembly>_<binner>_bins [...]" >&2
    exit 1
fi

printf "assembly\tbinner\tbin\tn_contigs\ttotal_length\tn50\tgc_pct\n"

for dir in "$@"; do
    dir_name=$(basename "$dir")
    dir_name=${dir_name%_bins}
    assembly=${dir_name%_*}
    binner=${dir_name##*_}

    for fa in "$dir"/*.fa; do
        [ -e "$fa" ] || continue
        bin=$(basename "$fa" .fa)

        # GC content over all A/C/G/T bases
        gc=$(awk '
            !/^>/ {
                seq = toupper($0)
                gc += gsub(/[GC]/, "", seq)
                at += gsub(/[AT]/, "", seq)
            }
            END { printf "%.2f", (gc + at > 0 ? 100 * gc / (gc + at) : 0) }
        ' "$fa")

        # Contig lengths, then number of contigs, total length and N50
        awk '
            /^>/ { if (len > 0) print len; len = 0; next }
            { gsub(/[ \t\r]/, ""); len += length($0) }
            END { if (len > 0) print len }
        ' "$fa" | sort -nr | \
        awk -v a="$assembly" -v b="$binner" -v n="$bin" -v gc="$gc" '
            { len[NR] = $1; total += $1 }
            END {
                cum = 0; n50 = 0
                for (i = 1; i <= NR; i++) {
                    cum += len[i]
                    if (cum >= total / 2) { n50 = len[i]; break }
                }
                printf "%s\t%s\t%s\t%d\t%d\t%d\t%s\n", a, b, n, NR, total, n50, gc
            }'
    done
done
