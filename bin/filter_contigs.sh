#!/usr/bin/env bash
# Filter assembled contigs by minimum length and report assembly statistics.
# Usage: filter_contigs.sh <contigs.fa[.gz]> <min_length> <output.fa.gz> <stats.tsv> <assembly_id>
set -euo pipefail

if [ "$#" -ne 5 ]; then
    echo "Usage: $(basename "$0") <contigs.fa[.gz]> <min_length> <output.fa.gz> <stats.tsv> <assembly_id>" >&2
    exit 1
fi

input=$1
min_len=$2
output=$3
stats=$4
id=$5

lengths=$(mktemp)
trap 'rm -f "$lengths"' EXIT

# Filter contigs and record the length of every contig (all / kept)
gzip -dcf "$input" | awk -v min="$min_len" -v lens="$lengths" '
    function flush() {
        if (name == "") return
        len = length(seq)
        print "all\t" len > lens
        if (len >= min) {
            print "kept\t" len > lens
            print name
            for (i = 1; i <= len; i += 80) print substr(seq, i, 80)
        }
    }
    /^>/ { flush(); name = $0; seq = ""; next }
    { gsub(/[ \t\r]/, ""); seq = seq $0 }
    END { flush() }
' | gzip > "$output"

# Number of contigs, total length, N50 and longest contig for one set
summarise() {
    local set=$1 label=$2
    awk -v set="$set" '$1 == set { print $2 }' "$lengths" | sort -nr | \
    awk -v id="$id" -v label="$label" '
        { len[NR] = $1; total += $1 }
        END {
            n50 = 0; cum = 0
            for (i = 1; i <= NR; i++) {
                cum += len[i]
                if (cum >= total / 2) { n50 = len[i]; break }
            }
            printf "%s\t%s\t%d\t%d\t%d\t%d\n", id, label, NR, total, n50, (NR > 0 ? len[1] : 0)
        }'
}

{
    printf "assembly\tcontig_set\tn_contigs\ttotal_length\tn50\tmax_length\n"
    summarise all  "all"
    summarise kept ">=${min_len}bp"
} > "$stats"
