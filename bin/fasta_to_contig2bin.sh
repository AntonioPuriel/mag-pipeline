#!/usr/bin/env bash
# Build a contig-to-bin table (contig <tab> bin) from a directory of bin FASTA files.
# Usage: fasta_to_contig2bin.sh <bin_dir> [extension, default: fa]
set -euo pipefail

if [ "$#" -lt 1 ]; then
    echo "Usage: $(basename "$0") <bin_dir> [extension]" >&2
    exit 1
fi

bin_dir=$1
ext=${2:-fa}

for fa in "$bin_dir"/*."$ext"; do
    [ -e "$fa" ] || continue
    bin=$(basename "$fa" ."$ext")
    awk -v bin="$bin" '/^>/ { sub(/^>/, ""); sub(/[ \t].*/, ""); print $0 "\t" bin }' "$fa"
done
