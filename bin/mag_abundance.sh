#!/usr/bin/env bash
# Length-weighted mean depth of every MAG in every sample.
#
# Each contig contributes to its bin in proportion to its length, so the value
# is the average coverage of the genome rather than the average over contigs:
#     mean_depth = sum(depth_i * length_i) / sum(length_i)
#
# Contig names repeat across assemblies (every MEGAHIT run has a k141_1), so
# contigs are matched within their assembly, taken from the depth file name.
# Usage: mag_abundance.sh <contig2bin.tsv> <assembly>.depth.txt [...]
set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $(basename "$0") <contig2bin.tsv> <depth.txt> [...]" >&2
    exit 1
fi

contig2bin=$1
shift

printf "assembly\tbinner\tbin\tsample\tmean_depth\tbin_length\n"

for depth in "$@"; do
    this_assembly=$(basename "$depth" .depth.txt)

    awk -F '\t' -v OFS='\t' -v map="$contig2bin" -v this_assembly="$this_assembly" '
        # contig -> bin, read first
        FILENAME == map {
            if (FNR == 1) next
            if ($1 != this_assembly) next
            binner_of[$4, ++n_bins_of[$4]] = $2
            bin_of[$4, n_bins_of[$4]]      = $3
            next
        }

        # Depth table: work out which columns hold the per-sample depths
        FNR == 1 {
            for (i = 4; i <= NF; i++) {
                if ($i ~ /[.-]var$/) continue
                sample = $i
                sub(/\.bam$/, "", sample)
                sub(/\.sorted$/, "", sample)
                col[++n_samples] = i
                name[n_samples] = sample
            }
            next
        }

        {
            contig = $1
            if (!(contig in n_bins_of)) next
            len = $2

            # a contig belongs to one bin per binner
            for (b = 1; b <= n_bins_of[contig]; b++) {
                key = this_assembly SUBSEP binner_of[contig, b] SUBSEP bin_of[contig, b]
                length_sum[key] += len
                for (s = 1; s <= n_samples; s++) {
                    weighted[key, s] += $(col[s]) * len
                }
            }
        }

        END {
            for (key in length_sum) {
                split(key, k, SUBSEP)
                for (s = 1; s <= n_samples; s++) {
                    printf "%s\t%s\t%s\t%s\t%.4f\t%d\n",
                        k[1], k[2], k[3], name[s],
                        weighted[key, s] / length_sum[key], length_sum[key]
                }
            }
        }
    ' "$contig2bin" "$depth"
done | sort -t $'\t' -k1,1 -k2,2 -k3,3 -k4,4
