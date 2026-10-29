#!/usr/bin/env bash
# Count tRNA genes per genome from tRNAscan-SE output (--brief format).
#
# MIMAG asks for "at least 18 tRNAs" in a high-quality draft, read as tRNAs for
# at least 18 of the 20 standard amino acids. trna_types counts those amino
# acids, trna_total every tRNA gene for them, and trna_missing lists the absent
# ones. Pseudogenes, selenocysteine, suppressor and undetermined tRNAs are left
# out; fMet, iMet and Ile2 count as Met and Ile.
#
# Each genome is searched with the bacterial (-B) and the archaeal (-A) models
# and both results are reported: the domain is chosen later, in mimag_quality.py.
#
# Usage: trna_summary.sh <genomes.txt> <trnascan_dir>
#   genomes.txt   one FASTA path per line (genome name = file name without .fa)
#   trnascan_dir  <genome>.bac.txt and <genome>.arc.txt from tRNAscan-SE --brief
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $(basename "$0") <genomes.txt> <trnascan_dir>" >&2
    exit 1
fi

genomes=$1
trnascan_dir=$2

printf "genome\tmodel\ttrna_types\ttrna_total\ttrna_missing\n"

while read -r fasta; do
    [ -n "$fasta" ] || continue
    genome=$(basename "$fasta" .fa)
    for model in bac arc; do
        txt="$trnascan_dir/$genome.$model.txt"
        [ -e "$txt" ] || txt=/dev/null
        awk -F '\t' -v OFS='\t' -v genome="$genome" -v model="$model" '
            BEGIN {
                n_aa = split("Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val", aa, " ")
                for (i = 1; i <= n_aa; i++) standard[aa[i]] = 1
            }
            {
                tag = $2; gsub(/ /, "", tag)
                if (NF < 9 || tag !~ /^[0-9]+$/) next     # not a tRNA row
                if (tolower($10) ~ /pseudo/) next
                type = $5; gsub(/ /, "", type)
                if (type == "fMet" || type == "iMet") type = "Met"
                sub(/[0-9]+$/, "", type)                   # Ile2 -> Ile
                if (type in standard) count[type]++
            }
            END {
                types = total = 0; missing = ""
                for (i = 1; i <= n_aa; i++) {
                    if (aa[i] in count) { types++; total += count[aa[i]] }
                    else missing = missing (missing == "" ? "" : ",") aa[i]
                }
                print genome, model, types, total, (missing == "" ? "-" : missing)
            }
        ' "$txt"
    done
done < <(sort "$genomes")
