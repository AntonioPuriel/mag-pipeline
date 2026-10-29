#!/usr/bin/env bash
# Count 5S, 16S and 23S rRNA genes per genome from barrnap GFF files.
#
# Each genome is searched twice, with the bacterial and the archaeal models, and
# both results are reported: the domain is chosen later, in mimag_quality.py.
# A gene that barrnap marks as partial (shorter than the full model) still
# counts; the number of partial genes is reported in its own column.
#
# Usage: rrna_summary.sh <genomes.txt> <gff_dir>
#   genomes.txt  one FASTA path per line (genome name = file name without .fa)
#   gff_dir      <genome>.bac.gff and <genome>.arc.gff from barrnap
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $(basename "$0") <genomes.txt> <gff_dir>" >&2
    exit 1
fi

genomes=$1
gff_dir=$2

printf "genome\tmodel\trrna_5S\trrna_16S\trrna_23S\trrna_partial\n"

while read -r fasta; do
    [ -n "$fasta" ] || continue
    genome=$(basename "$fasta" .fa)
    for model in bac arc; do
        gff="$gff_dir/$genome.$model.gff"
        [ -e "$gff" ] || gff=/dev/null
        awk -F '\t' -v OFS='\t' -v genome="$genome" -v model="$model" '
            !/^#/ && NF >= 9 {
                name = ""; product = ""
                n = split($9, attributes, ";")
                for (i = 1; i <= n; i++) {
                    if (attributes[i] ~ /^Name=/)    name    = substr(attributes[i], 6)
                    if (attributes[i] ~ /^product=/) product = substr(attributes[i], 9)
                }
                sub(/_rRNA$/, "", name)
                if (name == "5S" || name == "16S" || name == "23S") {
                    count[name]++
                    if (tolower(product) ~ /partial/) partial++
                }
            }
            END { print genome, model, count["5S"] + 0, count["16S"] + 0, count["23S"] + 0, partial + 0 }
        ' "$gff"
    done
done < <(sort "$genomes")
