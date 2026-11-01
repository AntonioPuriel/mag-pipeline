#!/usr/bin/env bash
# Annotate DIAMOND hits against HADEG with the gene, pathway and compound of
# each HADEG protein, one row per query protein (its best hit).
#
# HADEG sequence names start with the UniProt accession (<accession>|...), which
# is the Protein_ID column of the HADEG pathway table (Scripts/7_All_pathways.csv).
# Proteins absent from that table are reported with pathway NA.
#
# Protein names carry their genome as <genome>|<gene> (see PRODIGAL).
#
# Usage: hadeg_summary.sh <diamond_hits.tsv> <7_All_pathways.csv>
#   diamond_hits.tsv  --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $(basename "$0") <diamond_hits.tsv> <7_All_pathways.csv>" >&2
    exit 1
fi

hits=$1
pathways=$2

printf "genome\tprotein\tdatabase\tmarker\tcompound\tpathway\tsubpathway\thadeg_id\tidentity\tbitscore\tevalue\n"

awk -v OFS='\t' '
    # Pathway table (CSV with Windows line endings, no quoted fields)
    FNR == NR {
        sub(/\r$/, "")
        if (FNR == 1) {
            n = split($0, h, ",")
            for (i = 1; i <= n; i++) col[h[i]] = i
            next
        }
        split($0, f, ",")
        id = f[col["Protein_ID"]]
        gene[id] = f[col["Gene"]]; compound[id] = f[col["Compound"]]
        pathway[id] = f[col["Pathway"]]; subpathway[id] = f[col["Subpathway"]]
        next
    }
    {
        protein = $1; subject = $2; bits = $8 + 0
        if (protein in best && bits <= best[protein]) next
        best[protein] = bits
        id = subject; sub(/\|.*/, "", id)
        row[protein] = id SUBSEP $3 SUBSEP $7
    }
    END {
        for (protein in best) {
            split(row[protein], r, SUBSEP)
            id = r[1]
            genome = protein; sub(/\|.*/, "", genome)
            known = (id in gene)
            print genome, protein, "HADEG",
                  (known ? gene[id] : "NA"), (known ? compound[id] : "NA"),
                  (known ? pathway[id] : "NA"), (known ? subpathway[id] : "NA"),
                  id, r[2], best[protein], r[3]
        }
    }
' FS=',' "$pathways" FS='\t' "$hits" | sort -t $'\t' -k1,1 -k5,5 -k4,4 -k2,2
