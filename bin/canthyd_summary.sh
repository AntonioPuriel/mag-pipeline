#!/usr/bin/env bash
# Turn the hmmsearch table of CANT-HYD hits into one row per protein.
#
# hmmsearch runs with --cut_tc, so every hit is already above the trusted
# cutoff of its model. A protein that passes more than one model (for example
# NdoB and MAH_alpha, both dioxygenase alpha subunits) keeps only the model with
# the highest score.
#
# Protein names carry their genome as <genome>|<gene> (see PRODIGAL).
#
# Usage: canthyd_summary.sh <hits.tblout> <CANT-HYD.hmm>
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $(basename "$0") <hits.tblout> <CANT-HYD.hmm>" >&2
    exit 1
fi

tblout=$1
hmm=$2

printf "genome\tprotein\tdatabase\tmarker\tdescription\tscore\tevalue\n"

awk -v OFS='\t' '
    # Model descriptions from the HMM file
    FNR == NR {
        if ($1 == "NAME") name = $2
        if ($1 == "DESC") { d = $0; sub(/^DESC +/, "", d); desc[name] = d }
        next
    }
    /^#/ { next }
    {
        protein = $1; model = $3; evalue = $5; score = $6
        if (!(protein in best) || score + 0 > best[protein] + 0) {
            best[protein] = score; hit[protein] = model; ev[protein] = evalue
        }
    }
    END {
        for (protein in best) {
            genome = protein; sub(/\|.*/, "", genome)
            model = hit[protein]
            print genome, protein, "CANT-HYD", model, (model in desc ? desc[model] : "NA"), best[protein], ev[protein]
        }
    }
' "$hmm" "$tblout" | sort -t $'\t' -k1,1 -k4,4 -k2,2
