#!/bin/bash
#SBATCH --job-name=magpipe-gtdbtk
#SBATCH --partition=bigmem
#SBATCH --account=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.error

# Prueba del pipeline con clasificación taxonómica (GTDB-Tk).
# Lanzar DESDE la carpeta mag-pipeline con:  sbatch run_taxonomy_test.sh

set -euo pipefail

module load bioconda-tools/3
source activate "$HOME/envs/nextflow"

GTDB_DB=/softs/contrib/apps/anaconda/3/envs/gtdbtk-2.7.2/share/gtdbtk-2.7.2/db

nextflow run . \
    -profile test,pyrene \
    -resume \
    --refine \
    --dastool_score_threshold 0.05 \
    --gtdbtk_db "$GTDB_DB" \
    --outdir results_taxonomy

echo "Terminado. Tabla final:"
cat results_taxonomy/05_quality/mag_quality.tsv
