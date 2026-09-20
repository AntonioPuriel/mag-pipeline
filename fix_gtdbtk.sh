#!/usr/bin/env bash
# Corrige GTDB-Tk: la opcion --skip_ani_screen no existe en todas las versiones (v1.1.1)
# Ejecutar DENTRO de la carpeta mag-pipeline
set -euo pipefail

if [ ! -f modules/gtdbtk.nf ]; then echo "Ejecuta este script dentro de la carpeta mag-pipeline"; exit 1; fi

cat > modules/gtdbtk.nf << 'EOF_MAGPIPE'
process GTDB_TK {
    tag "${id} ${binner}"
    label 'process_high_memory'
    container 'quay.io/biocontainers/gtdbtk:2.4.0--pyhdfd78af_2'
    publishDir "${params.outdir}/06_taxonomy", mode: 'copy'

    input:
    tuple val(id), val(binner), path(bins_dir)
    path db

    output:
    tuple val(id), val(binner), path("${id}_${binner}.gtdbtk.tsv"), emit: summary

    script:
    """
    export GTDBTK_DATA_PATH=\$(readlink -f ${db})

    n_bins=\$(find -L ${bins_dir} -maxdepth 1 -name '*.fa' | wc -l)

    if [ "\$n_bins" -eq 0 ]; then
        echo "No bins in ${bins_dir}: writing an empty classification" >&2
        printf "user_genome\\tclassification\\tclosest_genome_ani\\n" > ${id}_${binner}.gtdbtk.tsv
        exit 0
    fi

    # --skip_ani_screen exists only in some GTDB-Tk versions: it was required in
    # 2.3-2.4 and dropped once skani replaced the Mash pre-screen
    ani_screen=""
    if gtdbtk classify_wf --help 2>&1 | grep -q -- "--skip_ani_screen"; then
        ani_screen="--skip_ani_screen"
    fi

    gtdbtk classify_wf \\
        --genome_dir ${bins_dir} \\
        --extension fa \\
        --out_dir gtdbtk_out \\
        --cpus ${task.cpus} \\
        --pplacer_cpus ${params.gtdbtk_pplacer_cpus} \\
        \$ani_screen \\
        ${params.gtdbtk_args}

    gtdbtk_summary.sh gtdbtk_out > ${id}_${binner}.gtdbtk.tsv
    """

    stub:
    """
    printf "user_genome\\tclassification\\tclosest_genome_ani\\n${id}.${binner}.1\\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\\t98.5\\n" > ${id}_${binner}.gtdbtk.tsv
    """
}
EOF_MAGPIPE

cat > nextflow.config << 'EOF_MAGPIPE'
manifest {
    name            = 'AntonioPuriel/mag-pipeline'
    author          = 'Antonio Puriel Hernández'
    description     = 'MAG recovery from paired-end shotgun metagenomes'
    nextflowVersion = '>=24.04.0'
    version         = '1.1.1'
}

params {
    input                   = null
    outdir                  = 'results'
    assembly_mode           = 'coassembly'         // 'coassembly' or 'per_sample'
    min_read_length         = 50
    min_contig_length       = 1500
    binners                 = 'concoct,metabat2'   // any of: concoct, metabat2
    refine                  = false                // combine bins with DAS_Tool
    dastool_score_threshold = 0.5
    checkm2_db              = null                 // path to uniref100.KO.1.dmnd
    skip_checkm2            = false
    gtdbtk_db               = null                 // GTDB-Tk reference data directory
    gtdbtk_pplacer_cpus     = 1                    // more threads multiply pplacer memory
    gtdbtk_args             = ''                   // extra flags for gtdbtk classify_wf
}

includeConfig 'conf/base.config'

profiles {
    test        { includeConfig 'conf/test.config' }
    docker      { docker.enabled = true }
    singularity {
        singularity.enabled    = true
        singularity.autoMounts = true
    }
    slurm       { includeConfig 'conf/slurm.config' }
    pyrene      { includeConfig 'conf/pyrene.config' }
}
EOF_MAGPIPE

echo "Corregido. Vuelve a lanzar la prueba con -resume"
