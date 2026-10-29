process GUNC {
    tag "${id} ${binner}"
    label 'process_medium'
    // Check the exact biocontainers tag at https://quay.io/repository/biocontainers/gunc?tab=tags
    container 'quay.io/biocontainers/gunc:1.0.6--pyhdfd78af_0'
    publishDir "${params.outdir}/05_quality/gunc", mode: 'copy'

    input:
    tuple val(id), val(binner), path(bins_dir)
    path db

    output:
    tuple val(id), val(binner), path("${id}_${binner}.gunc.tsv"), emit: report

    script:
    """
    n_bins=\$(find -L ${bins_dir} -maxdepth 1 -name '*.fa' | wc -l)

    if [ "\$n_bins" -eq 0 ]; then
        echo "No bins in ${bins_dir}: writing an empty report" >&2
        printf "genome\\tpass.GUNC\\n" > ${id}_${binner}.gunc.tsv
        exit 0
    fi

    mkdir -p gunc_out gunc_tmp

    # GUNC flags chimeric genomes: contigs whose genes point to different
    # lineages. It only reports, it does not change any bin.
    gunc run \\
        --input_dir ${bins_dir} \\
        --file_suffix .fa \\
        --db_file ${db} \\
        --threads ${task.cpus} \\
        --temp_dir gunc_tmp \\
        --out_dir gunc_out

    # One row per genome, at the taxonomic level with the highest clade
    # separation score; pass.GUNC is False for chimeric genomes
    cp gunc_out/GUNC.*.maxCSS_level.tsv ${id}_${binner}.gunc.tsv
    rm -rf gunc_tmp gunc_out/diamond_output
    """

    stub:
    """
    printf "genome\\tpass.GUNC\\n" > ${id}_${binner}.gunc.tsv
    for f in \$(find -L ${bins_dir} -maxdepth 1 -name '*.fa'); do
        printf "\$(basename \$f .fa)\\tTrue\\n" >> ${id}_${binner}.gunc.tsv
    done
    """
}
