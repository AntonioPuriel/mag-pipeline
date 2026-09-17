process CHECKM2 {
    tag "${id} ${binner}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/checkm2:1.1.0--60f287bc25d7a10d'
    publishDir "${params.outdir}/05_quality/checkm2", mode: 'copy'

    input:
    tuple val(id), val(binner), path(bins_dir)
    path db

    output:
    tuple val(id), val(binner), path("${id}_${binner}.checkm2.tsv"), emit: report

    script:
    """
    n_bins=\$(find -L ${bins_dir} -maxdepth 1 -name '*.fa' | wc -l)

    if [ "\$n_bins" -eq 0 ]; then
        echo "No bins in ${bins_dir}: writing an empty report" >&2
        printf "Name\\tCompleteness\\tContamination\\n" > ${id}_${binner}.checkm2.tsv
        exit 0
    fi

    checkm2 predict \\
        --input ${bins_dir} \\
        --extension fa \\
        --output-directory checkm2_out \\
        --database_path ${db} \\
        --threads ${task.cpus} \\
        --remove_intermediates

    cp checkm2_out/quality_report.tsv ${id}_${binner}.checkm2.tsv
    """

    stub:
    """
    printf "Name\\tCompleteness\\tContamination\\n${id}.${binner}.1\\t95.0\\t1.0\\n" > ${id}_${binner}.checkm2.tsv
    """
}
