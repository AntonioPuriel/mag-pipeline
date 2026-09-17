process FILTER_CONTIGS {
    tag "${id}"
    label 'process_low'
    publishDir "${params.outdir}/02_assembly", mode: 'copy'

    input:
    tuple val(id), path(contigs)

    output:
    tuple val(id), path("${id}.contigs.min${params.min_contig_length}.fa.gz"), emit: contigs
    path "${id}.contigs.stats.tsv",                                          emit: stats

    script:
    """
    filter_contigs.sh \\
        ${contigs} \\
        ${params.min_contig_length} \\
        ${id}.contigs.min${params.min_contig_length}.fa.gz \\
        ${id}.contigs.stats.tsv \\
        ${id}
    """

    stub:
    """
    printf ">k141_1\\nACGT\\n" | gzip > ${id}.contigs.min${params.min_contig_length}.fa.gz
    touch ${id}.contigs.stats.tsv
    """
}
