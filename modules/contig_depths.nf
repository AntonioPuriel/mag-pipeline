process CONTIG_DEPTHS {
    tag "${id}"
    label 'process_low'
    container 'quay.io/biocontainers/metabat2:2.17--hd498684_0'
    publishDir "${params.outdir}/03_mapping", mode: 'copy'

    input:
    tuple val(id), path(bams), path(bais)

    output:
    tuple val(id), path("${id}.depth.txt"), emit: depth

    script:
    """
    jgi_summarize_bam_contig_depths \\
        --outputDepth ${id}.depth.txt \\
        ${bams}
    """

    stub:
    """
    touch ${id}.depth.txt
    """
}
