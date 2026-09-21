process BOWTIE2_ALIGN {
    tag "${sample} -> ${id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6'
    publishDir "${params.outdir}/03_mapping/logs", mode: 'copy', pattern: '*.log'

    input:
    tuple val(id), path(index), val(sample), path(r1), path(r2)

    output:
    tuple val(id), val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: bam
    path "${id}__${sample}.bowtie2.log",                                                             emit: log

    script:
    """
    set -o pipefail

    bowtie2 \\
        -x ${index}/${id} \\
        -1 ${r1} \\
        -2 ${r2} \\
        --threads ${task.cpus} \\
        2> ${id}__${sample}.bowtie2.log \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -o ${sample}.sorted.bam \\
        -

    samtools index ${sample}.sorted.bam
    """

    stub:
    """
    touch ${sample}.sorted.bam ${sample}.sorted.bam.bai ${id}__${sample}.bowtie2.log
    """
}
