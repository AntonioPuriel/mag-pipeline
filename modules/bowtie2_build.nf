process BOWTIE2_BUILD {
    tag "${id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6'

    input:
    tuple val(id), path(contigs)

    output:
    tuple val(id), path("bowtie2_index"), emit: index

    script:
    """
    mkdir -p bowtie2_index
    bowtie2-build \\
        --threads ${task.cpus} \\
        ${contigs} \\
        bowtie2_index/${id}
    """

    stub:
    """
    mkdir -p bowtie2_index
    touch bowtie2_index/${id}.1.bt2
    """
}
