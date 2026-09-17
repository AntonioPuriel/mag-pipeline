process METABAT2 {
    tag "${id}"
    label 'process_medium'
    container 'quay.io/biocontainers/metabat2:2.17--hd498684_0'
    publishDir "${params.outdir}/04_binning/metabat2", mode: 'copy'

    input:
    tuple val(id), path(contigs), path(depth)

    output:
    tuple val(id), val('metabat2'), path("${id}_metabat2_bins"), emit: bins

    script:
    // MetaBAT2 does not accept contigs shorter than 1500 bp
    def min_len = Math.max(1500, params.min_contig_length as int)
    """
    mkdir -p ${id}_metabat2_bins

    metabat2 \\
        --inFile ${contigs} \\
        --abdFile ${depth} \\
        --outFile ${id}_metabat2_bins/${id}.metabat2 \\
        --minContig ${min_len} \\
        --numThreads ${task.cpus} \\
        --seed 42
    """

    stub:
    """
    mkdir -p ${id}_metabat2_bins
    printf ">k141_1\\nACGT\\n" > ${id}_metabat2_bins/${id}.metabat2.1.fa
    """
}
