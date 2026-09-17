process FASTP {
    tag "${sample}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/fastp:1.3.6--4df8d6c11b471bde'
    publishDir "${params.outdir}/01_qc", mode: 'copy', pattern: '*.{json,html}'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_R1.trimmed.fastq.gz"), path("${sample}_R2.trimmed.fastq.gz"), emit: reads
    path "${sample}.fastp.json", emit: json
    path "${sample}.fastp.html", emit: html

    script:
    """
    fastp \\
        --in1 ${r1} \\
        --in2 ${r2} \\
        --out1 ${sample}_R1.trimmed.fastq.gz \\
        --out2 ${sample}_R2.trimmed.fastq.gz \\
        --detect_adapter_for_pe \\
        --length_required ${params.min_read_length} \\
        --thread ${task.cpus} \\
        --json ${sample}.fastp.json \\
        --html ${sample}.fastp.html
    """

    stub:
    """
    echo -n | gzip > ${sample}_R1.trimmed.fastq.gz
    echo -n | gzip > ${sample}_R2.trimmed.fastq.gz
    touch ${sample}.fastp.json ${sample}.fastp.html
    """
}
