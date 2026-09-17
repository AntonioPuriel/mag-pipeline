process MEGAHIT {
    tag "${id}"
    label 'process_high'
    container 'community.wave.seqera.io/library/megahit_pigz:87a590163e594224'
    publishDir "${params.outdir}/02_assembly/raw", mode: 'copy', pattern: '*.{fa.gz,log}'

    input:
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), path("${id}.contigs.fa.gz"), emit: contigs
    path "${id}.megahit.log",                   emit: log

    script:
    def r1_list = [r1].flatten().join(',')
    def r2_list = [r2].flatten().join(',')
    """
    megahit \\
        -1 ${r1_list} \\
        -2 ${r2_list} \\
        --num-cpu-threads ${task.cpus} \\
        --memory ${task.memory.toBytes()} \\
        --out-dir megahit_out \\
        --out-prefix ${id}

    gzip megahit_out/${id}.contigs.fa
    mv megahit_out/${id}.contigs.fa.gz .
    mv megahit_out/${id}.log ${id}.megahit.log
    """

    stub:
    """
    printf ">k141_1\\nACGT\\n" | gzip > ${id}.contigs.fa.gz
    touch ${id}.megahit.log
    """
}
