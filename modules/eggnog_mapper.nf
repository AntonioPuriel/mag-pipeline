process EGGNOG_MAPPER {
    tag "${id} ${binner}"
    label 'process_high'
    container 'quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0'
    publishDir "${params.outdir}/07_function", mode: 'copy', pattern: '*.emapper.annotations'

    input:
    tuple val(id), val(binner), path(proteins)
    path db

    output:
    tuple val(id), val(binner), path("${id}_${binner}.emapper.annotations"), emit: annotations

    script:
    def dbmem = params.eggnog_dbmem ? '--dbmem' : ''
    """
    emapper.py \\
        -i ${proteins} \\
        --itype proteins \\
        -m diamond \\
        --data_dir ${db} \\
        --output ${id}_${binner} \\
        --output_dir . \\
        --cpu ${task.cpus} \\
        --temp_dir . \\
        --override \\
        ${dbmem} \\
        ${params.eggnog_args}
    """

    stub:
    """
    printf "## emapper stub\\n#query\\tseed_ortholog\\tevalue\\tscore\\teggNOG_OGs\\tmax_annot_lvl\\tCOG_category\\tDescription\\tPreferred_name\\tGOs\\tEC\\tKEGG_ko\\tKEGG_Pathway\\n${id}.${binner}.1|k141_1_1\\t123.abc\\t1e-50\\t200\\tCOG0001\\tBacteria\\tC\\talkane 1-monooxygenase\\talkB\\t-\\t1.14.15.3\\tko:K00496\\tko00071\\n" > ${id}_${binner}.emapper.annotations
    """
}
