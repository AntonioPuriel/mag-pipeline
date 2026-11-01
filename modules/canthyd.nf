process CANTHYD {
    tag "${id} ${binner}"
    label 'process_medium'
    // Check the exact biocontainers tag at https://quay.io/repository/biocontainers/hmmer?tab=tags
    container 'quay.io/biocontainers/hmmer:3.4--h7d74f8d_5'
    publishDir "${params.outdir}/07_function/hydrocarbons", mode: 'copy'

    input:
    tuple val(id), val(binner), path(proteins)
    path hmm

    output:
    tuple val(id), val(binner), path("${id}_${binner}.canthyd.tsv"), emit: hits

    script:
    """
    # Hydrocarbon degradation marker genes (CANT-HYD), aerobic and anaerobic.
    # --cut_tc keeps only hits above the trusted cutoff of each model.
    hmmsearch \\
        --cut_tc \\
        --cpu ${task.cpus} \\
        --noali \\
        --tblout canthyd.tblout \\
        ${hmm} \\
        ${proteins} > /dev/null

    canthyd_summary.sh canthyd.tblout ${hmm} > ${id}_${binner}.canthyd.tsv
    """

    stub:
    """
    printf "genome\\tprotein\\tdatabase\\tmarker\\tdescription\\tscore\\tevalue\\n" > ${id}_${binner}.canthyd.tsv
    printf "${id}.${binner}.1\\t${id}.${binner}.1|k141_1_1\\tCANT-HYD\\tAlkB\\talkane hydrolase\\t600.0\\t1e-180\\n" >> ${id}_${binner}.canthyd.tsv
    """
}
