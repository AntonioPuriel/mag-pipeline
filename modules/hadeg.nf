process HADEG {
    tag "${id} ${binner}"
    label 'process_medium'
    // Check the exact biocontainers tag at https://quay.io/repository/biocontainers/diamond?tab=tags
    container 'quay.io/biocontainers/diamond:2.1.10--h5ca1c30_3'
    publishDir "${params.outdir}/07_function/hydrocarbons", mode: 'copy'

    input:
    tuple val(id), val(binner), path(proteins)
    path hadeg_dir

    output:
    tuple val(id), val(binner), path("${id}_${binner}.hadeg.tsv"), emit: hits

    script:
    """
    # The HADEG repository: the newest protein database and the pathway table
    faa=\$(ls ${hadeg_dir}/HADEG_protein_database_*.faa | sort | tail -n 1)
    pathways=${hadeg_dir}/Scripts/7_All_pathways.csv
    echo "HADEG database: \$(basename \$faa)" >&2

    # Aerobic hydrocarbon degradation, plastic degradation and biosurfactant
    # genes (HADEG). The identity threshold follows the HADEG authors (50 %);
    # both proteins must be covered over most of their length, so a shared
    # domain alone is not a hit.
    diamond makedb --in \$faa --db hadeg --threads ${task.cpus} --quiet

    diamond blastp \\
        --query ${proteins} \\
        --db hadeg \\
        --out hadeg_hits.tsv \\
        --id ${params.hadeg_identity} \\
        --query-cover ${params.hadeg_coverage} \\
        --subject-cover ${params.hadeg_coverage} \\
        --max-target-seqs 1 \\
        --sensitive \\
        --threads ${task.cpus} \\
        --quiet \\
        --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore

    hadeg_summary.sh hadeg_hits.tsv \$pathways > ${id}_${binner}.hadeg.tsv
    """

    stub:
    """
    printf "genome\\tprotein\\tdatabase\\tmarker\\tcompound\\tpathway\\tsubpathway\\thadeg_id\\tidentity\\tbitscore\\tevalue\\n" > ${id}_${binner}.hadeg.tsv
    printf "${id}.${binner}.1\\t${id}.${binner}.1|k141_1_1\\tHADEG\\tAeAE_alkB\\tA_Alkanes\\tA_Terminal/biterminal_oxidation\\tA_Terminal/biterminal_oxidation\\tQ9I0R2\\t85.0\\t700\\t1e-200\\n" >> ${id}_${binner}.hadeg.tsv
    """
}
