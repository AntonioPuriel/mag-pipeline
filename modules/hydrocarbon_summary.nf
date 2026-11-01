process HYDROCARBON_SUMMARY {
    label 'process_low'
    publishDir "${params.outdir}/07_function", mode: 'copy'

    input:
    path canthyd, stageAs: 'canthyd/*'
    path hadeg,   stageAs: 'hadeg/*'

    output:
    path "hydrocarbon_per_mag.tsv", emit: tsv

    // No stub block: built from the stub outputs of CANTHYD and HADEG
    script:
    """
    hydrocarbon_summary.sh canthyd/* hadeg/* > hydrocarbon_per_mag.tsv
    """
}
