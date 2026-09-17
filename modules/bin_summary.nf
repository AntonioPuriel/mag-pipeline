process BIN_SUMMARY {
    label 'process_low'
    publishDir "${params.outdir}/04_binning", mode: 'copy'

    input:
    path bin_dirs

    output:
    path "bin_summary.tsv", emit: tsv

    script:
    """
    bin_summary.sh ${bin_dirs} > bin_summary.tsv
    """

    stub:
    """
    touch bin_summary.tsv
    """
}
