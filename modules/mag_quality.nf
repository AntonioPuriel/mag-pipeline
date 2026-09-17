process MAG_QUALITY {
    label 'process_low'
    publishDir "${params.outdir}/05_quality", mode: 'copy'

    input:
    path bin_summary
    path checkm2_reports

    output:
    path "mag_quality.tsv", emit: tsv

    script:
    """
    mag_quality.sh ${bin_summary} ${checkm2_reports} > mag_quality.tsv
    """

    stub:
    """
    touch mag_quality.tsv
    """
}
