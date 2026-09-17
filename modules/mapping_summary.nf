process MAPPING_SUMMARY {
    label 'process_low'
    publishDir "${params.outdir}/03_mapping", mode: 'copy'

    input:
    path logs

    output:
    path "mapping_summary.tsv", emit: tsv

    script:
    """
    mapping_summary.sh ${logs} > mapping_summary.tsv
    """

    stub:
    """
    touch mapping_summary.tsv
    """
}
