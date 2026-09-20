process FUNCTIONAL_SUMMARY {
    label 'process_low'
    publishDir "${params.outdir}/07_function", mode: 'copy'

    input:
    path annotations
    path gene_counts

    output:
    path "ko_per_mag.tsv",       emit: ko
    path "annotation_stats.tsv", emit: stats

    script:
    """
    functional_summary.sh ko ${annotations} > ko_per_mag.tsv
    functional_summary.sh stats ${annotations} > annotation_stats.tsv
    """

    stub:
    """
    touch ko_per_mag.tsv annotation_stats.tsv
    """
}
