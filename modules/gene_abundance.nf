/*
 * Depth and TPM of every catalog gene in every sample (--workflow contigs).
 * Each gene inherits the depth of its contig from the binning depth tables:
 * no new mapping. See bin/gene_abundance.py for how clusters and assemblies
 * are combined.
 */

process GENE_ABUNDANCE {
    label 'process_medium'
    // Check the exact tag before using containers: any image with Python 3, numpy and pandas
    container 'quay.io/biocontainers/pandas:2.2.1'
    publishDir "${params.outdir}/08_contigs/abundance", mode: 'copy'

    input:
    path clusters
    path depths

    output:
    path "gene_depth.tsv.gz",          emit: depth
    path "gene_tpm.tsv.gz",            emit: tpm
    path "gene_abundance_stats.tsv",   emit: stats

    script:
    """
    gene_abundance.py \\
        --clusters ${clusters} \\
        --depths ${depths} \\
        --out-depth gene_depth.tsv.gz \\
        --out-tpm gene_tpm.tsv.gz \\
        --out-stats gene_abundance_stats.tsv
    """

    stub:
    """
    printf "gene\\ts1\\na|k141_1_1\\t10\\n" | gzip > gene_depth.tsv.gz
    printf "gene\\ts1\\na|k141_1_1\\t1000000\\n" | gzip > gene_tpm.tsv.gz
    printf "sample\\tn_genes_detected\\ttotal_depth\\ns1\\t1\\t10\\n" > gene_abundance_stats.tsv
    """
}
