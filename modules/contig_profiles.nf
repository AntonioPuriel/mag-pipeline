/*
 * Community profiles from the gene catalog (--workflow contigs): TPM of every
 * KEGG orthologue, hydrocarbon marker and taxon in every sample, and one
 * annotation row per catalog gene. Every annotation is optional; absent ones
 * arrive as NO_FILE and their tables are not written.
 */

process CONTIG_PROFILES {
    label 'process_medium'
    // Check the exact tag before using containers: any image with Python 3, numpy and pandas
    container 'quay.io/biocontainers/pandas:2.2.1'
    publishDir "${params.outdir}/08_contigs/profiles", mode: 'copy'

    input:
    path tpm
    path eggnog,  stageAs: 'eggnog/*'
    path canthyd, stageAs: 'canthyd/*'
    path hadeg,   stageAs: 'hadeg/*'
    path cat,     stageAs: 'cat/*'

    output:
    path "*.tsv",                      emit: tables
    path "catalog_annotation.tsv.gz",  emit: annotation

    // No stub block: runs on the stub outputs of the steps before it
    script:
    """
    contig_profiles.py \\
        --tpm ${tpm} \\
        --eggnog eggnog/* \\
        --canthyd canthyd/* \\
        --hadeg hadeg/* \\
        --cat cat/* \\
        --outdir .
    """
}
