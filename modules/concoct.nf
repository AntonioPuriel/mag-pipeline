process CONCOCT {
    tag "${id}"
    label 'process_medium'
    container 'quay.io/biocontainers/concoct:1.1.0--py39h8907335_8'
    publishDir "${params.outdir}/04_binning/concoct", mode: 'copy'

    input:
    // cutup and coverage come from CONCOCT_CUTUP and MERGE_CONCOCT_COV (modules/map_depth.nf)
    tuple val(id), path(contigs), path(cutup_fa), path(coverage)

    output:
    tuple val(id), val('concoct'), path("${id}_concoct_bins"), emit: bins

    script:
    """
    gzip -dc ${contigs} > contigs.fa

    # 1. Clustering of the 10 kb pieces
    concoct \\
        --composition_file ${cutup_fa} \\
        --coverage_file ${coverage} \\
        --length_threshold ${params.min_contig_length} \\
        --threads ${task.cpus} \\
        --seed 42 \\
        --basename concoct_out/

    # 2. Merge the pieces back into their original contigs
    merge_cutup_clustering.py concoct_out/clustering_gt${params.min_contig_length}.csv \\
        > clustering_merged.csv

    # 3. One FASTA file per bin, named <assembly>.concoct.<cluster>.fa
    mkdir -p ${id}_concoct_bins
    extract_fasta_bins.py contigs.fa clustering_merged.csv --output_path ${id}_concoct_bins

    for fa in ${id}_concoct_bins/*.fa; do
        [ -e "\$fa" ] || continue
        mv "\$fa" "${id}_concoct_bins/${id}.concoct.\$(basename "\$fa")"
    done
    """

    stub:
    """
    mkdir -p ${id}_concoct_bins
    printf ">k141_1\\nACGT\\n" > ${id}_concoct_bins/${id}.concoct.0.fa
    """
}
