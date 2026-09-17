process CONCOCT {
    tag "${id}"
    label 'process_medium'
    container 'quay.io/biocontainers/concoct:1.1.0--py39h8907335_8'
    publishDir "${params.outdir}/04_binning/concoct", mode: 'copy'

    input:
    tuple val(id), path(contigs), path(bams), path(bais)

    output:
    tuple val(id), val('concoct'), path("${id}_concoct_bins"), emit: bins

    script:
    """
    gzip -dc ${contigs} > contigs.fa

    # 1. Cut contigs into 10 kb pieces
    cut_up_fasta.py contigs.fa \\
        --chunk_size 10000 \\
        --overlap_size 0 \\
        --merge_last \\
        --bedfile contigs_10K.bed \\
        > contigs_10K.fa

    # 2. Coverage of each piece in each sample
    concoct_coverage_table.py contigs_10K.bed ${bams} > coverage_table.tsv

    # 3. Clustering
    concoct \\
        --composition_file contigs_10K.fa \\
        --coverage_file coverage_table.tsv \\
        --length_threshold ${params.min_contig_length} \\
        --threads ${task.cpus} \\
        --seed 42 \\
        --basename concoct_out/

    # 4. Merge the pieces back into their original contigs
    merge_cutup_clustering.py concoct_out/clustering_gt${params.min_contig_length}.csv \\
        > clustering_merged.csv

    # 5. One FASTA file per bin, named <assembly>.concoct.<cluster>.fa
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
