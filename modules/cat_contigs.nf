/*
 * Taxonomy of every contig with CAT (CAT_pack contigs).
 *
 * CAT aligns the proteins of each contig with DIAMOND against the reference
 * proteins and classifies the contig by a vote of its genes, weighted by bit
 * score (LCA). It works on contigs of unknown organisms, which GTDB-Tk cannot
 * classify: GTDB-Tk needs a genome with its marker genes.
 *
 * --cat_db is the folder made by `CAT_pack prepare` (subfolders db/ and tax/).
 * Built from GTDB, the taxonomy matches that of GTDB-Tk for the MAGs.
 *
 * The proteins from PRODIGAL_CONTIGS are reused (no second gene prediction).
 * CAT reads the contig of each protein from its name (<contig>_<n>), so the
 * <assembly>| prefix is removed first. The DIAMOND alignment is deleted at the
 * end: it is large and only needed to rerun CAT with other parameters.
 */

process CAT_CONTIGS {
    tag "${id}"
    label 'process_high'
    container 'quay.io/biocontainers/cat:6.0.1--hdfd78af_1'
    publishDir "${params.outdir}/08_contigs/taxonomy", mode: 'copy', pattern: '*.{named.txt,log}'

    input:
    tuple val(id), path(contigs), path(proteins)
    path db

    output:
    tuple val(id), path("${id}.CAT.named.txt"), emit: classification
    path "${id}.CAT.log",                        emit: log

    script:
    """
    set -euo pipefail

    gzip -dc ${contigs} > contigs.fa
    gzip -dc ${proteins} | sed 's/^>${id}|/>/' > proteins.faa

    CAT_pack contigs \\
        -c contigs.fa \\
        -p proteins.faa \\
        -d ${db}/db \\
        -t ${db}/tax \\
        -o ${id}.CAT \\
        -n ${task.cpus} \\
        --block_size ${params.cat_block_size} \\
        --index_chunks ${params.cat_index_chunks} \\
        --tmpdir . \\
        --no_stars

    # superkingdom (domain), phylum, class, order, family, genus, species,
    # each with the fraction of the bit score that supports it
    CAT_pack add_names \\
        -i ${id}.CAT.contig2classification.txt \\
        -o ${id}.CAT.named.txt \\
        -t ${db}/tax \\
        --only_official

    rm -f contigs.fa proteins.faa ${id}.CAT.alignment.diamond*
    """

    stub:
    """
    printf "# contig\\tclassification\\treason\\tlineage\\tlineage scores\\tsuperkingdom\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\n" > ${id}.CAT.named.txt
    printf "k141_1\\ttaxid assigned\\tbased on 1/1 ORFs\\t1;2\\t1.00;1.00\\td__Bacteria: 1.00\\tno support\\tno support\\tno support\\tno support\\tno support\\tno support\\n" >> ${id}.CAT.named.txt
    touch ${id}.CAT.log
    """
}
