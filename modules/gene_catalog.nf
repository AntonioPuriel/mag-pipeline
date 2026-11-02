/*
 * Non-redundant gene catalog of all assemblies (--workflow contigs).
 *
 * Genes are clustered on their nucleotide sequence with MMseqs2 linclust, the
 * criteria of reference gene catalogs: >= 95 % identity over >= 90 % of the
 * shorter gene. --cluster-mode 2 makes the longest gene of each cluster its
 * representative, and --cov-mode 1 asks that the representative cover 90 % of
 * every member. linclust runs in linear time and memory, which keeps several
 * million genes within a standard node.
 *
 * Outputs: the representatives as nucleotide and protein sequences, and the
 * representative -> member table that GENE_ABUNDANCE uses to give each
 * representative the coverage of its members.
 */

process GENE_CATALOG {
    label 'process_high'
    container 'quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0'
    publishDir "${params.outdir}/08_contigs/catalog", mode: 'copy'

    input:
    path fna, stageAs: 'fna/*'
    path faa, stageAs: 'faa/*'

    output:
    path "gene_catalog.fna.gz",          emit: fna
    path "gene_catalog.faa.gz",          emit: proteins
    path "gene_catalog.clusters.tsv.gz", emit: clusters
    path "gene_catalog.stats.tsv",       emit: stats

    script:
    // leave a margin below the task memory for what MMseqs2 does not count
    def mem_limit = task.memory ? "--split-memory-limit ${(task.memory.toMega() * 0.8) as long}M" : ''
    """
    set -euo pipefail

    # Genes of every assembly in a fixed order, keeping those of at least
    # ${params.catalog_min_gene_length} nt (gene IDs only)
    zcat \$(ls fna/*.fna.gz | sort) \\
        | awk -v min=${params.catalog_min_gene_length} '
              function flush() { if (h != "") { n++; if (length(s) >= min) { k++; print h; print s } } }
              /^>/ { flush(); h = \$1; s = ""; next }
              { s = s \$0 }
              END { flush(); print n + 0 > "n_called.txt"; print k + 0 > "n_kept.txt" }' \\
        > genes.fna

    mmseqs easy-linclust genes.fna catalog tmp \\
        --min-seq-id ${params.catalog_min_id} \\
        -c ${params.catalog_min_cov} \\
        --cov-mode 1 \\
        --cluster-mode 2 \\
        --threads ${task.cpus} \\
        ${mem_limit} \\
        > mmseqs.log

    # Proteins of the representatives, in the order of the input files
    cut -f1 catalog_cluster.tsv | sort -u > reps.txt
    zcat \$(ls faa/*.faa.gz | sort) \\
        | awk 'NR == FNR { keep[\$1]; next }
               /^>/ { id = substr(\$1, 2); p = (id in keep); if (p) print ">" id; next }
               p' reps.txt - \\
        | gzip > gene_catalog.faa.gz

    n_reps=\$(wc -l < reps.txt)
    n_prot=\$(zcat gene_catalog.faa.gz | grep -c '^>' || true)
    if [ "\$n_reps" -ne "\$n_prot" ]; then
        echo "\$n_reps representatives but \$n_prot proteins" >&2
        exit 1
    fi

    gzip -c catalog_rep_seq.fasta > gene_catalog.fna.gz
    gzip -c catalog_cluster.tsv > gene_catalog.clusters.tsv.gz

    n_single=\$(cut -f1 catalog_cluster.tsv | uniq -c | awk '\$1 == 1' | wc -l)
    printf "genes_called\\tgenes_kept\\tcatalog_genes\\tsingletons\\n%s\\t%s\\t%s\\t%s\\n" \\
        "\$(cat n_called.txt)" "\$(cat n_kept.txt)" "\$n_reps" "\$n_single" \\
        > gene_catalog.stats.tsv

    rm -rf tmp genes.fna catalog_all_seqs.fasta catalog_rep_seq.fasta catalog_cluster.tsv
    """

    stub:
    """
    printf ">a|k141_1_1\\nATGGCAAAAGTA\\n" | gzip > gene_catalog.fna.gz
    printf ">a|k141_1_1\\nMAKV\\n" | gzip > gene_catalog.faa.gz
    printf "a|k141_1_1\\ta|k141_1_1\\n" | gzip > gene_catalog.clusters.tsv.gz
    printf "genes_called\\tgenes_kept\\tcatalog_genes\\tsingletons\\n1\\t1\\t1\\t1\\n" > gene_catalog.stats.tsv
    """
}

/*
 * The catalog proteins in pieces, so that eggNOG-mapper, CANT-HYD and HADEG
 * annotate them in parallel jobs. Pieces are named catalog_chunkNNN.faa.
 */
process SPLIT_CATALOG {
    label 'process_low'

    input:
    path proteins

    output:
    path "catalog_chunk*.faa", emit: chunks

    script:
    """
    set -euo pipefail
    zcat ${proteins} \\
        | awk -v n=${params.catalog_chunk_size} '
              /^>/ { if (c % n == 0) { if (f) close(f); f = sprintf("catalog_chunk%03d.faa", c / n + 1) }; c++ }
              { print > f }'
    """

    stub:
    """
    printf ">a|k141_1_1\\nMAKV\\n" > catalog_chunk001.faa
    """
}
