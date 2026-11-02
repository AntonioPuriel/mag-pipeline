/*
 * Genes of every filtered contig of one assembly (--workflow contigs).
 *
 * Prodigal runs in metagenomic mode (-p meta): contigs come from many
 * organisms, so there is no single genome to train on. Each contig is
 * predicted on its own in that mode, so the contigs are split in as many
 * pieces as CPUs and predicted in parallel; the genes are the same as in a
 * single run.
 *
 * Gene names become <assembly>|<contig>_<n>: unique across assemblies (MEGAHIT
 * reuses contig names), and the contig and assembly of every gene can be read
 * back from its name.
 */

process PRODIGAL_CONTIGS {
    tag "${id}"
    label 'process_medium'
    container 'quay.io/biocontainers/prodigal:2.6.3--h577a1d6_11'
    publishDir "${params.outdir}/08_contigs/genes", mode: 'copy'

    input:
    tuple val(id), path(contigs)

    output:
    tuple val(id), path("${id}.genes.faa.gz"), path("${id}.genes.fna.gz"), emit: genes
    path "${id}.genes.stats.tsv",                                           emit: stats

    script:
    """
    set -euo pipefail

    gzip -dc ${contigs} > contigs.fa

    # contigs dealt in turn to ${task.cpus} pieces
    awk -v n=${task.cpus} '/^>/ { if (f) close(f); i = (i + 1) % n; f = sprintf("part_%03d.fa", i) }
                           { print >> f }' contigs.fa

    ls part_*.fa | xargs -P ${task.cpus} -I {} \\
        prodigal -i {} -a {}.faa -d {}.fna -o /dev/null -p meta -q

    cat part_*.fa.faa | sed 's/^>/>${id}|/' | gzip > ${id}.genes.faa.gz
    cat part_*.fa.fna | sed 's/^>/>${id}|/' | gzip > ${id}.genes.fna.gz

    # contigs, contigs with genes, genes and complete genes (partial=00)
    n_contigs=\$(grep -c '^>' contigs.fa)
    gzip -dc ${id}.genes.faa.gz \\
        | awk -v a=${id} -v nc=\$n_contigs 'BEGIN { OFS = "\\t" }
              /^>/ { n++; if (\$0 ~ /partial=00/) full++
                     c = substr(\$1, 2); sub(/_[0-9]+\$/, "", c); seen[c] = 1 }
              END { for (c in seen) nwith++
                    print "assembly", "n_contigs", "n_contigs_with_genes", "n_genes", "n_complete_genes"
                    print a, nc, nwith + 0, n + 0, full + 0 }' \\
        > ${id}.genes.stats.tsv

    rm -f contigs.fa part_*
    """

    stub:
    """
    printf ">${id}|k141_1_1 # 1 # 12 # 1 # ID=1_1;partial=00\\nMAKV\\n" | gzip > ${id}.genes.faa.gz
    printf ">${id}|k141_1_1 # 1 # 12 # 1 # ID=1_1;partial=00\\nATGGCAAAAGTA\\n" | gzip > ${id}.genes.fna.gz
    printf "assembly\\tn_contigs\\tn_contigs_with_genes\\tn_genes\\tn_complete_genes\\n${id}\\t1\\t1\\t1\\t1\\n" > ${id}.genes.stats.tsv
    """
}
