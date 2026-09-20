process PRODIGAL {
    tag "${id} ${binner}"
    label 'process_low'
    container 'quay.io/biocontainers/prodigal:2.6.3--h7b50bb2_10'
    publishDir "${params.outdir}/07_function/genes", mode: 'copy'

    input:
    tuple val(id), val(binner), path(bins_dir)

    output:
    tuple val(id), val(binner), path("${id}_${binner}.faa"), emit: proteins
    path "${id}_${binner}.gene_counts.tsv",                  emit: counts

    script:
    """
    : > ${id}_${binner}.faa
    printf "bin\\tn_genes\\n" > ${id}_${binner}.gene_counts.tsv

    for fa in ${bins_dir}/*.fa; do
        [ -e "\$fa" ] || continue
        bin=\$(basename "\$fa" .fa)

        # -p single: the bins are genomes, not raw metagenomic contigs
        prodigal \\
            -i "\$fa" \\
            -a "\$bin.faa" \\
            -d "\$bin.fna" \\
            -o /dev/null \\
            -p single \\
            -q

        # Keep the bin of origin in the sequence name: <bin>|<gene>
        sed "s|^>|>\$bin\\||" "\$bin.faa" >> ${id}_${binner}.faa
        printf "%s\\t%d\\n" "\$bin" "\$(grep -c '^>' "\$bin.faa")" >> ${id}_${binner}.gene_counts.tsv
        rm -f "\$bin.faa" "\$bin.fna"
    done
    """

    stub:
    """
    printf ">%s|k141_1_1 # gene\\nMAKV\\n" "${id}.${binner}.1" > ${id}_${binner}.faa
    printf "bin\\tn_genes\\n${id}.${binner}.1\\t1\\n" > ${id}_${binner}.gene_counts.tsv
    """
}
