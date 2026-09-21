process MAG_ABUNDANCE {
    label 'process_low'
    publishDir "${params.outdir}/04_binning", mode: 'copy'

    input:
    path bin_dirs
    path depths

    output:
    path "contig2bin.tsv",     emit: contig2bin
    path "mag_abundance.tsv",  emit: abundance

    script:
    """
    # contig -> bin for every bin set, keeping the assembly and binner
    printf "assembly\\tbinner\\tbin\\tcontig\\n" > contig2bin.tsv

    for dir in ${bin_dirs}; do
        name=\$(basename "\$dir")
        name=\${name%_bins}
        binner=\${name##*_}
        assembly=\${name%_*}

        fasta_to_contig2bin.sh "\$dir" \\
            | awk -v a="\$assembly" -v b="\$binner" -v OFS='\\t' '{ print a, b, \$2, \$1 }' \\
            >> contig2bin.tsv
    done

    # Length-weighted mean depth of every MAG in every sample
    mag_abundance.sh contig2bin.tsv ${depths} > mag_abundance.tsv
    """

    stub:
    """
    printf "assembly\\tbinner\\tbin\\tcontig\\n" > contig2bin.tsv
    printf "assembly\\tbinner\\tbin\\tsample\\tmean_depth\\tbin_length\\n" > mag_abundance.tsv
    """
}
