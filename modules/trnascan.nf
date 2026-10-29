process TRNASCAN {
    tag "${id} ${binner}"
    label 'process_medium'
    // Check the exact biocontainers tag at https://quay.io/repository/biocontainers/trnascan-se?tab=tags
    container 'quay.io/biocontainers/trnascan-se:2.0.12--pl5321h7b50bb2_2'
    publishDir "${params.outdir}/05_quality/mimag", mode: 'copy'

    input:
    tuple val(id), val(binner), path(bins_dir)

    output:
    tuple val(id), val(binner), path("${id}_${binner}.trna.tsv"), emit: summary
    path "trnascan_out",                                          emit: raw

    script:
    """
    mkdir -p trnascan_out
    find -L ${bins_dir} -maxdepth 1 -name '*.fa' | sort > genomes.txt

    # tRNA genes with the bacterial (-B) and the archaeal (-A) models, for the
    # same reason as barrnap. One thread per search and ${task.cpus} genomes at
    # once: Infernal scales better across genomes than within one
    xargs -r -P ${task.cpus} -I{} sh -c 'set -e
        g=\$(basename {} .fa)
        tRNAscan-SE -B --brief --quiet --forceow --thread 1 -o trnascan_out/\$g.bac.txt {} > /dev/null
        tRNAscan-SE -A --brief --quiet --forceow --thread 1 -o trnascan_out/\$g.arc.txt {} > /dev/null
    ' < genomes.txt

    trna_summary.sh genomes.txt trnascan_out > ${id}_${binner}.trna.tsv
    """

    stub:
    """
    mkdir -p trnascan_out
    printf "genome\\tmodel\\ttrna_types\\ttrna_total\\ttrna_missing\\n" > ${id}_${binner}.trna.tsv
    for f in \$(find -L ${bins_dir} -maxdepth 1 -name '*.fa'); do
        g=\$(basename \$f .fa)
        printf "\$g\\tbac\\t20\\t45\\t-\\n\$g\\tarc\\t17\\t30\\tAsn,Gln,Val\\n" >> ${id}_${binner}.trna.tsv
    done
    """
}
