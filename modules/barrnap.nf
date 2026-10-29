process BARRNAP {
    tag "${id} ${binner}"
    label 'process_medium'
    // Check the exact biocontainers tag at https://quay.io/repository/biocontainers/barrnap?tab=tags
    container 'quay.io/biocontainers/barrnap:0.9--hdfd78af_4'
    publishDir "${params.outdir}/05_quality/mimag", mode: 'copy'

    input:
    tuple val(id), val(binner), path(bins_dir)

    output:
    tuple val(id), val(binner), path("${id}_${binner}.rrna.tsv"), emit: summary
    path "barrnap_gff",                                           emit: gff

    script:
    """
    mkdir -p barrnap_gff
    find -L ${bins_dir} -maxdepth 1 -name '*.fa' | sort > genomes.txt

    # rRNA genes (5S, 16S, 23S) with the bacterial and the archaeal models:
    # the genomes have no taxonomy yet, so the domain is chosen later, in the
    # MIMAG table. barrnap is fast, so each genome uses one thread and
    # ${task.cpus} genomes run at once
    xargs -r -P ${task.cpus} -I{} sh -c 'set -e
        g=\$(basename {} .fa)
        for k in bac arc; do
            barrnap --kingdom \$k --threads 1 --quiet {} > barrnap_gff/\$g.\$k.gff
        done
    ' < genomes.txt

    rrna_summary.sh genomes.txt barrnap_gff > ${id}_${binner}.rrna.tsv
    """

    stub:
    """
    mkdir -p barrnap_gff
    printf "genome\\tmodel\\trrna_5S\\trrna_16S\\trrna_23S\\trrna_partial\\n" > ${id}_${binner}.rrna.tsv
    for f in \$(find -L ${bins_dir} -maxdepth 1 -name '*.fa'); do
        g=\$(basename \$f .fa)
        printf "\$g\\tbac\\t1\\t1\\t1\\t0\\n\$g\\tarc\\t0\\t1\\t0\\t0\\n" >> ${id}_${binner}.rrna.tsv
    done
    """
}
