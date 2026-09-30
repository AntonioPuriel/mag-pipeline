process GALAH {
    tag "all assemblies"
    label 'process_medium'
    // Check the exact biocontainers tag at https://quay.io/repository/biocontainers/galah?tab=tags
    container 'quay.io/biocontainers/galah:0.5.2'
    publishDir "${params.outdir}/05_quality/dereplication", mode: 'copy'

    input:
    path bins_dirs, stageAs: 'bins_*'
    path reports,   stageAs: 'checkm2_*.tsv'

    output:
    path "representatives",            emit: representatives
    path "dereplication_clusters.tsv", emit: clusters

    script:
    """
    # Every bin from every assembly. galah matches each genome to its CheckM2
    # row by file name, so names must be unique across assemblies.
    find -L bins_* -maxdepth 1 -name '*.fa' | sort > genomes.txt

    dups=\$(xargs -r -n1 basename < genomes.txt | sort | uniq -d)
    if [ -n "\$dups" ]; then
        echo "Bin names are repeated across assemblies, cannot dereplicate:" >&2
        echo "\$dups" >&2
        exit 1
    fi

    printf "representative\\tmember\\n" > dereplication_clusters.tsv

    if [ ! -s genomes.txt ]; then
        echo "No bins to dereplicate" >&2
        mkdir -p representatives
        exit 0
    fi

    # One CheckM2 table for all assemblies, keeping only the first header
    awk 'FNR == 1 && NR != 1 { next } { print }' checkm2_*.tsv > checkm2_all.tsv

    # Bins below the quality thresholds are left out; within each cluster at
    # the ANI threshold, the best bin by CheckM2 quality is the representative
    galah cluster \\
        --genome-fasta-list genomes.txt \\
        --checkm2-quality-report checkm2_all.tsv \\
        --min-completeness ${params.derep_min_completeness} \\
        --max-contamination ${params.derep_max_contamination} \\
        --ani ${params.derep_ani} \\
        --threads ${task.cpus} \\
        --output-cluster-definition galah_clusters.tsv \\
        --output-representative-fasta-directory-copy representatives

    mkdir -p representatives

    # From paths to bin names: representative<TAB>member
    awk -F'\\t' -v OFS='\\t' '{
        for (i = 1; i <= 2; i++) { sub(/.*\\//, "", \$i); sub(/\\.fa\$/, "", \$i) }
        print
    }' galah_clusters.tsv >> dereplication_clusters.tsv
    """

    stub:
    """
    mkdir -p representatives
    printf "representative\\tmember\\n" > dereplication_clusters.tsv
    for f in \$(find -L bins_* -maxdepth 1 -name '*.fa'); do
        cp -L \$f representatives/
        n=\$(basename \$f .fa)
        printf "\$n\\t\$n\\n" >> dereplication_clusters.tsv
    done
    """
}
