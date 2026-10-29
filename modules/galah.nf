process GALAH {
    tag "all assemblies"
    label 'process_medium'
    // Check the exact biocontainers tag at https://quay.io/repository/biocontainers/galah?tab=tags
    container 'quay.io/biocontainers/galah:0.5.2'
    publishDir "${params.outdir}/05_quality/dereplication", mode: 'copy'

    input:
    path bins_dirs, stageAs: 'bins_*'
    path reports,   stageAs: 'checkm2_*.tsv'
    path gunc,      stageAs: 'gunc/*'

    output:
    path "representatives",            emit: representatives
    path "dereplication_clusters.tsv", emit: clusters
    path "gunc_excluded.txt",          emit: gunc_excluded

    script:
    def exclude = params.gunc_exclude ? 'true' : 'false'
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

    # With --gunc_exclude, genomes that GUNC flags as chimeric cannot be
    # chosen as representatives: they are removed before clustering
    : > gunc_excluded.txt
    if [ "${exclude}" = "true" ] && ls gunc/*.gunc.tsv >/dev/null 2>&1; then
        awk -F'\\t' 'FNR == 1 { g = p = 0; for (i = 1; i <= NF; i++) { if (\$i == "genome") g = i; if (\$i == "pass.GUNC") p = i }; next }
                     g && p && \$p == "False" { print \$g }' gunc/*.gunc.tsv | sort -u > gunc_excluded.txt
        if [ -s gunc_excluded.txt ]; then
            echo "Left out of dereplication as chimeric (GUNC): \$(wc -l < gunc_excluded.txt) genomes" >&2
            awk 'NR == FNR { bad[\$1]; next } { n = \$0; sub(/.*\\//, "", n); sub(/\\.fa\$/, "", n) } !(n in bad)' \\
                gunc_excluded.txt genomes.txt > genomes_kept.txt
            mv genomes_kept.txt genomes.txt
        fi
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
    : > gunc_excluded.txt
    printf "representative\\tmember\\n" > dereplication_clusters.tsv
    for f in \$(find -L bins_* -maxdepth 1 -name '*.fa'); do
        cp -L \$f representatives/
        n=\$(basename \$f .fa)
        printf "\$n\\t\$n\\n" >> dereplication_clusters.tsv
    done
    """
}
