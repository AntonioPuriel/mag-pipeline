/*
 * De novo phylogenomic tree of the final genomes (--mag_tree).
 *
 * GTDB-Tk classify_wf does not align the genomes it assigns to a species by
 * ANI (skani), so its alignment and its placement trees miss every genome of
 * a known species. Here the GTDB-Tk marker genes of all genomes are found
 * (identify), aligned without the GTDB reference genomes (align
 * --skip_gtdb_refs), and a tree is inferred with FastTree (infer: WAG + Gamma).
 *
 * One tree per marker set: bac120 for bacteria, ar53 for archaea. A tree needs
 * at least 3 genomes, so a domain with fewer is skipped (Marine Honour has a
 * single archaeon).
 *
 * Published in 06_taxonomy/mag_tree/:
 *   <marker>.tree                 unrooted tree, leaves named after the genomes
 *   <marker>.user_msa.fasta.gz    the alignment it was built from
 *   <marker>.filtered.tsv         genomes left out by GTDB-Tk (too few markers)
 *   mag_tree.log                  genomes per marker set, and trees skipped
 */

process MAG_TREE {
    label 'process_medium'
    container 'quay.io/biocontainers/gtdbtk:2.4.0--pyhdfd78af_2'
    publishDir "${params.outdir}/06_taxonomy", mode: 'copy'

    input:
    path bin_dirs, stageAs: 'bins/set??'
    path db

    output:
    path "mag_tree/*.tree",              emit: trees,    optional: true
    path "mag_tree/*.user_msa.fasta.gz", emit: msa,      optional: true
    path "mag_tree/*.filtered.tsv",      emit: filtered, optional: true
    path "mag_tree/mag_tree.log",        emit: log

    script:
    """
    set -euo pipefail
    export GTDBTK_DATA_PATH=\$(readlink -f ${db})

    # Every genome in one folder. Bin names are unique across assemblies
    # (they carry the assembly), so they never collide.
    mkdir -p genomes mag_tree
    for fa in bins/*/*.fa; do
        [ -e "\$fa" ] || continue
        ln -s "\$(readlink -f "\$fa")" genomes/
    done
    n_genomes=\$(find -L genomes -name '*.fa' | wc -l)
    echo "Genomes: \$n_genomes" > mag_tree/mag_tree.log

    if [ "\$n_genomes" -lt 3 ]; then
        echo "Fewer than 3 genomes: no tree" | tee -a mag_tree/mag_tree.log >&2
        exit 0
    fi

    gtdbtk identify \\
        --genome_dir genomes \\
        --extension fa \\
        --out_dir denovo \\
        --cpus ${task.cpus}

    gtdbtk align \\
        --identify_dir denovo \\
        --out_dir denovo \\
        --skip_gtdb_refs \\
        --cpus ${task.cpus}

    for marker in bac120 ar53; do
        shopt -s nullglob
        msa=( denovo/align/gtdbtk.\${marker}.user_msa.fasta* )
        shopt -u nullglob

        filtered=denovo/align/gtdbtk.\${marker}.filtered.tsv
        if [ -f "\$filtered" ]; then
            cp "\$filtered" mag_tree/\${marker}.filtered.tsv
        fi
        if [ "\${#msa[@]}" -eq 0 ]; then
            echo "\${marker}: no genomes" >> mag_tree/mag_tree.log
            continue
        fi

        n_seq=\$(zcat -f "\${msa[0]}" | grep -c '^>' || true)
        echo "\${marker}: \${n_seq} genomes in the alignment" >> mag_tree/mag_tree.log
        if [ "\$n_seq" -lt 3 ]; then
            echo "\${marker}: fewer than 3 genomes, tree skipped" >> mag_tree/mag_tree.log
            continue
        fi

        gtdbtk infer \\
            --msa_file "\${msa[0]}" \\
            --out_dir infer_\${marker} \\
            --gamma \\
            --cpus ${task.cpus}

        # infer writes the tree in intermediate_results and links it from the
        # top of its folder: publish one real copy
        cp -L infer_\${marker}/infer/intermediate_results/gtdbtk.unrooted.tree mag_tree/\${marker}.tree
        zcat -f "\${msa[0]}" | gzip > mag_tree/\${marker}.user_msa.fasta.gz
    done

    cat mag_tree/mag_tree.log >&2
    """

    stub:
    """
    mkdir -p mag_tree
    printf "(a:0.1,b:0.2,c:0.3);\\n" > mag_tree/bac120.tree
    printf ">a\\nMK\\n>b\\nMK\\n>c\\nMK\\n" | gzip > mag_tree/bac120.user_msa.fasta.gz
    printf "Genomes: 3\\nbac120: 3 genomes in the alignment\\nar53: no genomes\\n" > mag_tree/mag_tree.log
    """
}
