/*
 * Taxonomic profile of every sample from its reads, with SingleM.
 *
 * Independent of assembly and binning: SingleM finds reads of single-copy
 * marker genes and classifies them against its metapackage (GTDB). It sees
 * the organisms that were never assembled or binned, which is what the MAGs
 * miss (in Marine Honour, almost all archaea). It is also the profile the
 * original study used to check its MAGs.
 *
 * --singlem_db is the folder written by `singlem data --output-directory`, or
 * the metapackage folder inside it (the one with CONTENTS.json).
 */

process SINGLEM {
    tag "${sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/singlem:0.21.3--pyhdfd78af_0'
    publishDir "${params.outdir}/09_singlem/profiles", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)
    path db

    output:
    path "${sample}.singlem.tsv", emit: profile

    script:
    """
    set -euo pipefail

    # SingleM wants the folder that holds CONTENTS.json. `singlem data` puts it
    # one level down (<dir>/S<version>.GTDB_r<release>.metapackage_<date>.smpkg.zb),
    # so accept either folder
    db_dir=\$(readlink -f ${db})
    if [ -f "\$db_dir/CONTENTS.json" ]; then
        metapackage=\$db_dir
    else
        shopt -s nullglob
        found=( "\$db_dir"/*/CONTENTS.json )
        shopt -u nullglob
        if [ "\${#found[@]}" -ne 1 ]; then
            echo "--singlem_db: expected one metapackage (a folder with CONTENTS.json) in \$db_dir, found \${#found[@]}" >&2
            exit 1
        fi
        metapackage=\$(dirname "\${found[0]}")
    fi
    echo "SingleM metapackage: \$metapackage" >&2
    export SINGLEM_METAPACKAGE_PATH=\$metapackage

    # Temporary files in the task folder, not in the node's /tmp: on some
    # clusters /tmp is in memory, and SingleM writes large chunks there
    mkdir -p tmp
    export TMPDIR=\$PWD/tmp

    # SingleM runs one DIAMOND per thread, each with its own memory: half as
    # many threads as cores leaves about 4 GB per DIAMOND (with 16 GB and 8
    # threads they ran out of memory on 12-16 Gbase samples)
    singlem pipe \\
        -1 ${r1} \\
        -2 ${r2} \\
        -p raw_profile.tsv \\
        --threads ${Math.max(1, task.cpus.intdiv(2))}

    rm -rf tmp

    # SingleM names the sample after the read file: use the samplesheet name
    awk -v s="${sample}" 'BEGIN { FS = OFS = "\\t" } NR == 1 { print; next } { \$1 = s; print }' \\
        raw_profile.tsv > ${sample}.singlem.tsv
    """

    stub:
    """
    printf "sample\\tcoverage\\ttaxonomy\\n${sample}\\t10.0\\tRoot; d__Bacteria\\n${sample}\\t1.0\\tRoot; d__Archaea\\n" > ${sample}.singlem.tsv
    """
}

/*
 * The profiles of all samples together: one table with coverage and relative
 * abundance of every taxon, and site-by-taxon tables of relative abundance at
 * each rank from domain to species.
 */
process SINGLEM_SUMMARY {
    label 'process_low'
    container 'quay.io/biocontainers/singlem:0.21.3--pyhdfd78af_0'
    publishDir "${params.outdir}/09_singlem", mode: 'copy'

    input:
    path profiles

    output:
    path "singlem_*.tsv", emit: tables

    script:
    """
    set -euo pipefail
    files=\$(ls *.singlem.tsv | sort)

    singlem summarise \\
        --input-taxonomic-profiles \$files \\
        --output-taxonomic-profile-with-extras singlem_profiles_with_extras.tsv

    for level in domain phylum class order family genus species; do
        singlem summarise \\
            --input-taxonomic-profiles \$files \\
            --output-species-by-site-relative-abundance singlem_\${level}.tsv \\
            --output-species-by-site-level \$level
    done
    """

    stub:
    """
    for t in profiles_with_extras domain phylum class order family genus species; do
        printf "sample\\ttaxonomy\\n" > singlem_\$t.tsv
    done
    """
}
