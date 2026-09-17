process DAS_TOOL {
    tag "${id}"
    label 'process_medium'
    container 'quay.io/biocontainers/das_tool:1.1.7--r44hdfd78af_1'
    publishDir "${params.outdir}/04_binning/dastool", mode: 'copy'

    input:
    tuple val(id), path(contigs), path(bin_dirs)

    output:
    tuple val(id), val('dastool'), path("${id}_dastool_bins"), emit: bins
    path "${id}_DASTool_summary.tsv", optional: true,         emit: summary
    path "${id}_DASTool.log",                                 emit: log

    script:
    """
    gzip -dc ${contigs} > contigs.fa

    # One contig-to-bin table per binner, skipping binners without bins
    inputs=()
    labels=()
    for dir in ${bin_dirs}; do
        name=\$(basename "\$dir")
        name=\${name%_bins}
        binner=\${name##*_}
        fasta_to_contig2bin.sh "\$dir" > "\${binner}.contigs2bin.tsv"
        if [ -s "\${binner}.contigs2bin.tsv" ]; then
            inputs+=("\${binner}.contigs2bin.tsv")
            labels+=("\$binner")
        fi
    done

    mkdir -p ${id}_dastool_bins

    if [ "\${#inputs[@]}" -eq 0 ]; then
        echo "No bins available for refinement" | tee ${id}_DASTool.log >&2
        exit 0
    fi

    # DAS_Tool exits with an error when no bin passes the score threshold:
    # treat that case as "no refined bins" and fail on any other error
    set +e
    DAS_Tool \\
        --bins "\$(IFS=,; echo "\${inputs[*]}")" \\
        --labels "\$(IFS=,; echo "\${labels[*]}")" \\
        --contigs contigs.fa \\
        --outputbasename dastool/${id} \\
        --score_threshold ${params.dastool_score_threshold} \\
        --write_bins \\
        --threads ${task.cpus}
    status=\$?
    set -e

    if [ -f dastool/${id}_DASTool.log ]; then
        cp dastool/${id}_DASTool.log .
    else
        touch ${id}_DASTool.log
    fi

    if [ "\$status" -ne 0 ]; then
        if grep -q "No bins with bin-score" ${id}_DASTool.log; then
            echo "DAS_Tool: no bins above score threshold ${params.dastool_score_threshold}" >&2
        else
            exit "\$status"
        fi
    fi

    if [ -f dastool/${id}_DASTool_summary.tsv ]; then
        cp dastool/${id}_DASTool_summary.tsv .
    fi
    if [ -d dastool/${id}_DASTool_bins ]; then
        find dastool/${id}_DASTool_bins -name '*.fa' -exec mv {} ${id}_dastool_bins/ \\;
    fi
    """

    stub:
    """
    mkdir -p ${id}_dastool_bins
    printf ">k141_1\\nACGT\\n" > ${id}_dastool_bins/${id}.metabat2.1.fa
    touch ${id}_DASTool_summary.tsv ${id}_DASTool.log
    """
}
