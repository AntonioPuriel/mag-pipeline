process MIMAG_QUALITY {
    label 'process_low'
    // Only needs Python 3 (standard library)
    container 'python:3.12-slim'
    publishDir "${params.outdir}/05_quality", mode: 'copy'

    input:
    path bin_summary
    path checkm2, stageAs: 'checkm2/*'
    path rrna,    stageAs: 'rrna/*'
    path trna,    stageAs: 'trna/*'
    path gunc,    stageAs: 'gunc/*'
    path gtdbtk,  stageAs: 'gtdbtk/*'

    output:
    path "mimag_quality.tsv", emit: tsv
    tuple val("${task.process}"), val('python'), eval("python3 --version 2>&1 | grep -oE '[0-9]+(\\.[0-9]+)+[0-9a-z.-]*' | head -n1 || true"), topic: versions

    // No stub block: the table is built from the stub outputs of the steps
    // before it, so stub runs test the real script
    script:
    """
    mimag_quality.py \\
        --bin-summary ${bin_summary} \\
        --checkm2 checkm2/* \\
        --rrna rrna/* \\
        --trna trna/* \\
        --gunc gunc/* \\
        --gtdbtk gtdbtk/* \\
        > mimag_quality.tsv
    """
}
