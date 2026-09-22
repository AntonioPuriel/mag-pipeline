/*
 * Mapping without persistent BAMs.
 *
 * MAP_DEPTH maps one sample against one assembly and computes everything that
 * needs the BAM while it exists: the jgi depth for MetaBAT2 and the bedcov
 * coverage for CONCOCT. The BAM is then deleted. It is never declared as an
 * output, so -resume stays safe.
 *
 * MERGE_DEPTHS and MERGE_CONCOCT_COV rebuild the per-assembly tables that
 * CONTIG_DEPTHS and concoct_coverage_table.py produced in v1.4.0.
 */

process MAP_DEPTH {
    tag "${sample} -> ${id}"
    label 'process_medium'
    maxForks params.max_mapping_jobs
    // needs bowtie2, samtools and jgi_summarize_bam_contig_depths (metabat2)
    container 'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6'
    publishDir "${params.outdir}/03_mapping/logs", mode: 'copy', pattern: '*.log'

    input:
    tuple val(id), path(index), path(cutup_bed), val(sample), path(r1), path(r2)

    output:
    tuple val(id), val(sample), path("${id}__${sample}.depth.txt"),  emit: depth
    tuple val(id), val(sample), path("${id}__${sample}.concoct.tsv"), emit: concoct_cov, optional: true
    path "${id}__${sample}.bowtie2.log",                              emit: log

    script:
    """
    set -o pipefail

    # -F 4 drops unmapped reads before sorting: smaller temporary files,
    # same depths (unmapped reads never count towards coverage)
    bowtie2 \\
        -x ${index}/${id} \\
        -1 ${r1} \\
        -2 ${r2} \\
        --threads ${task.cpus} \\
        2> ${id}__${sample}.bowtie2.log \\
    | samtools view -u -F 4 - \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -T ${sample}.sorttmp \\
        -o ${sample}.sorted.bam \\
        -

    samtools index ${sample}.sorted.bam

    # MetaBAT2 depth of this sample. The BAM keeps its v1.4.0 name, so the
    # column headers (<sample>.sorted.bam, <sample>.sorted.bam-var) do not change
    jgi_summarize_bam_contig_depths \\
        --outputDepth ${id}__${sample}.depth.txt \\
        ${sample}.sorted.bam

    # CONCOCT coverage of this sample: what concoct_coverage_table.py computes
    # (samtools bedcov sum divided by the length of each 10 kb piece)
    if [ "${cutup_bed.name}" != "NO_FILE" ]; then
        samtools bedcov ${cutup_bed} ${sample}.sorted.bam \\
            | awk -v s="${sample}" 'BEGIN { OFS = "\\t"; print "contig", "cov_mean_sample_" s }
                                    { printf "%s\\t%.3f\\n", \$4, \$NF / (\$3 - \$2) }' \\
            > ${id}__${sample}.concoct.tsv
    fi

    rm -f ${sample}.sorted.bam ${sample}.sorted.bam.bai
    """

    stub:
    """
    printf "contigName\\tcontigLen\\ttotalAvgDepth\\t${sample}.sorted.bam\\t${sample}.sorted.bam-var\\n" > ${id}__${sample}.depth.txt
    touch ${id}__${sample}.bowtie2.log
    if [ "${cutup_bed.name}" != "NO_FILE" ]; then
        printf "contig\\tcov_mean_sample_${sample}\\n" > ${id}__${sample}.concoct.tsv
    fi
    """
}

process MERGE_DEPTHS {
    tag "${id}"
    label 'process_low'
    publishDir "${params.outdir}/03_mapping", mode: 'copy'

    input:
    tuple val(id), path(depths)

    output:
    tuple val(id), path("${id}.depth.txt"), emit: depth

    script:
    """
    set -o pipefail

    # sort by file name so the sample columns always come in the same order
    files=( \$(printf '%s\\n' ${depths} | sort) )

    cut -f1,2 "\${files[0]}" > base.tsv
    for f in "\${files[@]}"; do
        if ! cmp -s <(cut -f1 "\$f") <(cut -f1 "\${files[0]}"); then
            echo "Contig order differs in \$f" >&2
            exit 1
        fi
        cut -f4,5 "\$f" > "\$f.cols"
    done

    # totalAvgDepth = sum of the per-sample mean depths
    paste base.tsv \$(printf '%s.cols ' "\${files[@]}") \\
        | awk 'BEGIN { FS = OFS = "\\t" }
               NR == 1 { printf "%s\\t%s\\ttotalAvgDepth", \$1, \$2
                         for (i = 3; i <= NF; i++) printf "\\t%s", \$i
                         print ""; next }
               { t = 0
                 for (i = 3; i <= NF; i += 2) t += \$i
                 printf "%s\\t%s\\t%.4f", \$1, \$2, t
                 for (i = 3; i <= NF; i++) printf "\\t%s", \$i
                 print "" }' \\
        > ${id}.depth.txt
    """

    stub:
    """
    touch ${id}.depth.txt
    """
}

process MERGE_CONCOCT_COV {
    tag "${id}"
    label 'process_low'

    input:
    tuple val(id), path(covs)

    output:
    tuple val(id), path("${id}.concoct_coverage.tsv"), emit: coverage

    script:
    """
    set -o pipefail

    files=( \$(printf '%s\\n' ${covs} | sort) )

    cp "\${files[0]}" merged.tsv
    for f in "\${files[@]:1}"; do
        if ! cmp -s <(cut -f1 "\$f") <(cut -f1 "\${files[0]}"); then
            echo "Contig order differs in \$f" >&2
            exit 1
        fi
        paste merged.tsv <(cut -f2 "\$f") > tmp.tsv
        mv tmp.tsv merged.tsv
    done

    mv merged.tsv ${id}.concoct_coverage.tsv
    """

    stub:
    """
    touch ${id}.concoct_coverage.tsv
    """
}

process CONCOCT_CUTUP {
    tag "${id}"
    label 'process_low'
    container 'quay.io/biocontainers/concoct:1.1.0--py39h8907335_8'

    input:
    tuple val(id), path(contigs)

    output:
    tuple val(id), path("${id}_contigs_10K.fa"), path("${id}_contigs_10K.bed"), emit: cutup

    script:
    """
    gzip -dc ${contigs} > contigs.fa

    cut_up_fasta.py contigs.fa \\
        --chunk_size 10000 \\
        --overlap_size 0 \\
        --merge_last \\
        --bedfile ${id}_contigs_10K.bed \\
        > ${id}_contigs_10K.fa

    rm contigs.fa
    """

    stub:
    """
    printf ">k141_1.concoct_part_0\\nACGT\\n" > ${id}_contigs_10K.fa
    printf "k141_1\\t0\\t4\\tk141_1.concoct_part_0\\n" > ${id}_contigs_10K.bed
    """
}
