#!/usr/bin/env nextflow
/*
 * Runs only FASTP (same module and config as main.nf) and writes
 * trimmed_paths.csv with the trimmed reads of every sample.
 */

include { FASTP } from './modules/fastp'

workflow {
    samplesheet = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true), row.group ?: 'all')
        }

    FASTP(samplesheet.map { sample, r1, r2, group -> tuple(sample, r1, r2) })

    FASTP.out.reads
        .join(samplesheet.map { sample, r1, r2, group -> tuple(sample, group) })
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { rows ->
            'sample,fastq_1,fastq_2,group\n' + rows.collect { it.join(',') }.join('\n') + '\n'
        }
        .collectFile(name: 'trimmed_paths.csv', storeDir: "${launchDir}")
}
