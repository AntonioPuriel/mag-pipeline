#!/usr/bin/env nextflow
/*
 * mag-pipeline: MAG recovery from paired-end shotgun metagenomes
 * Author: Antonio Puriel Hernández
 */

include { FASTP          } from './modules/fastp'
include { MEGAHIT        } from './modules/megahit'
include { FILTER_CONTIGS } from './modules/filter_contigs'

workflow {

    // ---- Parameter checks ----
    if (!params.input) {
        error "Please provide a samplesheet with --input"
    }
    if (!(params.assembly_mode in ['coassembly', 'per_sample'])) {
        error "--assembly_mode must be 'coassembly' or 'per_sample' (got '${params.assembly_mode}')"
    }

    // ---- Read samplesheet: sample,fastq_1,fastq_2 ----
    reads_ch = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true))
        }

    // ---- Quality control ----
    FASTP(reads_ch)

    // ---- Assembly input ----
    if (params.assembly_mode == 'coassembly') {
        assembly_input = FASTP.out.reads
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { samples ->
                tuple('coassembly', samples.collect { it[1] }, samples.collect { it[2] })
            }
    } else {
        assembly_input = FASTP.out.reads
            .map { sample, r1, r2 -> tuple(sample, [r1], [r2]) }
    }

    // ---- Assembly and contig filtering ----
    MEGAHIT(assembly_input)
    FILTER_CONTIGS(MEGAHIT.out.contigs)
}
