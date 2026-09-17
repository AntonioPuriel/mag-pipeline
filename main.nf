#!/usr/bin/env nextflow
/*
 * mag-pipeline: MAG recovery from paired-end shotgun metagenomes
 * Author: Antonio Puriel Hernández
 */

include { FASTP           } from './modules/fastp'
include { MEGAHIT         } from './modules/megahit'
include { FILTER_CONTIGS  } from './modules/filter_contigs'
include { BOWTIE2_BUILD   } from './modules/bowtie2_build'
include { BOWTIE2_ALIGN   } from './modules/bowtie2_align'
include { CONTIG_DEPTHS   } from './modules/contig_depths'
include { MAPPING_SUMMARY } from './modules/mapping_summary'

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

    // ---- Read mapping ----
    BOWTIE2_BUILD(FILTER_CONTIGS.out.contigs)

    if (params.assembly_mode == 'coassembly') {
        // every sample is mapped back to the co-assembly
        mapping_input = BOWTIE2_BUILD.out.index
            .combine(FASTP.out.reads)
    } else {
        // each sample is mapped back to its own assembly
        mapping_input = BOWTIE2_BUILD.out.index
            .join(FASTP.out.reads)
            .map { id, index, r1, r2 -> tuple(id, index, id, r1, r2) }
    }

    BOWTIE2_ALIGN(mapping_input)
    MAPPING_SUMMARY(BOWTIE2_ALIGN.out.log.collect())

    // ---- Contig coverage per assembly (input for binning) ----
    bams_per_assembly = BOWTIE2_ALIGN.out.bam
        .map { id, sample, bam, bai -> tuple(id, bam, bai) }
        .groupTuple()

    CONTIG_DEPTHS(bams_per_assembly)
}
