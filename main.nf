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
include { METABAT2        } from './modules/metabat2'
include { CONCOCT         } from './modules/concoct'
include { BIN_SUMMARY     } from './modules/bin_summary'
include { DAS_TOOL        } from './modules/dastool'
include { CHECKM2         } from './modules/checkm2'
include { MAG_QUALITY     } from './modules/mag_quality'

workflow {

    // ---- Parameter checks ----
    if (!params.input) {
        error "Please provide a samplesheet with --input"
    }
    if (!(params.assembly_mode in ['coassembly', 'per_sample'])) {
        error "--assembly_mode must be 'coassembly' or 'per_sample' (got '${params.assembly_mode}')"
    }

    def binners = params.binners.toString().tokenize(',').collect { it.trim() }
    def unknown_binners = binners - ['concoct', 'metabat2']
    if (!binners || unknown_binners) {
        error "--binners must be a comma-separated list of: concoct, metabat2 (got '${params.binners}')"
    }
    if (params.refine && binners.size() < 2) {
        log.warn "--refine with a single binner: DAS_Tool will only filter bins, not combine binners"
    }
    if (!params.skip_checkm2 && !params.checkm2_db && !workflow.stubRun) {
        error "Please provide the CheckM2 database with --checkm2_db, or skip MAG quality assessment with --skip_checkm2"
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

    // ---- Binning ----
    bins_ch = channel.empty()

    if ('metabat2' in binners) {
        METABAT2(FILTER_CONTIGS.out.contigs.join(CONTIG_DEPTHS.out.depth))
        bins_ch = bins_ch.mix(METABAT2.out.bins)
    }

    if ('concoct' in binners) {
        CONCOCT(FILTER_CONTIGS.out.contigs.join(bams_per_assembly))
        bins_ch = bins_ch.mix(CONCOCT.out.bins)
    }

    // ---- Bin refinement (optional) ----
    if (params.refine) {
        dastool_input = FILTER_CONTIGS.out.contigs
            .join(bins_ch.map { id, binner, dir -> tuple(id, dir) }.groupTuple())

        DAS_TOOL(dastool_input)
        final_bins = DAS_TOOL.out.bins
        all_bins   = bins_ch.mix(DAS_TOOL.out.bins)
    } else {
        final_bins = bins_ch
        all_bins   = bins_ch
    }

    BIN_SUMMARY(all_bins.map { id, binner, dir -> dir }.collect())

    // ---- MAG quality ----
    if (!params.skip_checkm2) {
        checkm2_db = params.checkm2_db
            ? file(params.checkm2_db, checkIfExists: !workflow.stubRun)
            : file("${projectDir}/assets/NO_FILE")

        CHECKM2(final_bins, checkm2_db)
        MAG_QUALITY(
            BIN_SUMMARY.out.tsv,
            CHECKM2.out.report.map { id, binner, report -> report }.collect()
        )
    }
}
