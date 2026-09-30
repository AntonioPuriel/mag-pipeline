#!/usr/bin/env nextflow
/*
 * mag-pipeline: MAG recovery from paired-end shotgun metagenomes
 * Author: Antonio Puriel Hernández
 */

include { FASTP           } from './modules/fastp'
include { MEGAHIT         } from './modules/megahit'
include { FILTER_CONTIGS  } from './modules/filter_contigs'
include { BOWTIE2_BUILD   } from './modules/bowtie2_build'
include { MAP_DEPTH; MERGE_DEPTHS; MERGE_CONCOCT_COV; CONCOCT_CUTUP } from './modules/map_depth'
include { MAPPING_SUMMARY } from './modules/mapping_summary'
include { METABAT2        } from './modules/metabat2'
include { CONCOCT         } from './modules/concoct'
include { BIN_SUMMARY     } from './modules/bin_summary'
include { DAS_TOOL        } from './modules/dastool'
include { CHECKM2         } from './modules/checkm2'
include { MAG_QUALITY     } from './modules/mag_quality'
include { GTDB_TK         } from './modules/gtdbtk'
include { PRODIGAL           } from './modules/prodigal'
include { EGGNOG_MAPPER      } from './modules/eggnog_mapper'
include { FUNCTIONAL_SUMMARY } from './modules/functional_summary'
include { MAG_ABUNDANCE      } from './modules/mag_abundance'
include { GALAH              } from './modules/galah'

workflow {

    // ---- Parameter checks ----
    if (!params.input) {
        error "Please provide a samplesheet with --input"
    }
    if (!(params.assembly_mode in ['coassembly', 'group', 'per_sample'])) {
        error "--assembly_mode must be 'coassembly', 'group' or 'per_sample' (got '${params.assembly_mode}'). See docs/assembly_strategies.md"
    }

    if (params.assembly_mode in ['group', 'per_sample'] && !params.dereplicate) {
        log.warn "--assembly_mode ${params.assembly_mode} builds more than one assembly: the same genome can be " +
                 "recovered once per assembly, so MAGs may be redundant across assemblies. " +
                 "Use --dereplicate before counting genomes (see docs/assembly_strategies.md)."
    }
    if (params.dereplicate && params.skip_checkm2) {
        error "--dereplicate needs CheckM2 to choose the best genome of each cluster: remove --skip_checkm2"
    }
    if (params.map_all_samples && params.assembly_mode == 'coassembly') {
        log.info "--map_all_samples has no effect with --assembly_mode coassembly: every sample is already mapped"
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

    // ---- Read samplesheet: sample,fastq_1,fastq_2[,group] ----
    samplesheet = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (params.assembly_mode == 'group' && !row.group) {
                error "--assembly_mode group needs a 'group' column in the samplesheet (sample '${row.sample}' has none)"
            }
            tuple(row.sample, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true), row.group ?: 'all')
        }

    reads_ch      = samplesheet.map { sample, r1, r2, group -> tuple(sample, r1, r2) }
    sample_groups = samplesheet.map { sample, r1, r2, group -> tuple(sample, group) }

    // ---- Quality control ----
    // --skip_fastp: the samplesheet already points to trimmed reads
    if (params.skip_fastp) {
        trimmed_reads = reads_ch
    } else {
        FASTP(reads_ch)
        trimmed_reads = FASTP.out.reads
    }

    // Trimmed reads with the group of each sample: sample, r1, r2, group
    trimmed = trimmed_reads.join(sample_groups)

    // ---- Assembly input (see docs/assembly_strategies.md) ----
    if (params.assembly_mode == 'coassembly') {
        // one assembly with every sample
        assembly_input = trimmed
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { samples ->
                tuple('coassembly', samples.collect { it[1] }, samples.collect { it[2] })
            }
    } else if (params.assembly_mode == 'group') {
        // one co-assembly per value of the group column
        assembly_input = trimmed
            .map { sample, r1, r2, group -> tuple(group, sample, r1, r2) }
            .groupTuple(sort: { a, b -> a.toString() <=> b.toString() })
            .map { group, samples, r1s, r2s -> tuple(group, r1s, r2s) }
    } else {
        // one assembly per sample
        assembly_input = trimmed
            .map { sample, r1, r2, group -> tuple(sample, [r1], [r2]) }
    }

    // ---- Assembly and contig filtering ----
    MEGAHIT(assembly_input)
    FILTER_CONTIGS(MEGAHIT.out.contigs)

    // ---- Read mapping (BAMs are deleted inside MAP_DEPTH) ----
    BOWTIE2_BUILD(FILTER_CONTIGS.out.contigs)

    // CONCOCT needs its 10 kb pieces before mapping, to compute its coverage
    // while each BAM still exists
    if ('concoct' in binners) {
        CONCOCT_CUTUP(FILTER_CONTIGS.out.contigs)
        beds = CONCOCT_CUTUP.out.cutup.map { id, fa, bed -> tuple(id, bed) }
    } else {
        no_bed = file("${projectDir}/assets/NO_FILE")
        beds   = FILTER_CONTIGS.out.contigs.map { id, contigs -> tuple(id, no_bed) }
    }

    // id, index, bed
    index_bed = BOWTIE2_BUILD.out.index.join(beds)

    reads_for_mapping = trimmed.map { sample, r1, r2, group -> tuple(sample, r1, r2) }

    if (params.assembly_mode == 'coassembly' || params.map_all_samples) {
        // every sample is mapped against every assembly: more coverage
        // profiles for binning, at the cost of more mapping jobs
        mapping_input = index_bed
            .combine(reads_for_mapping)
    } else if (params.assembly_mode == 'group') {
        // each sample is mapped against the co-assembly of its own group
        mapping_input = index_bed
            .combine(trimmed.map { sample, r1, r2, group -> tuple(group, sample, r1, r2) }, by: 0)
    } else {
        // each sample is mapped against its own assembly
        mapping_input = index_bed
            .join(reads_for_mapping)
            .map { id, index, bed, r1, r2 -> tuple(id, index, bed, id, r1, r2) }
    }

    // id, index, bed, sample, r1, r2
    MAP_DEPTH(mapping_input)
    MAPPING_SUMMARY(MAP_DEPTH.out.log.collect())

    // Number of alignments expected per assembly. With it, the coverage of an
    // assembly is merged (and its binning starts) as soon as its own alignments
    // finish, instead of waiting for every alignment of every assembly.
    if (params.assembly_mode == 'coassembly' || params.map_all_samples) {
        // every sample against every assembly
        n_expected = FILTER_CONTIGS.out.contigs
            .map { id, contigs -> id }
            .combine(samplesheet.count())
    } else if (params.assembly_mode == 'group') {
        // the samples of each group against their own co-assembly
        n_expected = sample_groups
            .map { sample, group -> tuple(group, sample) }
            .groupTuple()
            .map { group, samples -> tuple(group, samples.size()) }
    } else {
        // each sample against its own assembly
        n_expected = FILTER_CONTIGS.out.contigs.map { id, contigs -> tuple(id, 1) }
    }

    // ---- Contig coverage per assembly (input for binning) ----
    // The files of each assembly are grouped as soon as all have arrived, and
    // sorted by name so that merged columns keep the same order in every run
    MERGE_DEPTHS(
        MAP_DEPTH.out.depth
            .map { id, sample, depth -> tuple(id, depth) }
            .combine(n_expected, by: 0)
            .map { id, depth, n -> tuple(groupKey(id, n), depth) }
            .groupTuple(sort: true)
            .map { key, depths -> tuple(key.getGroupTarget(), depths) }
    )
    depth_ch = MERGE_DEPTHS.out.depth

    // ---- Binning ----
    bins_ch = channel.empty()

    if ('metabat2' in binners) {
        METABAT2(FILTER_CONTIGS.out.contigs.join(depth_ch))
        bins_ch = bins_ch.mix(METABAT2.out.bins)
    }

    if ('concoct' in binners) {
        MERGE_CONCOCT_COV(
            MAP_DEPTH.out.concoct_cov
                .map { id, sample, cov -> tuple(id, cov) }
                .combine(n_expected, by: 0)
                .map { id, cov, n -> tuple(groupKey(id, n), cov) }
                .groupTuple(sort: true)
                .map { key, covs -> tuple(key.getGroupTarget(), covs) }
        )
        // id, contigs, cutup_fa, coverage
        CONCOCT(
            FILTER_CONTIGS.out.contigs
                .join(CONCOCT_CUTUP.out.cutup.map { id, fa, bed -> tuple(id, fa) })
                .join(MERGE_CONCOCT_COV.out.coverage)
        )
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

    // Contig-to-bin map and coverage of every MAG in every sample, the input
    // for abundance figures downstream
    MAG_ABUNDANCE(
        all_bins.map { id, binner, dir -> dir }.collect(),
        depth_ch.map { id, depth -> depth }.collect()
    )

    // Bins passed to taxonomy and annotation: all final bins, or only the
    // dereplicated representatives with --dereplicate
    annotation_bins = final_bins

    // ---- MAG quality ----
    if (!params.skip_checkm2) {
        checkm2_db = params.checkm2_db
            ? file(params.checkm2_db, checkIfExists: !workflow.stubRun)
            : file("${projectDir}/assets/NO_FILE")

        CHECKM2(final_bins, checkm2_db)
        reports = CHECKM2.out.report.map { id, binner, report -> report }

        // ---- Dereplication (optional): one representative per cluster ----
        // Runs once, when the bins and CheckM2 reports of every assembly are ready
        if (params.dereplicate) {
            GALAH(
                final_bins.map { id, binner, dir -> dir }.collect(),
                reports.collect()
            )
            annotation_bins = GALAH.out.representatives
                .map { dir -> tuple('dereplicated', 'galah', dir) }
        }

        // ---- Taxonomy (optional, needs the GTDB-Tk reference data) ----
        if (params.gtdbtk_db) {
            gtdbtk_db = file(params.gtdbtk_db, checkIfExists: !workflow.stubRun)
            GTDB_TK(annotation_bins, gtdbtk_db)
            reports = reports.mix(GTDB_TK.out.summary.map { id, binner, summary -> summary })
        }

        MAG_QUALITY(BIN_SUMMARY.out.tsv, reports.collect())
    }
    else if (params.gtdbtk_db) {
        log.warn "--gtdbtk_db is ignored when --skip_checkm2 is set"
    }

    // ---- Functional annotation (optional, needs the eggNOG database) ----
    if (params.eggnog_db) {
        eggnog_db = file(params.eggnog_db, checkIfExists: !workflow.stubRun)

        PRODIGAL(annotation_bins)
        EGGNOG_MAPPER(PRODIGAL.out.proteins, eggnog_db)
        FUNCTIONAL_SUMMARY(
            EGGNOG_MAPPER.out.annotations.map { id, binner, ann -> ann }.collect(),
            PRODIGAL.out.counts.collect()
        )
    }
}
