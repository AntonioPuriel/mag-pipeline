# Choosing an assembly strategy

How reads are pooled before assembly is the decision that most changes what the
pipeline recovers. There is no best choice in general: it depends on the study
design and on the question. This page explains the three options and what each
one implies downstream.

## The three strategies

### `--assembly_mode coassembly` (default)

All samples are assembled together into a single assembly, and every sample is
mapped back to it.

- **Good for:** a small number of samples from similar communities, where a
  genome that is rare in each sample is common in the pool.
- **Recovers:** more complete genomes of low-abundance organisms, because
  pooling adds up their coverage.
- **Costs:** memory and time grow with the total amount of data, and strains of
  the same species from different samples can collapse into one fragmented
  consensus.
- **Redundancy:** none. There is one assembly, so each genome appears once.

### `--assembly_mode group`

Samples are assembled by group, one co-assembly per value of the `group` column
of the samplesheet (for example one per sampling time, site or treatment).

```csv
sample,fastq_1,fastq_2,group
T0-oil-1,T0-oil-1_R1.fastq.gz,T0-oil-1_R2.fastq.gz,T0
T0-ctrl-1,T0-ctrl-1_R1.fastq.gz,T0-ctrl-1_R2.fastq.gz,T0
T1-oil-1,T1-oil-1_R1.fastq.gz,T1-oil-1_R2.fastq.gz,T1
```

- **Good for:** time series or several sites, where communities change between
  groups but are similar within each one.
- **Recovers:** genomes specific to each group, with a manageable assembly size.
- **Costs:** several medium-sized assemblies instead of one large one; they run
  in parallel.
- **Redundancy:** yes. See below.

### `--assembly_mode per_sample`

Each sample is assembled on its own.

- **Good for:** very different communities, or when strain resolution matters
  more than completeness.
- **Recovers:** sample-specific strains without mixing them.
- **Costs:** the lowest memory per assembly, but many assemblies, and genomes
  that are rare in every single sample may not assemble at all.
- **Redundancy:** yes, and usually the most.

## Mapping: `--map_all_samples`

Binning separates genomes using how the coverage of each contig varies across
samples. The more samples are mapped against an assembly, the more of that
signal the binners get.

- With `coassembly`, every sample is always mapped: the option has no effect.
- With `group` or `per_sample`, by default each sample is only mapped against its
  own assembly. With `--map_all_samples`, every sample is mapped against every
  assembly.

Mapping all samples usually gives cleaner bins, at the cost of more mapping jobs:
with 8 groups and 48 samples, 384 alignments instead of 48. Each one is
independent, so on a cluster they run in parallel.

## Why several assemblies need dereplication

When a genome is present in several groups, it is recovered once per assembly.
An organism that persists across eight sampling times can appear as eight MAGs,
one from each co-assembly. Counting them as eight genomes would inflate
diversity, and their abundances would be split across copies.

Dereplication groups MAGs that belong to the same species, usually at 95%
average nucleotide identity, and keeps one representative per group, chosen by
completeness and contamination. The pipeline warns when a run builds more than
one assembly for this reason.

Until dereplication is part of the pipeline, run it on the final bins with a tool
such as dRep or galah, using the CheckM2 results to choose representatives.

## Quick guide

| Study design | Suggested mode | `--map_all_samples` | Dereplicate |
|--------------|----------------|---------------------|-------------|
| Few samples, one environment | `coassembly` | not needed | no |
| Time series or several sites | `group` | yes | yes |
| Very different communities | `per_sample` | optional | yes |
| Strain-level questions | `per_sample` | yes | yes, at a stricter threshold |
