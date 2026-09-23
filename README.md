# mag-pipeline

Nextflow DSL2 pipeline for recovering metagenome-assembled genomes (MAGs) from paired-end shotgun metagenomes.

From raw paired-end reads to quality-assessed MAGs. Every step has been tested on an HPC cluster (SLURM) with conda environments and in CI with stub runs.

## Workflow

| Step | Tool | Status |
|------|------|--------|
| Read QC and trimming | fastp | ✅ |
| Assembly (co-assembly or per sample) | MEGAHIT | ✅ |
| Contig length filtering and assembly stats | custom awk script | ✅ |
| Read mapping and contig coverage | Bowtie2, samtools, jgi_summarize_bam_contig_depths | ✅ |
| Binning and bin statistics | CONCOCT, MetaBAT2 | ✅ |
| Bin refinement (optional) | DAS_Tool | ✅ |
| MAG quality | CheckM2 | ✅ |
| Taxonomic classification (optional) | GTDB-Tk | ✅ |
| Functional annotation (optional) | Prodigal, eggNOG-mapper | ✅ |

## Designed for constrained resources

Developed and tested on a shared SLURM cluster with a 1 TB per-user quota and no
administrator rights. Peak disk usage is a design constraint, not an afterthought:

| Choice | Effect |
|--------|--------|
| Alignments are never kept as persistent BAMs | Each sample is mapped, its depth and coverage computed, and the BAM deleted within the same task, so peak usage scales with `--max_mapping_jobs`, not with the number of samples |
| MEGAHIT intermediate files are removed after each assembly | Avoids keeping k-mer graphs for the whole run |
| `--skip_fastp` starts from trimmed reads | Lets you delete the raw FASTQ files once QC is done |
| Conda environments installed at user level | No root access required |

On 48 samples with 8 co-assemblies and all-vs-all mapping (384 alignments), peak
usage stays under XXX GB.

## Quick start

```bash
# Test run with small public data
nextflow run . -profile test,singularity

# Your own data on a SLURM cluster
nextflow run . -profile singularity,slurm \
    --input samplesheet.csv \
    --checkm2_db /path/to/CheckM2_database/uniref100.KO.1.dmnd \
    --refine \
    --outdir results
```

### Profiles

| Profile | Use |
|---------|-----|
| `test` | Small public test dataset |
| `docker` / `singularity` | Run tools from containers |
| `slurm` | Submit each process as a SLURM job |
| `pyrene` | Use pre-installed conda environments on the Pyrene cluster (UPPA) |

## Input

A CSV samplesheet with three columns, plus an optional `group` column used by
`--assembly_mode group`:

```csv
sample,fastq_1,fastq_2
sampleA,/path/sampleA_R1.fastq.gz,/path/sampleA_R2.fastq.gz
sampleB,/path/sampleB_R1.fastq.gz,/path/sampleB_R2.fastq.gz
```

## Choosing an assembly strategy

Pooling all samples, pooling them by group, or assembling each sample alone
changes which genomes are recovered and whether they are redundant across
assemblies. [docs/assembly_strategies.md](docs/assembly_strategies.md) explains
the trade-offs and when to dereplicate. Read it before a large run.

## Main parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | – | Samplesheet (CSV) |
| `--outdir` | `results` | Output directory |
| `--assembly_mode` | `coassembly` | `coassembly`, `group` or `per_sample`; see [assembly strategies](docs/assembly_strategies.md) |
| `--map_all_samples` | `false` | Map every sample against every assembly |
| `--min_read_length` | `50` | Minimum read length after trimming |
| `--min_contig_length` | `1500` | Minimum contig length kept for mapping and binning |
| `--binners` | `concoct,metabat2` | Binners to run: `concoct`, `metabat2` or both |
| `--refine` | `false` | Combine and refine bins with DAS_Tool |
| `--dastool_score_threshold` | `0.5` | Minimum DAS_Tool score for a bin to be kept |
| `--checkm2_db` | – | CheckM2 DIAMOND database (`uniref100.KO.1.dmnd`) |
| `--skip_checkm2` | `false` | Skip MAG quality assessment |
| `--gtdbtk_db` | – | GTDB-Tk reference data directory (enables taxonomy) |
| `--gtdbtk_pplacer_cpus` | `1` | Threads for pplacer; more threads multiply memory use |
| `--gtdbtk_args` | – | Extra flags for `gtdbtk classify_wf` |
| `--eggnog_db` | – | eggNOG-mapper data directory (enables functional annotation) |
| `--eggnog_dbmem` | `false` | Load the eggNOG database in RAM (~45 GB, much faster) |
| `--eggnog_args` | – | Extra flags for `emapper.py` |

## Output

```
results/
├── 01_qc/          # fastp HTML and JSON reports
├── 02_assembly/    # filtered contigs, assembly stats (n contigs, total length, N50)
│   └── raw/        # unfiltered MEGAHIT contigs and logs
├── 03_mapping/     # contig depth table, mapping summary (alignment rate per sample)
│   └── logs/       # Bowtie2 logs
├── 04_binning/     # bin_summary.tsv, contig2bin.tsv, mag_abundance.tsv
│   ├── concoct/    # one FASTA per bin
│   ├── metabat2/   # one FASTA per bin
│   └── dastool/    # refined bins (with --refine)
├── 05_quality/     # mag_quality.tsv: final table with quality and taxonomy per MAG
│   └── checkm2/    # CheckM2 reports
├── 06_taxonomy/    # GTDB-Tk classification (with --gtdbtk_db)
├── 07_function/    # ko_per_mag.tsv, annotation_stats.tsv (with --eggnog_db)
│   └── genes/      # protein sequences and gene counts per MAG
└── run_info.txt    # pipeline version, command line, parameters and database dates
```

## MAG quality categories

`mag_quality.tsv` classifies each bin from CheckM2 completeness and contamination, following the MIMAG thresholds:

| Category | Completeness | Contamination |
|----------|--------------|---------------|
| high | ≥ 90 % | < 5 % |
| medium | ≥ 50 % | < 10 % |
| low | anything below medium | |

rRNA and tRNA presence, also part of the MIMAG high-quality standard, is not assessed.

When `--refine` is used, CheckM2 and GTDB-Tk evaluate the DAS_Tool bins; otherwise they evaluate the bins of each binner.

## Functional annotation

With `--eggnog_db`, genes are called on the final bins with Prodigal (`-p single`,
since bins are genomes rather than raw metagenomic contigs) and annotated with
eggNOG-mapper. Protein names carry their bin (`<bin>|<gene>`), so the pipeline
can write two tidy tables: `ko_per_mag.tsv`, with the number of genes per MAG and
KEGG orthologue, and `annotation_stats.tsv`, with genes, annotated genes and
distinct KOs per MAG. eggNOG-mapper only reports the genes it could annotate, so
the total gene count comes from Prodigal: the annotated fraction is a real
number, not 100% by construction.

## Traceability

Every run writes `run_info.txt` next to the results: pipeline and Nextflow
versions, the exact command line, every parameter and the modification date of
each reference database. eggNOG and GTDB releases are not versioned by name, so
that date is what makes an annotation reproducible months later.

## Taxonomy

Taxonomic classification is optional because the GTDB-Tk reference data is ~100 GB and `pplacer` needs more than 100 GB of RAM. Pass `--gtdbtk_db /path/to/gtdb/release` to enable it; the `process_high_memory` label lets you send only this step to a large-memory partition. The resulting ranks are added to `mag_quality.tsv`, which reports `NA` when taxonomy is not run.

## License

MIT
