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
| MAG quality and final MAG table | CheckM2 | ✅ |

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

A CSV samplesheet with three columns:

```csv
sample,fastq_1,fastq_2
sampleA,/path/sampleA_R1.fastq.gz,/path/sampleA_R2.fastq.gz
sampleB,/path/sampleB_R1.fastq.gz,/path/sampleB_R2.fastq.gz
```

## Main parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | – | Samplesheet (CSV) |
| `--outdir` | `results` | Output directory |
| `--assembly_mode` | `coassembly` | `coassembly` or `per_sample` |
| `--min_read_length` | `50` | Minimum read length after trimming |
| `--min_contig_length` | `1500` | Minimum contig length kept for mapping and binning |
| `--binners` | `concoct,metabat2` | Binners to run: `concoct`, `metabat2` or both |
| `--refine` | `false` | Combine and refine bins with DAS_Tool |
| `--dastool_score_threshold` | `0.5` | Minimum DAS_Tool score for a bin to be kept |
| `--checkm2_db` | – | CheckM2 DIAMOND database (`uniref100.KO.1.dmnd`) |
| `--skip_checkm2` | `false` | Skip MAG quality assessment |

## Output

```
results/
├── 01_qc/          # fastp HTML and JSON reports
├── 02_assembly/    # filtered contigs, assembly stats (n contigs, total length, N50)
│   └── raw/        # unfiltered MEGAHIT contigs and logs
├── 03_mapping/     # contig depth table, mapping summary (alignment rate per sample)
│   └── logs/       # Bowtie2 logs
├── 04_binning/     # bin_summary.tsv (contigs, size, N50, GC per bin)
│   ├── concoct/    # one FASTA per bin
│   ├── metabat2/   # one FASTA per bin
│   └── dastool/    # refined bins (with --refine)
└── 05_quality/     # mag_quality.tsv: final table with completeness, contamination and quality
    └── checkm2/    # CheckM2 reports
```

## MAG quality categories

`mag_quality.tsv` classifies each bin from CheckM2 completeness and contamination, following the MIMAG thresholds:

| Category | Completeness | Contamination |
|----------|--------------|---------------|
| high | ≥ 90 % | < 5 % |
| medium | ≥ 50 % | < 10 % |
| low | anything below medium | |

rRNA and tRNA presence, also part of the MIMAG high-quality standard, is not assessed.

When `--refine` is used, CheckM2 evaluates the DAS_Tool bins; otherwise it evaluates the bins of each binner.

## License

MIT
