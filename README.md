# mag-pipeline

Nextflow DSL2 pipeline for recovering metagenome-assembled genomes (MAGs) from paired-end shotgun metagenomes.

> 🚧 Work in progress — QC, assembly and read mapping are implemented; binning and MAG quality assessment are next.

## Workflow

| Step | Tool | Status |
|------|------|--------|
| Read QC and trimming | fastp | ✅ |
| Assembly (co-assembly or per sample) | MEGAHIT | ✅ |
| Contig length filtering and assembly stats | custom awk script | ✅ |
| Read mapping and contig coverage | Bowtie2, samtools, jgi_summarize_bam_contig_depths | ✅ |
| Binning | CONCOCT, MetaBAT2 | 🔜 |
| Bin refinement | DAS_Tool | 🔜 |
| MAG quality | CheckM2 | 🔜 |

## Quick start

```bash
# Test run with small public data
nextflow run . -profile test,singularity

# Your own data on a SLURM cluster
nextflow run . -profile singularity,slurm \
    --input samplesheet.csv \
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

## Output

```
results/
├── 01_qc/          # fastp HTML and JSON reports
├── 02_assembly/    # filtered contigs, assembly stats (n contigs, total length, N50)
│   └── raw/        # unfiltered MEGAHIT contigs and logs
└── 03_mapping/     # contig depth table, mapping summary (alignment rate per sample)
    └── logs/       # Bowtie2 logs
```

## License

MIT
