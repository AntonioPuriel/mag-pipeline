# Contig branch: the community without genomes

MAGs only describe what could be reconstructed as a genome. Rare organisms,
very diverse populations and whole groups (in the Marine Honour case study,
all archaea but one) stay in the contigs or in the reads. `--workflow`
chooses what happens after the assemblies are filtered:

| `--workflow` | What runs |
|--------------|-----------|
| `mags` (default) | Binning, quality, dereplication, taxonomy and annotation of genomes, as in v1.x |
| `contigs` | A gene catalog of every filtered contig, with function, taxonomy and abundance per sample; no binning |
| `both` | Both branches on the same assemblies and the same depth tables |

```
reads → fastp → MEGAHIT → contig filter → mapping of every sample, depth per contig
                                  │
         ┌────────────────────────┴────────────────────────┐
   --workflow mags                                  --workflow contigs
   binning → MAGs (v1.x)                            Prodigal -p meta, every contig
                                                    → gene catalog (MMseqs2 linclust)
                                                    → eggNOG, CANT-HYD, HADEG (in chunks)
                                                    → CAT: taxonomy of every contig
                                                    → depth and TPM of every gene
                                                    → KO / marker / taxon x sample tables

   --singlem_db, with any --workflow: SingleM profile of every sample, from its reads
```

## Steps

**Genes.** Prodigal in metagenomic mode (`-p meta`) on the filtered contigs of
each assembly, split across CPUs. Genes are named `<assembly>|<contig>_<n>`,
so the assembly and contig of every gene can be read from its name.

**Gene catalog.** Genes of at least `--catalog_min_gene_length` nt from every
assembly are clustered with MMseqs2 `easy-linclust` on their nucleotide
sequence: ≥ `--catalog_min_id` identity (0.95) with the representative
covering ≥ `--catalog_min_cov` (0.9) of each member, the representative being
the longest gene (`--cluster-mode 2 --cov-mode 1`). These are the criteria of
reference gene catalogs. linclust scales linearly in time and memory, so
several million genes fit in a standard node; MMseqs2 is told to stay within
80 % of the task memory.

**Abundance, without mapping again.** Every gene takes the mean depth of its
contig from the depth tables used for binning. For a catalog gene:

1. within one assembly, the depths of its members are **summed** (the reads
   of a gene with two variants in the same assembly are split between them);
2. across assemblies, those sums are **averaged**. With `--map_all_samples`
   every sample is mapped against every assembly, so summing across
   assemblies would count the same reads once per assembly that contains the
   gene, and favour genes assembled many times. Only assemblies with a depth
   column for that sample count.

TPM is the depth of each gene over the summed depth of all genes in the
sample, × 10^6. Depth is already reads × read length / contig length, so no
further length correction is needed.

**Function.** The catalog proteins are split in pieces of
`--catalog_chunk_size` and annotated with the same eggNOG-mapper, CANT-HYD and
HADEG steps used for the MAGs (`--eggnog_db`, `--canthyd_db`, `--hadeg_db`).

**Taxonomy.** With `--cat_db`, [CAT](https://github.com/MGXlab/CAT_pack)
classifies each contig by a bit-score-weighted vote of the DIAMOND hits of its
genes, reusing the Prodigal proteins. GTDB-Tk cannot do this: it needs a
genome with its marker genes. Prepare the database from GTDB so that contigs
and MAGs share one taxonomy (see `conf/pyrene.config`). `--cat_block_size` and
`--cat_index_chunks` trade speed for memory in DIAMOND. Every gene takes the
classification of its contig.

**SingleM.** With `--singlem_db`, every sample's reads are profiled with
[SingleM](https://github.com/wwood/singlem). It needs neither assembly nor
binning, so it sees what was never assembled; SingleM 0.21.3 uses GTDB R232
by default. `09_singlem/singlem_domain.tsv` gives the bacterial and archaeal
fraction of each sample directly.

## Output

```
08_contigs/
├── genes/       <assembly>.genes.{faa,fna}.gz, gene counts per assembly
├── catalog/     gene_catalog.{faa,fna}.gz, gene_catalog.clusters.tsv.gz
│                (representative<TAB>member), gene_catalog.stats.tsv
├── abundance/   gene_depth.tsv.gz, gene_tpm.tsv.gz (gene x sample),
│                gene_abundance_stats.tsv
├── function/    catalog_chunkNNN.{emapper.annotations,canthyd.tsv,hadeg.tsv}
├── taxonomy/    <assembly>.CAT.named.txt
└── profiles/    ko_per_sample.tsv, canthyd_per_sample.tsv, hadeg_per_sample.tsv,
                 taxonomy_per_sample.tsv (rank, taxon, sample, n_genes, tpm),
                 contig_taxonomy.tsv, catalog_annotation.tsv.gz (one row per gene),
                 profile_stats.tsv (fraction of TPM annotated by each source)
09_singlem/
├── profiles/    <sample>.singlem.tsv
└── singlem_<rank>.tsv, singlem_profiles_with_extras.tsv
```

Profile tables are long: one row per feature and sample, with the number of
catalog genes behind the feature and the summed TPM. A gene with several KOs
counts towards each; genes on unclassified contigs are `Unclassified` at that
rank.

## Adding the branch to a finished MAG run

`--workflow both -resume` reuses everything up to the depth tables and every
MAG step: the MAG processes are unchanged in v2.0.0. Running `--workflow
contigs` alone skips CONCOCT's coverage, so its mapping tasks differ from
those of a MAG run and are repeated; use `both` to share them.
