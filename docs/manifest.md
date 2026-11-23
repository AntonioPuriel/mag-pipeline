# manifest.json

At the end of every run, successful or not, mag-pipeline writes
`manifest.json` in the results directory, next to `run_info.txt`. It records
what produced the results and lists the summary tables that exist, so that
[magreport](https://github.com/AntonioPuriel/magreport) (or any other program)
does not have to guess file names or folder layouts.

## Fields

| Field | Content |
|---|---|
| `schema_version` | Version of this format, `"major.minor"` (see below) |
| `pipeline` | `name`, `version`, `revision` and `commit` (both `null` when the pipeline is run from a local folder rather than with `nextflow run AntonioPuriel/mag-pipeline -r ...`), `nextflow` version |
| `run` | `name`, `session_id`, `started`, `completed`, `duration`, `success`, `exit_status`, `profile`, `container_engine`, `command_line` |
| `databases` | One entry per database given (`checkm2`, `gtdbtk`, `eggnog`, `canthyd`, `hadeg`, `cat`, `singlem`, `gunc`): `path` and `modified` (date of last modification, `null` if not readable from the head node) |
| `params` | Every parameter of the run |
| `tables` | Table name → list of paths, relative to `manifest.json` |

Only the tables that exist are listed: which ones depend on `--workflow`,
`--refine`, `--gunc_db`, `--mimag`, `--gtdbtk_db`, `--mag_tree`, `--eggnog_db`,
`--canthyd_db`/`--hadeg_db` and `--singlem_db`. A table is always a list, even
when there is one file, because some tables have one file per assembly or per
sample.

| Table | Files |
|---|---|
| `assembly_stats` | `02_assembly/<id>.contigs.stats.tsv` |
| `mapping_summary` | `03_mapping/mapping_summary.tsv` |
| `depth` | `03_mapping/<id>.depth.txt` |
| `bin_summary` | `04_binning/bin_summary.tsv` |
| `contig2bin` | `04_binning/contig2bin.tsv` |
| `mag_abundance` | `04_binning/mag_abundance.tsv` |
| `dastool_summary` | `04_binning/dastool/<id>_DASTool_summary.tsv` |
| `mag_quality` | `05_quality/mag_quality.tsv` |
| `mimag_quality` | `05_quality/mimag_quality.tsv` |
| `checkm2` | `05_quality/checkm2/<id>_<set>.checkm2.tsv` |
| `gunc` | `05_quality/gunc/<id>_<set>.gunc.tsv` |
| `dereplication` | `05_quality/dereplication/dereplication_clusters.tsv` |
| `gtdbtk` | `06_taxonomy/<id>_<set>.gtdbtk.tsv` |
| `mag_tree` | `06_taxonomy/mag_tree/*.tree` |
| `ko_per_mag` | `07_function/ko_per_mag.tsv` |
| `annotation_stats` | `07_function/annotation_stats.tsv` |
| `hydrocarbon_per_mag` | `07_function/hydrocarbon_per_mag.tsv` |
| `gene_catalog_stats` | `08_contigs/catalog/gene_catalog.stats.tsv` |
| `gene_abundance_stats` | `08_contigs/abundance/gene_abundance_stats.tsv` |
| `profile_stats` | `08_contigs/profiles/profile_stats.tsv` |
| `taxonomy_per_sample` | `08_contigs/profiles/taxonomy_per_sample.tsv` |
| `ko_per_sample` | `08_contigs/profiles/ko_per_sample.tsv` |
| `canthyd_per_sample` | `08_contigs/profiles/canthyd_per_sample.tsv` |
| `hadeg_per_sample` | `08_contigs/profiles/hadeg_per_sample.tsv` |
| `singlem` | `09_singlem/singlem_profiles_with_extras.tsv` |
| `singlem_per_sample` | `09_singlem/profiles/<sample>.singlem.tsv` |

Any of these may also end in `.gz`.

## The report/ folder

Next to the manifest, every run writes `report/`: a copy of the summary tables
small enough to copy to a laptop and read with `magreport::read_pipeline()`.

- Same layout as the results (`report/05_quality/mag_quality.tsv.gz`, ...).
- `.tsv` and `.txt` files are gzipped; trees are copied as they are.
- The tables with one row per contig, `depth` and `contig2bin`, are left out:
  they reach several GB in large runs and the report does not use them.
- `report/manifest.json` is the manifest of the run with `tables` pointing to
  the copies, plus `bundle.excluded`, the tables left out.
- `report/run_info.txt` is a copy of `run_info.txt`.

`report/` is rebuilt at the end of every run, also with `-resume`. In the
Marine Honour case study (42 samples) the summary tables take about 100 MB
uncompressed, against 2.6 GB with `depth` and `contig2bin`.

## Versions of the format

- **Minor** (`1.0` → `1.1`): a table or a column is added. Readers written for
  `1.0` keep working.
- **Major** (`1.x` → `2.0`): a table or column is renamed or removed, or a path
  changes meaning. Readers should refuse a major version they do not know.

| schema_version | mag-pipeline | magreport |
|---|---|---|
| (none: `run_info.txt` only) | ≤ 2.0.0 | all (finds tables by name) |
| 1.0 | ≥ 2.1.0 | ≥ 0.9.0 |

## Example

```json
{
    "schema_version": "1.0",
    "pipeline": {
        "name": "AntonioPuriel/mag-pipeline",
        "version": "2.1.0",
        "revision": "v2.1.0",
        "commit": "…",
        "nextflow": "26.04.0"
    },
    "run": {
        "name": "nostalgic_curie",
        "success": true,
        "profile": "pyrene",
        "command_line": "nextflow run AntonioPuriel/mag-pipeline -r v2.1.0 -profile pyrene ..."
    },
    "databases": {
        "gtdbtk": { "path": "/path/to/gtdbtk/db", "modified": "2026-05-02" }
    },
    "params": { "workflow": "both", "binners": "concoct,metabat2" },
    "tables": {
        "mag_quality": ["05_quality/mag_quality.tsv"],
        "gunc": ["05_quality/gunc/coassembly_dastool.gunc.tsv"]
    }
}
```
