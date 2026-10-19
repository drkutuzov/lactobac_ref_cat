# Lactobacillus Reference Genome Catalogue

A reproducible Snakemake pipeline for downloading, quality-controlling,
taxonomically validating, and curating bacterial reference genomes.

The pipeline generates a versioned Lactobacillus reference catalogue
from NCBI RefSeq assemblies.

## Workflow

The pipeline consists of six stages:

1. **Acquisition** — Download reference genome assemblies from NCBI
   using NCBI Datasets.
2. **Assembly QC** — Calculate assembly quality statistics using QUAST.
3. **Taxonomic validation** — Classify genome sequences using Kraken2.
4. **Genome quality assessment** — Estimate completeness and
   contamination using CheckM2.
5. **Curation** — Apply quality and taxonomy thresholds to select genomes.
6. **Catalogue statistics** — Summarize the resulting FASTA catalogue.

Snakemake manages dependencies and executes the workflow in the
appropriate order.

## Dataset

The dataset contains **50 NCBI RefSeq assemblies**:

- 45 *Lactobacillus* reference genomes.
- 5 negative-control genomes from outside the expected genus.

Assembly accessions are specified in `config.yaml`.

The negative controls test whether the taxonomic validation step correctly rejects genomes outside the expected genus, including high-quality assemblies that pass assembly quality checks.

## Requirements

- Conda or Miniforge
- Snakemake 9
- Internet access for downloading genome assemblies
- Local storage for Kraken2 and CheckM2 reference databases

The workflow uses separate Conda environments for NCBI Datasets,
QUAST, Kraken2, and CheckM2.

Software dependencies are specified in `workflow/envs/`.

## Reference databases

Reference databases must be installed separately before running
the complete workflow.

### Kraken2

The pipeline supports two Kraken2 database configurations:

| Database | Intended use | Release |
|---|---|---|
| Standard-8 (8 GB) | Smaller-scale local testing | 2026-06-26 |
| Full Standard | Complete catalogue production and HPC execution | 2026-06-26 |

The full 50-genome catalogue was generated using the **full Kraken2 Standard database**.

#### Database configuration

Specify the database directory in `config.yaml`:

```yaml
databases:
  kraken2: "databases/kraken2"
```

For the full Standard database, use:

```yaml
databases:
  kraken2: "databases/kraken2_standard"
```

On HPC systems, the database path can be overridden in `workflow/config/hpc.yaml`.

The selected database must contain `hash.k2d`, `opts.k2d`, `taxo.k2d`, `nodes.dmp`, and `names.dmp`. The last two files are required by the pipeline for resolving taxonomic ancestry.

#### Database provenance

**Standard-8 database — local testing**

- Release: 2026-06-26
- Archive: `k2_standard_08_GB_20260626.tar.gz`
- SHA-256:

```text
b17d05ca1459564b49b63d014c4b2ee6ebe1ca0143e6c184e0dfd6d940a55981
```

**Full Standard database — HPC validation**

- Release: 2026-06-26
- Archive: `k2_standard_20260626.tar.gz`
- Download URL: `https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20260626.tar.gz`
- Extracted Kraken2 index (`hash.k2d`): approximately 104 GB
- Used for the successful 50-genome Hendrix SLURM run on 19 October 2026

A SHA-256 checksum for the full Standard archive has not yet been recorded.

Both databases are stored outside version control because of their size.

### Length-weighted taxonomy validation

Kraken2 is used to classify each genome assembly.

In addition to identifying the dominant genus, the pipeline calculates
the percentage of assembly base pairs belonging to contigs classified
within the expected genus.

For each contig, its Kraken2 taxonomic assignment is traced to the
genus level using the NCBI taxonomy tree.

The percentage is calculated as:

expected_genus_bp_percent =
    100 × expected_genus_bp / total_assembly_bp

An assembly passes taxonomy validation only if:

1. Its dominant genus matches the expected genus.
2. At least 95% of its assembly base pairs belong to contigs
   classified within the expected genus.

Both the expected genus and the threshold are configurable
in config.yaml.

The calculation also reports the percentage of assembly base pairs
belonging to unclassified contigs.

**Limitation:** This is a contig-length-weighted metric, not
a nucleotide-level classification. A contig assigned to a genus
contributes its entire length to that genus.

**Database requirement:** The configured Kraken2 database must
contain nodes.dmp and names.dmp for taxonomic ancestry resolution.


### CheckM2

CheckM2 requires a separate DIAMOND reference database.

Activate the CheckM2 environment and download the database:

```bash
checkm2 database --download
```

Alternatively, register an existing database:

```bash
checkm2 database --setdblocation /path/to/uniref100.KO.1.dmnd
```

Verify the configured location:

```bash
checkm2 database --current
```

CheckM2 also supports specifying the database through the `CHECKM2DB`
environment variable.

The database is not included in this repository.

## Configuration

Pipeline parameters are defined in `config.yaml`.

The curation criteria are:

| Criterion | Threshold |
|---|---|
| Assembly N50 | ≥ 50,000 bp |
| CheckM2 completeness | ≥ 90% |
| CheckM2 contamination | ≤ 5% |
| Expected dominant genus | Lactobacillus |

Genomes failing any criterion are excluded from the final catalogue.

## Running the pipeline

From the project directory, activate your Snakemake environment.

Run the complete workflow:

```bash
snakemake --cores 4 --software-deployment-method conda
```

Snakemake creates the required software environments automatically.

To preview the planned jobs without executing them:

```bash
snakemake -n --cores 4 --software-deployment-method conda
```

## HPC execution (SLURM)

The complete 50-genome workflow was successfully executed on the University of Copenhagen's **Hendrix HPC cluster** using Snakemake 9 and the SLURM executor.

The cluster configuration is defined in:

- `workflow/profiles/slurm/config.yaml` — SLURM executor, job limits, resource requests, and account/partition settings.
- `workflow/config/hpc.yaml` — HPC-specific database paths.

Before execution, ensure that Snakemake and Conda are available and that the Kraken2 and CheckM2 databases have been installed.

From the project root, launch the workflow with:

```bash
snakemake --profile workflow/profiles/slurm \
  --resources taxonomy_slots=2 \
  --rerun-triggers mtime input
```

This configuration permits up to 20 concurrently submitted jobs, including a maximum of two simultaneous Kraken2 taxonomy jobs. Each Kraken2 job requests 128 GB of RAM to accommodate the full Standard database.

### HPC validation results

The complete workflow finished successfully on **19 October 2026**:

- **202 of 202 Snakemake steps completed.**
- All 50 candidate genomes were processed.
- 45 *Lactobacillus* genomes passed the quality and taxonomy criteria.
- Five negative-control genomes were excluded.
- The final FASTA catalogue and accompanying statistics were generated.
- The automated `validate_catalogue` rule completed successfully.

This run demonstrates successful end-to-end execution of the pipeline on a SLURM-managed HPC cluster.

## Outputs

The main outputs are:

| File | Description |
|---|---|
| `results/lactobacillus_catalogue_v0.1.fasta` | Curated reference sequences |
| `results/catalogue_v0.1_report.tsv` | Per-genome QC, taxonomy, and curation decisions |
| `results/catalogue_v0.1_stats.tsv` | Summary statistics for the curated FASTA |

Additional intermediate outputs include:

- `results/qc_summary.tsv`
- `results/taxonomy_summary.tsv`
- `results/checkm2/quality_report.tsv`

FASTA headers are prefixed with their source NCBI assembly accession
to preserve provenance.

## Example results

Catalogue version 0.1 was generated from 50 candidate assemblies on the Hendrix HPC cluster on 19 October 2026.

| Metric | Result |
|---|---:|
| Candidate assemblies | 50 |
| Accepted assemblies | 45 |
| Rejected assemblies | 5 |
| Acceptance rate | 90% |
| Automated catalogue validation | PASS |
| Snakemake execution | 202/202 steps completed |

All 45 intended *Lactobacillus* genomes were retained, while all five negative controls were excluded.

The curated FASTA is available at `results/lactobacillus_catalogue_v0.1.fasta`, and individual inclusion or exclusion decisions are documented in `results/catalogue_v0.1_report.tsv`.

## Reproducibility

The pipeline uses:

- Snakemake for workflow orchestration.
- Conda environments for software dependencies.
- NCBI RefSeq assembly accessions for input provenance.
- Configurable quality thresholds.
- Versioned catalogue outputs.

Large reference databases and downloaded genome assemblies are excluded
from Git and must be obtained separately.