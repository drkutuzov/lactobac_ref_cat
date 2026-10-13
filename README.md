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

The example dataset contains eight NCBI RefSeq assemblies:

- Six Lactobacillus genomes.
- Two negative controls: Escherichia coli and Bacillus subtilis.

Assembly accessions are specified in `config.yaml`.

The negative controls test whether the taxonomic validation step
correctly rejects genomes outside the expected genus.

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

The pipeline expects a Kraken2 database at:

`databases/kraken2/`

This location is configurable in `config.yaml`:

```yaml
databases:
  kraken2: "databases/kraken2"
```

The example workflow was tested with the Kraken2 Standard-8 database
dated 2026-06-26.

Download and extract a compatible Kraken2 database into the configured
directory.

The database directory is excluded from Git because of its size.

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

Catalogue version 0.1 produced:

| Metric | Result |
|---|---:|
| Candidate assemblies | 8 |
| Accepted assemblies | 6 |
| Rejected assemblies | 2 |
| FASTA records | 8 |
| Total sequence length | 11,915,429 bp |

Both negative controls were rejected by the taxonomic validation gate.

## Reproducibility

The pipeline uses:

- Snakemake for workflow orchestration.
- Conda environments for software dependencies.
- NCBI RefSeq assembly accessions for input provenance.
- Configurable quality thresholds.
- Versioned catalogue outputs.

Large reference databases and downloaded genome assemblies are excluded
from Git and must be obtained separately.