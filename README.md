# Lactobacillus Genome Pipeline

A small bioinformatics pipeline for acquiring and analysing bacterial reference genomes.

The pipeline reads NCBI RefSeq assembly accessions from `config.yaml` and downloads the corresponding genomic FASTA files using NCBI Datasets.

## Curated catalogue

The pipeline combines assembly QC and taxonomic validation to produce a
versioned Lactobacillus reference catalogue.

A candidate genome is included only if it:

- passes the configured QUAST N50 threshold
- is classified by Kraken2 with Lactobacillus as the dominant genus

Version 0.1 contains 6 accepted genomes represented by 8 FASTA records.
Two planted negative controls, Escherichia coli and Bacillus subtilis,
are rejected by the taxonomic validation gate.

Outputs:

- `results/lactobacillus_catalogue_v0.1.fasta` — curated reference catalogue
- `results/catalogue_v0.1_report.tsv` — candidate-level curation report with QC,
  taxonomy, final status, and rejection reason

Catalogue FASTA headers are prefixed with the source NCBI assembly accession
to preserve provenance.

## Dataset

The example dataset contains:

- 6 Lactobacillus genomes
- 2 non-Lactobacillus genomes used as decoys

Genome assemblies are specified using NCBI RefSeq (`GCF_`) accessions.

## Requirements

- Python 3
- NCBI Datasets CLI
- PyYAML

Install the Python dependency:

```bash
python3 -m pip install -r requirements.txt