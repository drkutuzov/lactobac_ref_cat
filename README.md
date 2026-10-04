# Lactobacillus Genome Pipeline

A small bioinformatics pipeline for acquiring and analysing bacterial reference genomes.

The pipeline reads NCBI RefSeq assembly accessions from `config.yaml` and downloads the corresponding genomic FASTA files using NCBI Datasets.

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