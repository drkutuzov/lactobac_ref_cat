configfile: "config.yaml"

VERSION = config["catalogue"]["version"]
OUTPUT_DIR = config["catalogue"]["output_dir"]

ACCESSIONS = [
    genome["accession"]
    for genome in config["genomes"]
]


rule all:
    input:
        f"{OUTPUT_DIR}/catalogue_v{VERSION}_report.tsv",
        f"{OUTPUT_DIR}/lactobacillus_catalogue_v{VERSION}.fasta",
        f"{OUTPUT_DIR}/catalogue_v{VERSION}_stats.tsv",


rule acquire:
    output:
        expand(
            "data/raw/{accession}_genomic.fna",
            accession=ACCESSIONS
        )
    shell:
        """
        python3 src/acquire.py
        """


rule qc:
    input:
        expand(
            "data/raw/{accession}_genomic.fna",
            accession=ACCESSIONS
        )
    output:
        "results/qc_summary.tsv"
    shell:
        """
        python3 src/qc.py
        """


rule taxonomy:
    input:
        expand(
            "data/raw/{accession}_genomic.fna",
            accession=ACCESSIONS
        )
    output:
        "results/taxonomy_summary.tsv"
    shell:
        """
        python3 src/taxonomy.py
        """


rule curate:
    input:
        qc="results/qc_summary.tsv",
        taxonomy="results/taxonomy_summary.tsv"
    output:
        report=f"{OUTPUT_DIR}/catalogue_v{VERSION}_report.tsv",
        catalogue=f"{OUTPUT_DIR}/lactobacillus_catalogue_v{VERSION}.fasta"
    shell:
        """
        python3 src/curate.py
        """


rule fasta_stats:
    input:
        catalogue=f"{OUTPUT_DIR}/lactobacillus_catalogue_v{VERSION}.fasta"
    output:
        stats=f"{OUTPUT_DIR}/catalogue_v{VERSION}_stats.tsv"
    shell:
        """
        python3 scripts/fasta_stats.py
        """