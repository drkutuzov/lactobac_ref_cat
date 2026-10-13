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
        "data/raw/{accession}_genomic.fna"
    conda:
        "workflow/envs/acquire.yaml"
    shell:
        """
        python3 src/acquire.py {wildcards.accession}
        """


rule qc:
    input:
        expand(
            "data/raw/{accession}_genomic.fna",
            accession=ACCESSIONS
        )
    output:
        "results/qc_summary.tsv"
    conda: "workflow/envs/quast.yaml"
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
    conda:
        "workflow/envs/kraken2.yaml"
    shell:
        """
        python3 src/taxonomy.py
        """

rule checkm2:
    input:
        expand(
            "data/raw/{accession}_genomic.fna",
            accession=ACCESSIONS
        )
    output:
        "results/checkm2/quality_report.tsv"
    threads: 4
    conda:
        "workflow/envs/checkm2.yaml"
    shell:
        """
        checkm2 predict \
            --input data/raw \
            --output-directory results/checkm2 \
            --extension .fna \
            --threads {threads} \
            --force
        """


rule curate:
    input:
        qc="results/qc_summary.tsv",
        taxonomy="results/taxonomy_summary.tsv",
        checkm2="results/checkm2/quality_report.tsv"
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
