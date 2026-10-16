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
        f"{OUTPUT_DIR}/catalogue_v{VERSION}_validated.ok",



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
        "data/raw/{accession}_genomic.fna"
    output:
        "results/qc/{accession}.tsv"
    conda:
        "workflow/envs/quast.yaml"
    shell:
        """
        python3 src/qc.py {wildcards.accession}
        """


rule aggregate_qc:
    input:
        expand(
            "results/qc/{accession}.tsv",
            accession=ACCESSIONS
        )
    output:
        "results/qc_summary.tsv"
    shell:
        """
        python3 src/aggregate_qc.py {input}
        """


rule taxonomy:
    input:
        "data/raw/{accession}_genomic.fna"
    output:
        "results/taxonomy/{accession}.tsv"
    conda:
        "workflow/envs/kraken2.yaml"
    shell:
        "python3 src/taxonomy.py {wildcards.accession}"


rule aggregate_taxonomy:
    input:
        expand(
            "results/taxonomy/{accession}.tsv",
            accession=ACCESSIONS
        )
    output:
        "results/taxonomy_summary.tsv"
    shell:
        "python3 src/aggregate_taxonomy.py {input}"


rule checkm2:
    input:
        "data/raw/{accession}_genomic.fna"
    output:
        "results/checkm2/{accession}/quality_report.tsv"
    threads: 4
    conda:
        "workflow/envs/checkm2.yaml"
    shell:
        """
        checkm2 predict \
            --input {input} \
            --output-directory results/checkm2/{wildcards.accession} \
            --threads {threads} \
            --force
        """


rule aggregate_checkm2:
    input:
        expand(
            "results/checkm2/{accession}/quality_report.tsv",
            accession=ACCESSIONS
        )
    output:
        "results/checkm2/quality_report.tsv"
    shell:
        "python3 src/aggregate_checkm2.py {input}"


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


rule validate_catalogue:
    input:
        report=f"{OUTPUT_DIR}/catalogue_v{VERSION}_report.tsv",
        catalogue=f"{OUTPUT_DIR}/lactobacillus_catalogue_v{VERSION}.fasta",
        genomes=expand(
            "data/raw/{accession}_genomic.fna",
            accession=ACCESSIONS
        )
    output:
        touch(f"{OUTPUT_DIR}/catalogue_v{VERSION}_validated.ok")
    shell:
        """
        python3 tests/validate_catalogue.py \
            --report {input.report} \
            --fasta {input.catalogue}
        """