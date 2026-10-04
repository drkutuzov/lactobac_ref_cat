from pathlib import Path
import shutil
import subprocess
import zipfile

import yaml


CONFIG_PATH = "config.yaml"


def load_config(path):
    """Load genome accessions and settings from YAML."""
    with open(path, "r") as file:
        return yaml.safe_load(file)
    

def validate_fasta(fasta_path):
    """Check that a FASTA file exists, is non-empty, and has a header."""

    if not fasta_path.exists():
        raise RuntimeError(f"FASTA file does not exist: {fasta_path}")

    if fasta_path.stat().st_size == 0:
        raise RuntimeError(f"FASTA file is empty: {fasta_path}")

    with open(fasta_path, "r") as file:
        first_line = file.readline()

    if not first_line.startswith(">"):
        raise RuntimeError(f"Invalid FASTA file: {fasta_path}")
    

def download_genome(accession, output_dir):
    """Download and extract a genome FASTA from NCBI."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    zip_path = output_dir / f"{accession}.zip"
    extract_dir = output_dir / accession
    final_fasta = output_dir / f"{accession}_genomic.fna"

    if final_fasta.exists():
        print(f"Skipping {accession}: FASTA already exists.")
        return

    print(f"Downloading {accession}...")

    subprocess.run(
        [
            "datasets",
            "download",
            "genome",
            "accession",
            accession,
            "--include",
            "genome",
            "--filename",
            str(zip_path),
        ],
        check=True,
    )

    print(f"Extracting {accession}...")

    with zipfile.ZipFile(zip_path, "r") as archive:
        archive.extractall(extract_dir)

    fasta_files = list(extract_dir.rglob("*_genomic.fna"))

    if len(fasta_files) != 1:
        raise RuntimeError(
            f"Expected 1 genomic FASTA for {accession}, "
            f"but found {len(fasta_files)}."
        )

    shutil.copy(fasta_files[0], final_fasta)
    validate_fasta(final_fasta)

    print(f"Saved FASTA: {final_fasta}")

    # Remove temporary files
    zip_path.unlink()
    shutil.rmtree(extract_dir)


def main():
    config = load_config(CONFIG_PATH)

    for genome in config["genomes"]:
        accession = genome["accession"]
        download_genome(accession, config["output_dir"])


if __name__ == "__main__":
    main()