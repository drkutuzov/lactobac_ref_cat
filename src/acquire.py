from pathlib import Path
import shutil
import sys
import tempfile
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
    """Download one genome and validate it before saving."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    final_fasta = output_dir / f"{accession}_genomic.fna"

    if final_fasta.exists():
        validate_fasta(final_fasta)
        print(f"Skipping {accession}: valid FASTA already exists.")
        return

    # Temporary files are kept separate from the final output.
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        zip_path = tmp_dir / f"{accession}.zip"
        extract_dir = tmp_dir / "extracted"

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

        with zipfile.ZipFile(zip_path) as archive:
            archive.extractall(extract_dir)

        fasta_files = list(extract_dir.rglob("*_genomic.fna"))

        if len(fasta_files) != 1:
            raise RuntimeError(
                f"Expected 1 genomic FASTA for {accession}, "
                f"but found {len(fasta_files)}."
            )

        # Validate before moving into the final data directory.
        validate_fasta(fasta_files[0])
        shutil.move(str(fasta_files[0]), str(final_fasta))

    print(f"Saved FASTA: {final_fasta}")


def main():
    config = load_config(CONFIG_PATH)

    if len(sys.argv) != 2:
        raise SystemExit("Usage: python3 src/acquire.py ACCESSION")

    accession = sys.argv[1]

    allowed_accessions = {
        genome["accession"] for genome in config["genomes"]
    }

    if accession not in allowed_accessions:
        raise ValueError(f"Accession not found in config.yaml: {accession}")

    download_genome(accession, config["output_dir"])


if __name__ == "__main__":
    main()