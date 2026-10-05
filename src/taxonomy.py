from pathlib import Path
import subprocess
import yaml
import csv

CONFIG_PATH = "config.yaml"


def load_config(path):
    """Load pipeline configuration."""
    with open(path, "r") as file:
        return yaml.safe_load(file)


def run_kraken(fasta_path, database_dir, report_path, output_path):
    """Run Kraken2 on one genome FASTA."""
    print(f"Running Kraken2: {fasta_path.name}")

    subprocess.run(
        [
            "kraken2",
            "--db", str(database_dir),
            "--report", str(report_path),
            "--output", str(output_path),
            str(fasta_path),
        ],
        check=True,
    )

def parse_dominant_genus(report_path):
    """Return the genus with the highest clade percentage."""
    genera = []

    with open(report_path, "r") as file:
        for line in file:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 6:
                continue

            percentage = float(fields[0])
            rank = fields[3]
            taxid = fields[4]
            name = fields[5].strip()

            if rank == "G":
                genera.append(
                    {
                        "percentage": percentage,
                        "taxid": taxid,
                        "name": name,
                    }
                )

    if not genera:
        return None

    return max(genera, key=lambda x: x["percentage"])


def main():
    config = load_config(CONFIG_PATH)

    raw_dir = Path(config["output_dir"])
    database_dir = Path("databases/kraken2")
    results_dir = Path("results/kraken2")
    results_dir.mkdir(parents=True, exist_ok=True)

    expected_genus = config["taxonomy"]["expected_genus"]

    results = []

    for genome in config["genomes"]:
        accession = genome["accession"]
        organism = genome["organism"]

        fasta_path = raw_dir / f"{accession}_genomic.fna"
        report_path = results_dir / f"{accession}.report"
        output_path = results_dir / f"{accession}.output"

        run_kraken(
            fasta_path,
            database_dir,
            report_path,
            output_path,
        )

        genus = parse_dominant_genus(report_path)

        if genus is None:
            status = "FAIL"
            genus_name = "Unclassified"
            percentage = 0.0
        else:
            genus_name = genus["name"]
            percentage = genus["percentage"]

            status = (
                "PASS"
                if genus_name == expected_genus
                else "FAIL"
            )

        print(
            f"{accession}: "
            f"{organism} → "
            f"{genus_name} ({percentage:.2f}%) → "
            f"{status}"
        )

        results.append(
            {
                "accession": accession,
                "organism": organism,
                "dominant_genus": genus_name,
                "genus_percent": percentage,
                "status": status,
            }
        )

        summary_path = Path("results/taxonomy_summary.tsv")


    with open(summary_path, "w", newline="") as file:
        fieldnames = [
            "accession",
            "organism",
            "dominant_genus",
            "genus_percent",
            "status",
        ]

        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(results)

    print(f"\nTaxonomy summary written to {summary_path}")


if __name__ == "__main__":
    main()