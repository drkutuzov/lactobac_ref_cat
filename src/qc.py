from pathlib import Path
import csv
import subprocess
import yaml


CONFIG_PATH = "config.yaml"


def load_config(path):
    """Load pipeline configuration."""
    with open(path, "r") as file:
        return yaml.safe_load(file)


def run_quast(fasta_path, output_dir):
    """Run QUAST on one genome FASTA."""

    print(f"Running QUAST: {fasta_path.name}")

    subprocess.run(
        [
            "quast",
            str(fasta_path),
            "-o",
            str(output_dir),
            "--silent",
        ],
        check=True,
    )


def parse_quast_report(report_path):
    """Extract selected assembly statistics from QUAST report.tsv."""

    metrics = {}

    with open(report_path, "r") as file:
        reader = csv.reader(file, delimiter="\t")

        for row in reader:
            if len(row) < 2:
                continue

            metric = row[0]
            value = row[1]

            if metric == "# contigs":
                metrics["contigs"] = int(value)

            elif metric == "Total length":
                metrics["total_length"] = int(value)

            elif metric == "GC (%)":
                metrics["gc_percent"] = float(value)

            elif metric == "N50":
                metrics["n50"] = int(value)

    return metrics


def main():
    config = load_config(CONFIG_PATH)

    raw_dir = Path(config["output_dir"])
    quast_dir = Path("results/quast")
    quast_dir.mkdir(parents=True, exist_ok=True)

    min_n50 = config["qc"]["min_n50"]

    results = []

    for genome in config["genomes"]:
        accession = genome["accession"]

        fasta_path = raw_dir / f"{accession}_genomic.fna"
        genome_output = quast_dir / accession

        run_quast(fasta_path, genome_output)

        report_path = genome_output / "report.tsv"
        metrics = parse_quast_report(report_path)

        status = "PASS" if metrics["n50"] >= min_n50 else "FAIL"

        results.append(
            {
                "accession": accession,
                "organism": genome["organism"],
                "contigs": metrics["contigs"],
                "total_length": metrics["total_length"],
                "gc_percent": metrics["gc_percent"],
                "n50": metrics["n50"],
                "status": status,
            }
)

        print(
            f"{accession}: "
            f"N50={metrics['n50']:,} bp → {status}"
        )

    summary_path = Path("results/qc_summary.tsv")

    with open(summary_path, "w", newline="") as file:
        fieldnames = [
            "accession",
            "organism",
            "contigs",
            "total_length",
            "gc_percent",
            "n50",
            "status",
        ]

        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(results)

    print(f"\nQC summary written to {summary_path}")


if __name__ == "__main__":
    main()