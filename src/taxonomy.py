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


def parse_sequence_assignments(output_path):
    """Read Kraken2 classifications and lengths for individual sequences."""
    assignments = []

    with open(output_path, "r") as file:
        for line in file:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 4:
                continue

            assignments.append(
                {
                    "classified": fields[0] == "C",
                    "sequence_id": fields[1],
                    "taxid": fields[2],
                    "length": int(fields[3]),
                }
            )

    return assignments


def load_taxonomy_nodes(nodes_path):
    """Load taxonomic parent relationships and ranks."""
    taxonomy = {}

    with open(nodes_path, "r") as file:
        for line in file:
            fields = [field.strip() for field in line.split("|")]

            taxid = fields[0]
            parent_taxid = fields[1]
            rank = fields[2]

            taxonomy[taxid] = {
                "parent": parent_taxid,
                "rank": rank,
            }

    return taxonomy


def get_genus_taxid(taxid, taxonomy):
    """Find the genus ancestor of a taxonomic ID."""

    while taxid in taxonomy:
        node = taxonomy[taxid]

        if node["rank"] == "genus":
            return taxid

        parent = node["parent"]

        if parent == taxid:
            break

        taxid = parent

    return None


def calculate_genus_bp_percent(assignments, taxonomy, expected_genus_taxid):
    """Calculate assembly bp percentages assigned to the expected genus."""

    total_bp = 0
    expected_genus_bp = 0
    unclassified_bp = 0

    for assignment in assignments:
        length = assignment["length"]
        total_bp += length

        if not assignment["classified"]:
            unclassified_bp += length
            continue

        genus_taxid = get_genus_taxid(
            assignment["taxid"],
            taxonomy,
        )

        if genus_taxid == expected_genus_taxid:
            expected_genus_bp += length

    if total_bp == 0:
        return 0.0, 0.0

    expected_percent = 100 * expected_genus_bp / total_bp
    unclassified_percent = 100 * unclassified_bp / total_bp

    return expected_percent, unclassified_percent


def find_genus_taxid(genus_name, names_path, taxonomy):
    """Find the NCBI taxid for a genus using its scientific name."""

    with open(names_path, "r") as file:
        for line in file:
            fields = [field.strip() for field in line.split("|")]

            taxid = fields[0]
            name = fields[1]
            name_class = fields[3]

            if (
                name == genus_name
                and name_class == "scientific name"
                and taxid in taxonomy
                and taxonomy[taxid]["rank"] == "genus"
            ):
                return taxid

    raise ValueError(f"Genus not found in taxonomy: {genus_name}")


def main():
    config = load_config(CONFIG_PATH)

    raw_dir = Path(config["output_dir"])
    database_dir = Path(config["databases"]["kraken2"])
    results_dir = Path("results/kraken2")
    results_dir.mkdir(parents=True, exist_ok=True)

    expected_genus = config["taxonomy"]["expected_genus"]
    taxonomy = load_taxonomy_nodes(
        database_dir / "nodes.dmp"
    )

    expected_genus_taxid = find_genus_taxid(
        expected_genus,
        database_dir / "names.dmp",
        taxonomy,
    )

    min_genus_bp_percent = config["taxonomy"]["min_expected_genus_bp_percent"]

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
        assignments = parse_sequence_assignments(output_path)

        genus_bp_percent, unclassified_bp_percent = calculate_genus_bp_percent(
            assignments,
            taxonomy,
            expected_genus_taxid,
        )


        if genus is None:
            genus_name = "Unclassified"
            percentage = 0.0
        else:
            genus_name = genus["name"]
            percentage = genus["percentage"]

        status = (
            "PASS"
            if genus_name == expected_genus
            and genus_bp_percent >= min_genus_bp_percent
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
                "expected_genus_bp_percent": genus_bp_percent,
                "unclassified_bp_percent": unclassified_bp_percent,
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
            "expected_genus_bp_percent",
            "unclassified_bp_percent",
            "status",
        ]

        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(results)

    print(f"\nTaxonomy summary written to {summary_path}")


if __name__ == "__main__":
    main()