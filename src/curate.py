from pathlib import Path
import csv
import yaml

CONFIG_PATH = "config.yaml"


def load_config(path):
    """Load pipeline configuration."""
    with open(path, "r") as file:
        return yaml.safe_load(file)


def load_tsv(path):
    """Load a TSV file indexed by accession."""
    records = {}

    with open(path, "r") as file:
        reader = csv.DictReader(file, delimiter="\t")

        for row in reader:
            records[row["accession"]] = row

    return records


def main():
    config = load_config(CONFIG_PATH)

    qc_results = load_tsv("results/qc_summary.tsv")
    taxonomy_results = load_tsv("results/taxonomy_summary.tsv")
    checkm2_results = {}

    with open("results/checkm2/quality_report.tsv") as file:
        reader = csv.DictReader(file, delimiter="\t")

        for row in reader:
            accession = row["Name"].removesuffix("_genomic")
            checkm2_results[accession] = row

    curated_records = []

    for genome in config["genomes"]:
        accession = genome["accession"]
        organism = genome["organism"]

        qc = qc_results[accession]
        taxonomy = taxonomy_results[accession]
        checkm2 = checkm2_results[accession]

        completeness = float(checkm2["Completeness"])
        contamination = float(checkm2["Contamination"])

        reasons = []

        if qc["status"] != "PASS":
            reasons.append(
                f"N50 below threshold ({qc['n50']} bp)"
            )

        if taxonomy["status"] != "PASS":
            expected_genus = config["taxonomy"]["expected_genus"]
            min_bp = config["taxonomy"]["min_expected_genus_bp_percent"]
            actual_bp = float(taxonomy["expected_genus_bp_percent"])

            if taxonomy["dominant_genus"] != expected_genus:
                reasons.append(
                    f"Unexpected genus: {taxonomy['dominant_genus']}"
                )

            if actual_bp < min_bp:
                reasons.append(
                    f"Expected genus bp below threshold "
                    f"({actual_bp:.2f}% < {min_bp:.2f}%)"
                )

        if completeness < config["qc"]["min_completeness"]:
            reasons.append(
                f"Completeness below threshold ({completeness:.2f}%)"
            )

        if contamination > config["qc"]["max_contamination"]:
            reasons.append(
                f"Contamination above threshold ({contamination:.2f}%)"
            )

        if reasons:
            final_status = "FAIL"
            reason = "; ".join(reasons)
        else:
            final_status = "PASS"
            reason = "Passed all gates"

        curated_records.append(
            {
                "accession": accession,
                "organism": organism,
                "n50": qc["n50"],
                "completeness": completeness,
                "contamination": contamination,
                "qc_status": qc["status"],
                "dominant_genus": taxonomy["dominant_genus"],
                "genus_percent": taxonomy["genus_percent"],
                "expected_genus_bp_percent": taxonomy["expected_genus_bp_percent"],
                "unclassified_bp_percent": taxonomy["unclassified_bp_percent"],
                "taxonomy_status": taxonomy["status"],
                "final_status": final_status,
                "reason": reason,
            }
        )

        print(
            f"{accession}: {final_status} — {reason}"
        )

    version = config["catalogue"]["version"]
    output_dir = Path(config["catalogue"]["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    report_path = output_dir / f"catalogue_v{version}_report.tsv"

    fieldnames = [
        "accession",
        "organism",
        "n50",
        "completeness", 
        "contamination",
        "qc_status",
        "dominant_genus",
        "genus_percent",
        "expected_genus_bp_percent",
        "unclassified_bp_percent",
        "taxonomy_status",
        "final_status",
        "reason",
    ]

    with open(report_path, "w", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(curated_records)

    print(f"\nCatalogue report written to {report_path}")

    catalogue_path = (
        output_dir / f"lactobacillus_catalogue_v{version}.fasta"
    )

    raw_dir = Path(config["output_dir"])

    included_count = 0

    with open(catalogue_path, "w") as catalogue_file:
        for record in curated_records:
            if record["final_status"] != "PASS":
                continue

            accession = record["accession"]
            fasta_path = raw_dir / f"{accession}_genomic.fna"

            with open(fasta_path, "r") as fasta_file:
                for line in fasta_file:
                    if line.startswith(">"):
                        original_header = line[1:].rstrip("\n")
                        catalogue_file.write(
                            f">{accession}|{original_header}\n"
                        )
                    else:
                        catalogue_file.write(line)

            included_count += 1

    print(
        f"Catalogue FASTA written to {catalogue_path} "
        f"({included_count} genomes included)"
    )

if __name__ == "__main__":
    main()