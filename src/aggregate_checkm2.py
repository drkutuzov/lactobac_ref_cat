from pathlib import Path
import csv
import sys


def main():
    input_files = [Path(path) for path in sys.argv[1:]]

    if not input_files:
        raise ValueError("No CheckM2 reports provided")

    output_path = Path("results/checkm2/quality_report.tsv")

    rows = []
    fieldnames = None

    for path in input_files:
        with path.open(newline="") as file:
            reader = csv.DictReader(file, delimiter="\t")

            if fieldnames is None:
                fieldnames = reader.fieldnames
            elif reader.fieldnames != fieldnames:
                raise ValueError(f"Inconsistent columns in {path}")

            rows.extend(reader)

    rows.sort(key=lambda row: row["Name"])

    with output_path.open("w", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Aggregated {len(rows)} genomes into {output_path}")


if __name__ == "__main__":
    main()