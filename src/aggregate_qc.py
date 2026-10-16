from pathlib import Path
import csv
import sys


def main():
    input_files = [Path(path) for path in sys.argv[1:]]

    if not input_files:
        raise ValueError("No QC files provided")

    output_path = Path("results/qc_summary.tsv")

    with output_path.open("w", newline="") as output:
        writer = None

        for path in input_files:
            with path.open() as file:
                reader = csv.DictReader(file, delimiter="\t")

                if writer is None:
                    writer = csv.DictWriter(
                        output,
                        fieldnames=reader.fieldnames,
                        delimiter="\t",
                    )
                    writer.writeheader()

                writer.writerows(reader)

    print(f"Aggregated {len(input_files)} genomes into {output_path}")


if __name__ == "__main__":
    main()