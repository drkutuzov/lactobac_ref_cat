from pathlib import Path
import csv
import yaml

CONFIG_PATH = "config.yaml"


def load_config(path):
    """Load pipeline configuration."""
    with open(path, "r") as file:
        return yaml.safe_load(file)


def parse_fasta(path):
    """Yield FASTA header and sequence pairs."""
    header = None
    sequence_parts = []

    with open(path, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)

                header = line[1:]
                sequence_parts = []
            else:
                sequence_parts.append(line)

        if header is not None:
            yield header, "".join(sequence_parts)


def calculate_n50(lengths):
    """Calculate N50 from a list of sequence lengths."""
    total_length = sum(lengths)
    halfway = total_length / 2

    cumulative = 0

    for length in sorted(lengths, reverse=True):
        cumulative += length

        if cumulative >= halfway:
            return length


def calculate_gc(sequences):
    """Calculate GC percentage across a collection of sequences."""
    total_length = sum(len(seq) for seq in sequences)

    if total_length == 0:
        return 0.0

    gc_count = sum(
        seq.upper().count("G") + seq.upper().count("C")
        for seq in sequences
    )

    return 100 * gc_count / total_length


def main():
    
    config = load_config(CONFIG_PATH)

    version = config["catalogue"]["version"]
    output_dir = Path(config["catalogue"]["output_dir"])

    fasta_path = (
        output_dir / f"lactobacillus_catalogue_v{version}.fasta"
    )

    genomes = {}

    for header, sequence in parse_fasta(fasta_path):
        accession = header.split("|")[0]

        if accession not in genomes:
            genomes[accession] = []

        genomes[accession].append(sequence)

    print(f"Genomes: {len(genomes)}")

    catalogue_length = sum(
        len(sequence)
        for sequences in genomes.values()
        for sequence in sequences
    )

    print(f"Total catalogue length: {catalogue_length:,} bp")

    results = []

    for accession, sequences in genomes.items():
        lengths = [len(seq) for seq in sequences]

        total_length = sum(lengths)
        n50 = calculate_n50(lengths)
        gc_percent = calculate_gc(sequences)

        results.append(
            {
                "accession": accession,
                "sequences": len(sequences),
                "total_length": total_length,
                "n50": n50,
                "gc_percent": round(gc_percent, 2),
            }
        )

        print(
            f"{accession}: "
            f"{len(sequences)} sequences, "
            f"{total_length:,} bp, "
            f"N50={n50:,} bp, "
            f"GC={gc_percent:.2f}%"
        )

        output_path = output_dir / f"catalogue_v{version}_stats.tsv"

    with open(output_path, "w", newline="") as file:
        fieldnames = [
            "accession",
            "sequences",
            "total_length",
            "n50",
            "gc_percent",
        ]

        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(results)

    print(f"\nStatistics written to {output_path}")



if __name__ == "__main__":
    main()