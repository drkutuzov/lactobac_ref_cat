
from collections import Counter
from pathlib import Path
import csv

import argparse

parser = argparse.ArgumentParser(
    description="Validate a curated bacterial genome catalogue"
)
parser.add_argument("--report", type=Path, required=True)
parser.add_argument("--fasta", type=Path, required=True)
args = parser.parse_args()

REPORT = args.report
FASTA = args.fasta

# Read accepted and rejected accessions.
with REPORT.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

accepted = {row["accession"] for row in rows if row["final_status"] == "PASS"}
rejected = {row["accession"] for row in rows if row["final_status"] == "FAIL"}

# Read FASTA headers and count sequence lengths.
genome_counts = Counter()
sequence_ids = set()
total_length = 0
record_count = 0

with FASTA.open() as f:
    for line in f:
        if line.startswith(">"):
            record_count += 1
            identifier = line[1:].split()[0]
            genome_accession, sequence_accession = identifier.split("|", 1)

            assert identifier not in sequence_ids, f"Duplicate sequence: {identifier}"
            sequence_ids.add(identifier)
            genome_counts[genome_accession] += 1
        else:
            total_length += len(line.strip())

observed = set(genome_counts)

assert observed == accepted, (
    f"Missing genomes: {accepted - observed}; "
    f"Unexpected genomes: {observed - accepted}"
)
assert not observed & rejected, "Rejected genomes found in FASTA"
assert total_length > 0, "Catalogue FASTA is empty"

print("Catalogue validation PASSED")
print(f"Accepted genomes: {len(accepted)}")
print(f"FASTA records: {record_count}")
print(f"Total length: {total_length:,} bp")



import hashlib


def fasta_checksums(path):
    """Return sequence IDs and SHA-256 hashes of their nucleotide sequences."""
    checksums = {}
    sequence_id = None
    digest = None

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if sequence_id is not None:
                    assert sequence_id not in checksums, f"Duplicate ID: {sequence_id}"
                    checksums[sequence_id] = digest.hexdigest()

                sequence_id = line[1:].split()[0]
                digest = hashlib.sha256()
            else:
                assert digest is not None, "Sequence before FASTA header"
                digest.update(line.strip().encode())

    if sequence_id is not None:
        assert sequence_id not in checksums, f"Duplicate ID: {sequence_id}"
        checksums[sequence_id] = digest.hexdigest()

    return checksums


# Calculate checksums for all curated sequences.
curated = fasta_checksums(FASTA)

for accession in sorted(accepted):
    original_path = Path(f"data/raw/{accession}_genomic.fna")
    original = fasta_checksums(original_path)

    # Curated headers contain the assembly accession as a prefix.
    expected_checksums = {
        f"{accession}|{seq_id}": checksum
        for seq_id, checksum in original.items()
    }

    observed_checksums = {
        seq_id: checksum
        for seq_id, checksum in curated.items()
        if seq_id.startswith(f"{accession}|")
    }

    assert observed_checksums == expected_checksums, (
        f"Sequence mismatch for {accession}"
    )

print(f"Sequence integrity PASSED: {len(accepted)} genomes")
