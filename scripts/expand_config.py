import csv
from pathlib import Path

import yaml


config_path = Path("config.yaml")
manifest_path = Path("metadata/genome_selection_manifest.tsv")

with config_path.open() as file:
    config = yaml.safe_load(file)

existing_accessions = {
    genome["accession"] for genome in config["genomes"]
}

with manifest_path.open() as file:
    reader = csv.DictReader(file, delimiter="\t")

    for row in reader:
        accession = row["accession"]

        if accession in existing_accessions:
            raise ValueError(f"Duplicate accession: {accession}")

        organism = f'{row["species"]} {row["strain"]}'.strip()

        config["genomes"].append({
            "accession": accession,
            "organism": organism,
        })

        existing_accessions.add(accession)

print(f"Total genomes: {len(config['genomes'])}")

if len(config["genomes"]) != 50:
    raise ValueError("Expected exactly 50 genomes")

with config_path.open("w") as file:
    yaml.safe_dump(config, file, sort_keys=False)

print(f"Updated {config_path}")