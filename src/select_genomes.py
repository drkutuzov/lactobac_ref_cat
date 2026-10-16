import json

EXISTING_SPECIES = {
    "Lactobacillus acidophilus",
    "Lactobacillus johnsonii",
    "Lactobacillus gasseri",
    "Lactobacillus delbrueckii",
    "Lactobacillus crispatus",
    "Lactobacillus amylovorus",
}

with open("lactobacillus_candidates.json") as f:
    reports = json.load(f)["reports"]

available_species = {
    " ".join(r["organism"]["organism_name"].split()[:2])
    for r in reports
    if " sp. " not in r["organism"]["organism_name"]
}

new_species = available_species - EXISTING_SPECIES

print("Existing species:", len(EXISTING_SPECIES))
print("New species available:", len(new_species))
print("Total species possible:", len(available_species | EXISTING_SPECIES))


# Select one representative genome per new species
selected = []

for species in sorted(new_species):
    candidates = [
        r for r in reports
        if r["organism"]["organism_name"].split()[:2] == species.split()
    ]

    # Prefer NCBI reference genomes, then higher assembly N50
    candidates.sort(
        key=lambda r: (
            r.get("assembly_info", {}).get("refseq_category") == "reference genome",
            int(r.get("assembly_stats", {}).get("contig_n50", 0)),
            r["accession"],
        ),
        reverse=True,
    )

    best = candidates[0]
    selected.append(best)

    print(
        best["accession"],
        species,
        best["assembly_stats"]["total_sequence_length"],
    )

print("\nNew species representatives selected:", len(selected))


print("\nQuality assessment of selected genomes:")

for r in selected:
    accession = r["accession"]
    checkm = r.get("checkm_info", {})
    ani = r.get("average_nucleotide_identity", {})

    print(
        accession,
        "Completeness:", checkm.get("completeness", "NA"),
        "Contamination:", checkm.get("contamination", "NA"),
        "Taxonomy:", ani.get("taxonomy_check_status", "NA"),
    )


print("\nAdditional assemblies from existing species:")

for species in sorted(EXISTING_SPECIES):
    candidates = [
        r for r in reports
        if " ".join(r["organism"]["organism_name"].split()[:2]) == species
    ]

    print(f"{species}: {len(candidates)}")


from collections import defaultdict

EXISTING_ACCESSIONS = {
    "GCF_000011985.1",
    "GCF_000008065.1",
    "GCF_000014425.1",
    "GCF_000056065.1",
    "GCF_000091765.1",
    "GCF_000182855.2",
}

def species_name(report):
    return " ".join(report["organism"]["organism_name"].split()[:2])

def quality_rank(report):
    checkm = report.get("checkm_info", {})
    ani = report.get("average_nucleotide_identity", {})
    return (
        ani.get("taxonomy_check_status") == "OK",
        report.get("assembly_info", {}).get("refseq_category") == "reference genome",
        float(checkm.get("completeness") or 0),
        -float(checkm.get("contamination") or 100),
        report["accession"],
    )

# Group eligible assemblies by species
candidates_by_species = defaultdict(list)

used_accessions = EXISTING_ACCESSIONS | {
    r["accession"] for r in selected
}

for report in reports:
    name = report["organism"]["organism_name"]

    if " sp. " in name or report["accession"] in used_accessions:
        continue

    checkm = report.get("checkm_info", {})

    if (
        checkm.get("completeness") is None
        or checkm.get("contamination") is None
        or float(checkm["completeness"]) < 90
        or float(checkm["contamination"]) > 5
    ):
        continue

    candidates_by_species[species_name(report)].append(report)

for candidates in candidates_by_species.values():
    candidates.sort(key=quality_rank, reverse=True)

# Add one genome per species per round
additional = []

while len(additional) < 19:
    added_this_round = 0

    for species in sorted(candidates_by_species):
        if not candidates_by_species[species]:
            continue

        genome = candidates_by_species[species].pop(0)
        additional.append(genome)
        added_this_round += 1

        if len(additional) == 19:
            break

    if added_this_round == 0:
        raise RuntimeError("Not enough eligible genomes")

print("\nAdditional strains:")
for genome in additional:
    print(genome["accession"], species_name(genome))

print("\nTotal Lactobacillus genomes:", 6 + len(selected) + len(additional))


import csv
from pathlib import Path

manifest_path = Path("results/genome_selection_manifest.tsv")
manifest_path.parent.mkdir(parents=True, exist_ok=True)

fieldnames = [
    "accession",
    "species",
    "strain",
    "selection_group",
    "genome_length_bp",
    "completeness",
    "contamination",
    "taxonomy_check",
]

with manifest_path.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()

    for group, genomes in [
        ("new_species", selected),
        ("additional_strain", additional),
    ]:
        for genome in genomes:
            writer.writerow({
                "accession": genome["accession"],
                "species": species_name(genome),
                "strain": genome["organism"]
                    .get("infraspecific_names", {})
                    .get("strain", ""),
                "selection_group": group,
                "genome_length_bp": genome["assembly_stats"]
                    ["total_sequence_length"],
                "completeness": genome.get("checkm_info", {})
                    .get("completeness", ""),
                "contamination": genome.get("checkm_info", {})
                    .get("contamination", ""),
                "taxonomy_check": genome.get(
                    "average_nucleotide_identity", {}
                ).get("taxonomy_check_status", ""),
            })

print(f"\nManifest saved: {manifest_path}")


negative_controls = [
    ("GCF_009913655.1", "plantarum_candidates.json"),
    ("GCF_900636965.1", "rhamnosus_candidates.json"),
    ("GCF_001008015.1", "thermophilus_candidates.json"),
]

with manifest_path.open("a", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )

    for accession, filename in negative_controls:
        with open(filename) as source:
            candidates = json.load(source)["reports"]

        genome = next(
            r for r in candidates
            if r["accession"] == accession
        )

        writer.writerow({
            "accession": accession,
            "species": species_name(genome),
            "strain": genome["organism"]
                .get("infraspecific_names", {})
                .get("strain", ""),
            "selection_group": "negative_control",
            "genome_length_bp": genome["assembly_stats"]["total_sequence_length"],
            "completeness": genome.get("checkm_info", {}).get("completeness", ""),
            "contamination": genome.get("checkm_info", {}).get("contamination", ""),
            "taxonomy_check": genome.get(
                "average_nucleotide_identity", {}
            ).get("taxonomy_check_status", ""),
        })

print("Total genomes in manifest:", 6 + 2 + 39 + 3)