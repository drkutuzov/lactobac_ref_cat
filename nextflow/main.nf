nextflow.enable.dsl=2

params.catalogue = "results/lactobacillus_catalogue_v0.1.fasta"


process FASTA_STATS {

    publishDir "${projectDir}/output", mode: 'copy'

    input:
    path catalogue

    output:
    path "catalogue_stats.tsv"

    script:
    """
    python3 ${projectDir}/../scripts/fasta_stats.py \
        --input ${catalogue} \
        --output catalogue_stats.tsv
    """
}


workflow {

    catalogue_ch = Channel.fromPath(
        params.catalogue,
        checkIfExists: true
    )

    FASTA_STATS(catalogue_ch)
}