/*
 * Process: Filter amplicons based on size and other criteria
 * Removes likely circular genome artifacts and invalid amplicons
 */
process FILTER_AMPLICONS {
    tag "$genome_input"

    input:
    tuple val(genome_input), path(primersearch_out)

    output:
    tuple val(genome_input), path("${input_safe}_filtered.out"), emit: filtered
    tuple val(genome_input), path("${input_safe}_stats.json"), emit: stats

    script:
    input_safe = genome_input.replaceAll(/\s+/, '_').replaceAll(/[^a-zA-Z0-9_.-]/, '')
    """
    filter_amplicons.py \
        --input ${primersearch_out} \
        --output ${input_safe}_filtered.out \
        --stats ${input_safe}_stats.json \
        --max-size ${params.max_amplicon_size} \
        --min-size ${params.min_amplicon_size}
    """
}
