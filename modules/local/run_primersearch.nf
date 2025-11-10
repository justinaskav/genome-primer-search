/*
 * Process: Run primersearch on each genome
 */
process RUN_PRIMERSEARCH {
    tag "$genome_input"
    conda 'bioconda::emboss'

    input:
    tuple val(genome_input), path(genome_fasta)
    path primers

    output:
    tuple val(genome_input), path("${input_safe}_primersearch.out"), emit: results

    script:
    input_safe = genome_input.replaceAll(/\s+/, '_').replaceAll(/[^a-zA-Z0-9_.-]/, '')
    """
    primersearch \
        -seqall ${genome_fasta} \
        -infile ${primers} \
        -mismatchpercent ${params.mismatch} \
        -outfile ${input_safe}_primersearch.out
    """
}
