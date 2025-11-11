/*
 * Process: Thermodynamic analysis of primer binding sites
 */
process THERMODYNAMIC_ANALYSIS {
    tag "$genome_input"
    conda 'conda-forge::biopython'

    input:
    tuple val(genome_input), path(primersearch_file), path(genome_fasta)
    path primers
    path reference
    val master_mix
    val annealing_temp
    val delta_tm_high
    val delta_tm_medium

    output:
    path("${input_safe}_thermo.json"), emit: results

    script:
    input_safe = genome_input.replaceAll(/\s+/, '_').replaceAll(/[^a-zA-Z0-9_.-]/, '')
    reference_arg = reference.name != 'NO_FILE' ? "--reference ${reference}" : ""
    """
    thermodynamic_analysis.py \\
        ${primersearch_file} \\
        ${primers} \\
        --genome-fasta ${genome_fasta} \\
        ${reference_arg} \\
        --genome-name "${genome_input}" \\
        --master-mix ${master_mix} \\
        --annealing-temp ${annealing_temp} \\
        --delta-tm-high ${delta_tm_high} \\
        --delta-tm-medium ${delta_tm_medium} \\
        -o ${input_safe}_thermo.json
    """
}
