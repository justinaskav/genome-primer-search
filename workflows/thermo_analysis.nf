/*
 * Thermodynamic Analysis Workflow
 * Optional workflow for analyzing off-target binding thermodynamics
 */

// Import processes
include { THERMODYNAMIC_ANALYSIS } from '../modules/local/thermodynamic_analysis'
include { GENERATE_THERMO_REPORTS } from '../modules/local/generate_thermo_reports'


workflow THERMO_ANALYSIS {
    take:
    primersearch_results_ch  // Channel: tuple(val(genome_input), path(primersearch_file))
    primers_ch               // Channel: path(primers_file)
    reference_ch             // Channel: path(reference_sequences)

    main:
    // Run thermodynamic analysis on each primersearch result
    THERMODYNAMIC_ANALYSIS(
        primersearch_results_ch,
        primers_ch.collect(),
        reference_ch.collect(),
        params.thermo_master_mix,
        params.thermo_annealing_temp,
        params.thermo_delta_tm_high,
        params.thermo_delta_tm_medium
    )

    // Generate comprehensive reports from all results
    GENERATE_THERMO_REPORTS(
        THERMODYNAMIC_ANALYSIS.out.results.collect()
    )

    emit:
    detailed = GENERATE_THERMO_REPORTS.out.detailed
    offtarget = GENERATE_THERMO_REPORTS.out.offtarget
    specificity = GENERATE_THERMO_REPORTS.out.specificity
    html = GENERATE_THERMO_REPORTS.out.html
}
