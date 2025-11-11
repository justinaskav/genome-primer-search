/*
 * Process: Generate thermodynamic analysis reports
 */
process GENERATE_THERMO_REPORTS {
    input:
    path json_files

    output:
    path "thermo_analysis.tsv", emit: detailed
    path "offtarget_summary.tsv", emit: offtarget
    path "primer_specificity.tsv", emit: specificity
    path "thermo_analysis.html", emit: html

    script:
    """
    generate_thermo_report.py \\
        ${json_files} \\
        --outdir . \\
        --min-size ${params.min_amplicon_size} \\
        --max-size ${params.max_amplicon_size}
    """
}
