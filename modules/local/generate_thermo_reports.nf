/*
 * Process: Generate thermodynamic analysis reports
 */
process GENERATE_THERMO_REPORTS {
    conda 'conda-forge::python>=3.8'

    input:
    path json_files

    output:
    path "thermo_reports/thermo_analysis.tsv", emit: detailed
    path "thermo_reports/offtarget_summary.tsv", emit: offtarget
    path "thermo_reports/primer_specificity.tsv", emit: specificity
    path "thermo_reports/thermo_analysis.html", emit: html

    script:
    """
    mkdir -p thermo_reports

    generate_thermo_report.py \\
        ${json_files} \\
        --outdir thermo_reports
    """
}
