/*
 * Process: Generate comprehensive summary reports
 * Creates TSV, HTML, and statistics files from all filtered results
 */
process GENERATE_REPORTS {

    input:
    path(filtered_results)
    path(stats_files)

    output:
    path("summary.tsv"), emit: summary_tsv
    path("summary.html"), emit: summary_html
    path("amplicon_stats.tsv"), emit: amplicon_stats
    path("primer_stats.tsv"), emit: primer_stats

    script:
    """
    generate_html_report.py \
        --filtered-dir . \
        --stats-dir . \
        --output-dir . \
        --min-size ${params.min_amplicon_size} \
        --max-size ${params.max_amplicon_size}
    """
}
