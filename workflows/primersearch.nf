/*
 * Genome Primer Search Workflow
 * Downloads genomes, runs primersearch, filters results, and generates reports
 */

// Import processes
include { DOWNLOAD_GENOME } from '../modules/local/download_genome'
include { EXTRACT_GENOME } from '../modules/local/extract_genome'
include { RUN_PRIMERSEARCH } from '../modules/local/run_primersearch'
include { FILTER_AMPLICONS } from '../modules/local/filter_amplicons'
include { GENERATE_REPORTS } from '../modules/local/generate_reports'


workflow PRIMERSEARCH {
    take:
    genomes_ch    // Channel with genome names/accessions
    primers_ch    // Channel with primers file

    main:
    // Step 1: Download genomes from NCBI
    DOWNLOAD_GENOME(genomes_ch)

    // Step 2: Extract FASTA from downloaded archives (or copy if already FASTA)
    EXTRACT_GENOME(DOWNLOAD_GENOME.out.genome_file)

    // Step 3: Run primersearch on each genome
    RUN_PRIMERSEARCH(
        EXTRACT_GENOME.out.genome_fasta,
        primers_ch.collect()
    )

    // Step 4: Filter amplicons (remove circular artifacts and invalid sizes)
    FILTER_AMPLICONS(RUN_PRIMERSEARCH.out.results)

    // Step 5: Generate comprehensive reports
    GENERATE_REPORTS(
        FILTER_AMPLICONS.out.filtered.map { it[1] }.collect(),
        FILTER_AMPLICONS.out.stats.map { it[1] }.collect()
    )

    emit:
    genomes = EXTRACT_GENOME.out.genome_fasta
    primersearch_raw = RUN_PRIMERSEARCH.out.results
    filtered_results = FILTER_AMPLICONS.out.filtered
    filter_stats = FILTER_AMPLICONS.out.stats
    summary_tsv = GENERATE_REPORTS.out.summary_tsv
    summary_html = GENERATE_REPORTS.out.summary_html
    amplicon_stats = GENERATE_REPORTS.out.amplicon_stats
    primer_stats = GENERATE_REPORTS.out.primer_stats
}
