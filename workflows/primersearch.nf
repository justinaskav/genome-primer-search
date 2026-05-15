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

    // Track download failures: DOWNLOAD_GENOME uses retry-then-ignore, so failed
    // genomes silently drop from the output channel. Left-join requested vs.
    // retrieved and emit the unmatched ones so users see what was skipped.
    failed_genomes_ch = genomes_ch
        .map { g -> tuple(g) }
        .join(DOWNLOAD_GENOME.out.genome_file, remainder: true)
        .filter { it[1] == null }
        .map { it[0] }
        .collectFile(
            name: 'failed_genomes.tsv',
            newLine: true,
            seed: 'genome\tstage'
        ) { "${it}\tdownload" }

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
    failed_genomes = failed_genomes_ch
    summary_tsv = GENERATE_REPORTS.out.summary_tsv
    summary_html = GENERATE_REPORTS.out.summary_html
    amplicon_stats = GENERATE_REPORTS.out.amplicon_stats
    primer_stats = GENERATE_REPORTS.out.primer_stats
}
