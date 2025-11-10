#!/usr/bin/env nextflow

/*
 * Genome Primer Search Pipeline
 * Description: Downloads reference genomes and searches for primer binding sites
 * Author: Justinas Kavoliūnas
 * Version: 0.2.0
 */

// Import main workflow
include { PRIMERSEARCH } from './workflows/primersearch'

// Entry workflow
workflow {
    main:
    // Validate required parameters
    if (!params.genomes || !params.primers) {
        error "Required parameters missing. Please specify --genomes and --primers"
    }

    // Create channel from genome list
    // Supports both taxon names and assembly accessions (GCA_/GCF_)
    genomes_ch = Channel
        .fromPath(params.genomes)
        .splitText()
        .map { it.trim() }
        .filter { it != '' && !it.startsWith('#') }  // Filter empty lines and comments

    // Load primers file
    primers_ch = Channel.fromPath(params.primers)

    // Run main workflow
    PRIMERSEARCH(genomes_ch, primers_ch)

    publish:
    genomes = PRIMERSEARCH.out.genomes
    primersearch_raw = PRIMERSEARCH.out.primersearch_raw
    filtered_results = PRIMERSEARCH.out.filtered_results
    filter_stats = PRIMERSEARCH.out.filter_stats
    summary_tsv = PRIMERSEARCH.out.summary_tsv
    summary_html = PRIMERSEARCH.out.summary_html
    amplicon_stats = PRIMERSEARCH.out.amplicon_stats
    primer_stats = PRIMERSEARCH.out.primer_stats
}

// Output configuration
output {
    genomes {
        path 'results/genomes'
    }
    primersearch_raw {
        path 'results/primersearch'
    }
    filtered_results {
        path 'results/filtered'
    }
    filter_stats {
        path 'results/filter_stats'
    }
    summary_tsv {
        path 'reports'
    }
    summary_html {
        path 'reports'
    }
    amplicon_stats {
        path 'reports'
    }
    primer_stats {
        path 'reports'
    }
}
