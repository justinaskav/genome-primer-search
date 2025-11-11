#!/usr/bin/env nextflow

/*
 * Genome Primer Search Pipeline
 * Description: Downloads reference genomes and searches for primer binding sites
 * Author: Justinas Kavoliūnas
 * Version: 0.2.0
 */

// Import workflows
include { PRIMERSEARCH } from './workflows/primersearch'
include { THERMO_ANALYSIS } from './workflows/thermo_analysis'

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

    // Optional: Run thermodynamic analysis
    if (params.enable_thermo_analysis) {
        // Load reference sequences (if provided)
        if (params.thermo_references && file(params.thermo_references).exists()) {
            reference_ch = Channel.fromPath(params.thermo_references)
        } else {
            // Create empty file channel if no reference provided
            reference_ch = Channel.value(file('NO_FILE'))
        }

        THERMO_ANALYSIS(
            PRIMERSEARCH.out.filtered_results,
            PRIMERSEARCH.out.genomes,
            primers_ch,
            reference_ch
        )

        thermo_detailed_ch = THERMO_ANALYSIS.out.detailed
        thermo_offtarget_ch = THERMO_ANALYSIS.out.offtarget
        thermo_specificity_ch = THERMO_ANALYSIS.out.specificity
        thermo_html_ch = THERMO_ANALYSIS.out.html
    } else {
        thermo_detailed_ch = Channel.empty()
        thermo_offtarget_ch = Channel.empty()
        thermo_specificity_ch = Channel.empty()
        thermo_html_ch = Channel.empty()
    }

    publish:
    genomes = PRIMERSEARCH.out.genomes
    primersearch_raw = PRIMERSEARCH.out.primersearch_raw
    filtered_results = PRIMERSEARCH.out.filtered_results
    filter_stats = PRIMERSEARCH.out.filter_stats
    summary_tsv = PRIMERSEARCH.out.summary_tsv
    summary_html = PRIMERSEARCH.out.summary_html
    amplicon_stats = PRIMERSEARCH.out.amplicon_stats
    primer_stats = PRIMERSEARCH.out.primer_stats
    thermo_detailed = thermo_detailed_ch
    thermo_offtarget = thermo_offtarget_ch
    thermo_specificity = thermo_specificity_ch
    thermo_html = thermo_html_ch
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
    thermo_detailed {
        path 'thermo_reports'
    }
    thermo_offtarget {
        path 'thermo_reports'
    }
    thermo_specificity {
        path 'thermo_reports'
    }
    thermo_html {
        path 'thermo_reports'
    }
}
