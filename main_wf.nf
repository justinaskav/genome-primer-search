#!/usr/bin/env nextflow

/*
 * Pipeline: Genome Primer Search
 * Description: Downloads reference genomes and runs primersearch
 */

/*
 * Process 1: Download genome from NCBI for each organism
 * Supports both taxon names (e.g., "Escherichia coli") and assembly accessions (e.g., "GCA_000005845.2")
 */
process DOWNLOAD_GENOME {
    tag "$genome_input"

    conda 'conda-forge::ncbi-datasets-cli'

    input:
    val genome_input

    output:
    tuple val(genome_input), path("${input_safe}.zip"), emit: genome_zip

    script:
    input_safe = genome_input.replaceAll(/\s+/, '_').replaceAll(/[^a-zA-Z0-9_.-]/, '')

    // Detect if input is an accession ID (starts with GCA_ or GCF_) or a taxon name
    if (genome_input.startsWith('GCA_') || genome_input.startsWith('GCF_')) {
        // Input is an assembly accession
        """
        datasets download genome accession "${genome_input}" \
            --filename ${input_safe}.zip \
            --include genome
        """
    } else {
        // Input is a taxon name
        """
        datasets download genome taxon "${genome_input}" \
            --reference \
            --filename ${input_safe}.zip \
            --include genome
        """
    }
}

/*
 * Process 2: Extract genome FASTA from downloaded zip
 */
process EXTRACT_GENOME {
    tag "$genome_input"

    input:
    tuple val(genome_input), path(genome_zip)

    output:
    tuple val(genome_input), path("${input_safe}_genome.fasta"), emit: genome_fasta

    script:
    input_safe = genome_input.replaceAll(/\s+/, '_').replaceAll(/[^a-zA-Z0-9_.-]/, '')

    """
    unzip -q ${genome_zip}

    # Find the first genomic.fna file and rename it
    find ncbi_dataset/data -name "*.fna" -type f | head -n 1 | \
        xargs -I {} cp {} ${input_safe}_genome.fasta

    # Check if genome file is not empty
    if [ ! -s ${input_safe}_genome.fasta ]; then
        echo "Error: No genome file found for ${genome_input}" >&2
        exit 1
    fi
    """
}

/*
 * Process 3: Run primersearch on each genome
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

/*
 * Process 4: Collect and summarize all results
 */
process COLLECT_RESULTS {

    input:
    path(results)

    output:
    path("summary_report.txt")

    script:
    """
    echo "PRIMER SEARCH SUMMARY REPORT" > summary_report.txt
    echo "============================" >> summary_report.txt
    echo "" >> summary_report.txt
    echo "Generated: \$(date)" >> summary_report.txt
    echo "" >> summary_report.txt

    for result in ${results}; do
        echo "----------------------------------------" >> summary_report.txt
        echo "File: \$result" >> summary_report.txt
        echo "----------------------------------------" >> summary_report.txt

        # Count primers with hits (primers followed by Amplimer lines)
        primer_hit_count=\$(awk '/^Primer name / {primer=\$0; getline; if (/^Amplimer/) print primer}' \$result | wc -l | tr -d ' ')

        # Count total amplicons found
        amplicon_count=\$(grep -c "^Amplimer [0-9]" \$result || echo "0")

        if [ \$primer_hit_count -gt 0 ]; then
            echo "Primers with hits: \$primer_hit_count" >> summary_report.txt
            echo "Total amplicons found: \$amplicon_count" >> summary_report.txt
            echo "" >> summary_report.txt
            echo "Primer details:" >> summary_report.txt
            # Use awk to find "Primer name" lines that are followed by "Amplimer" and count amplicons per primer
            awk '/^Primer name / {
                primer=\$0;
                count=0;
                while (getline > 0) {
                    if (/^Amplimer/) count++;
                    else if (/^Primer name /) {
                        if (count > 0) print primer " (" count " amplicon" (count>1?"s":"") ")";
                        primer=\$0;
                        count=0;
                    }
                }
                if (count > 0) print primer " (" count " amplicon" (count>1?"s":"") ")";
            }' \$result | sed 's/Primer name /  - /' >> summary_report.txt
        else
            echo "No primers with hits" >> summary_report.txt
        fi

        echo "" >> summary_report.txt
    done

    echo "========================================" >> summary_report.txt
    echo "Pipeline completed successfully" >> summary_report.txt
    """
}

/*
 * Main workflow
 */
workflow {
    main:

    // Create channel from genome list (supports both taxon names and accession IDs)
    genomes_ch = Channel
        .fromPath(params.genomes)
        .splitText()
        .map { it.trim() }
        .filter { it != '' && !it.startsWith('#') }  // Filter empty lines and comments

    // Load primers file
    primers_ch = Channel.fromPath(params.primers)

    // Download genomes
    DOWNLOAD_GENOME(genomes_ch)

    // Extract genomes
    EXTRACT_GENOME(DOWNLOAD_GENOME.out.genome_zip)

    // Run primersearch
    RUN_PRIMERSEARCH(
        EXTRACT_GENOME.out.genome_fasta,
        primers_ch.collect()
    )

    // Collect all results
    COLLECT_RESULTS(RUN_PRIMERSEARCH.out.results.map { it[1] }.collect())

    publish:
    genome = EXTRACT_GENOME.out.genome_fasta
    primersearch_results = RUN_PRIMERSEARCH.out.results
    collect_results = COLLECT_RESULTS.out
}

output {
    genome {
        path 'genome'
    }
    primersearch_results {
        path 'primersearch_results'
    }
    collect_results {
        path '.'
    }
}
