/*
 * Process: Extract genome FASTA from downloaded file
 * Handles both zip files (from datasets) and FASTA files (from efetch)
 */
process EXTRACT_GENOME {
    tag "$genome_input"

    input:
    tuple val(genome_input), path(genome_file)

    output:
    tuple val(genome_input), path("${input_safe}_genome.fasta"), emit: genome_fasta

    script:
    input_safe = genome_input.replaceAll(/\s+/, '_').replaceAll(/[^a-zA-Z0-9_.-]/, '')

    // Check if input is already FASTA (from nucleotide download) or zip (from assembly download)
    if (genome_file.name.endsWith('.fasta')) {
        // Already FASTA from nucleotide sequence download - just copy
        """
        cp ${genome_file} ${input_safe}_genome.fasta

        # Verify file is not empty
        if [ ! -s ${input_safe}_genome.fasta ]; then
            echo "Error: Empty genome file for ${genome_input}" >&2
            exit 1
        fi
        """
    } else {
        // Zip file from assembly download - extract and rename
        """
        unzip -q ${genome_file}

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
}
