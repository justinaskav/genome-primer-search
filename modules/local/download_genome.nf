/*
 * Process: Download genome/sequence from NCBI
 * Supports:
 *  - Taxon names (e.g., "Escherichia coli")
 *  - Assembly accessions (e.g., "GCA_000005845.2", "GCF_000001405.40")
 *  - Nucleotide sequence accessions (e.g., "MN986463.1", "NC_000913.3", "U00096")
 */
process DOWNLOAD_GENOME {
    tag "$genome_input"

    conda 'conda-forge::ncbi-datasets-cli bioconda::entrez-direct'

    input:
    val genome_input

    output:
    tuple val(genome_input), path("${input_safe}.*"), emit: genome_file

    script:
    input_safe = genome_input.replaceAll(/\s+/, '_').replaceAll(/[^a-zA-Z0-9_.-]/, '')

    // Detect input type using pattern matching
    // INSDC nucleotide accessions: 1-2 letters + 5-8 digits (e.g., MN986463, AF123456, U00096)
    // RefSeq nucleotide accessions: 2 letters + underscore (e.g., NC_000913, NM_001234)
    // Assembly accessions: GCA_ or GCF_ prefix
    // Everything else: treated as taxon name

    def isNucleotideAccession = genome_input.matches(/^[A-Z]{1,2}[0-9]{5,8}(\.[0-9]+)?$/) ||
                               genome_input.matches(/^[A-Z]{2}_[0-9]+(\.[0-9]+)?$/)
    def isAssemblyAccession = genome_input.startsWith('GCA_') || genome_input.startsWith('GCF_')

    if (isNucleotideAccession) {
        // Use efetch for nucleotide sequences - outputs FASTA directly
        """
        efetch -db nuccore -id "${genome_input}" -format fasta > ${input_safe}.fasta

        # Verify download succeeded
        if [ ! -s ${input_safe}.fasta ]; then
            echo "Error: Failed to download nucleotide sequence ${genome_input}" >&2
            exit 1
        fi
        """
    } else if (isAssemblyAccession) {
        // Use datasets for assembly accessions - outputs zip
        """
        datasets download genome accession "${genome_input}" \
            --filename ${input_safe}.zip \
            --include genome
        """
    } else {
        // Treat as taxon name - outputs zip
        """
        datasets download genome taxon "${genome_input}" \
            --reference \
            --filename ${input_safe}.zip \
            --include genome
        """
    }
}
