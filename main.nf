#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.output_folder = false
params.output_prefix = false
params.genomes = false
params.genes = false
params.min_coverage = 50
params.min_identity = 50
params.ftp_threads = 25

// Docker containers reused across processes
container__wget = "quay.io/fhcrc-microbiome/wget:latest"
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
container__mmseqs = "quay.io/fhcrc-microbiome/mmseqs2:version-12"
container__mashtree = "quay.io/hdc-workflows/mashtree:1.2.0"

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/gig-map <ARGUMENTS>

    Required Arguments:
      --genomes             Genome sequences in FASTA format (see note below)
      --genome_tables       Tables of NCBI genomes to analyze (see note below)
      --genes               Amino acid sequences to search for (multi-FASTA format)
      --output_folder       Folder to write output files to
      --output_prefix       Prefix to use for output file names

    Optional Arguments:
      --min_identity        Percent identity threshold used for alignment (default: 50)
      --min_coverage        Percent coverage threshold used for alignment (default: 50)
      --ftp_threads         Number of FTP downloads to execute concurrently (default: 25)

    
    Specifing Genomes for Alignment:

    Genomes may be provided for alignment in two different ways, which can be used
    individually or in combination.
    
    In the most straightforward way, genomes are provided in gzip-compressed FASTA format
    with a wildcard expression (e.g. 'genomes/*.fasta.gz'). When a single wildcard
    character ('*') is used, all files are considered within the specified folder.
    When a double wildcard character ('**') is used, all subdirectories are also
    traversed to find matches for the wildcard expression. One or more set of paths
    can be provided to --genomes using a comma to delimit multiple paths, for example:
        --genomes local_genomes/*.fasta.gz,other_genomes/**.fna.gz

    In addition, you may import genomes directly from the NCBI Prokaryotic Genome Browser
    found at (https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/). After selecting
    your genomes of interest, click on the "Download" button to save a CSV listing all
    of the genomes for alignment. That CSV file may also be specified with the --genome_tables
    flag. More than one table of genomes may be specified using the comma delimiter as above.

    """.stripIndent()
}

