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
params.aligner = 'diamond'
params.ftp_threads = 25
params.query_gencode = 11
params.max_evalue = 0.001
params.culling_limit = 5
params.max_target_seqs = 100000

// Import the processes to run in this workflow
include {
    parse_genome_csv;
    fetchFTP;
    mashtree;
    makedb_blast;
    align_blast;
    align_diamond;
    makedb_diamond;
} from './modules' params(
    output_folder: params.output_folder,
    output_prefix: params.output_prefix,
    min_identity: params.min_identity,
    min_coverage: params.min_coverage,
    ftp_threads: params.ftp_threads,
    query_gencode: params.query_gencode,
    max_evalue: params.max_evalue,
    culling_limit: params.culling_limit,
    max_target_seqs: params.max_target_seqs,
)

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
      --aligner             Alignment algorithm to use (options: diamond, blast; default: diamond)
      --min_identity        Percent identity threshold used for alignment (default: 50)
      --min_coverage        Percent coverage threshold used for alignment (default: 50)
      --ftp_threads         Number of FTP downloads to execute concurrently (default: 25)
      --query_gencode       Genetic code to use for conceptual translation (default: 11)
      --max_evalue          Maximum E-value for any alignment (default: 0.001)
      --culling_limit       If the query range of a hit is enveloped by that of at least
                            this many higher-scoring hits, delete the hit (default: 5, for BLAST)
      --max_target_seqs     Maximum number of alignments to keep, per genome (default: 100000, for BLAST)

    
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

    NOTE: All genomes must have a unique filename

    """.stripIndent()
}


workflow {

    // Show help message if the user specifies the --help flag at runtime
    if (params.help){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }

    // The user must specify each of the required arguments
    if (!params.genes || !params.output_folder || !params.output_prefix){
        log.info"""

        -----------------------
        MISSING REQUIRED INPUTS
        -----------------------

        """.stripIndent()
        helpMessage()

        // Exit out and do not run anything else
        exit 0
    }

    // The user must specify genomes either by files or from NCBI
    if (!params.genomes && !params.genome_tables){
        log.info"""

        ---------------------
        MISSING INPUT GENOMES
        ---------------------

        """.stripIndent()
        helpMessage()

        // Exit out and do not run anything else
        exit 0
    }

    // Parse the set of genomes specified by the user

    // Local files
    Channel
        .fromPath(
            params.genomes.split(',').toList()
        ).set {
            local_genomes
        }

    // NCBI genomes
    Channel
        .fromPath(
            params.genome_tables.split(",").toList()
        )
        .set {
            genome_manifests
        }

    // Read the contents of each manifest file
    parse_genome_csv(
        genome_manifests
    )

    // Download each of the files
    fetchFTP(
        parse_genome_csv
            .out
            .splitText()
    )

    // Join together the genomes from both sources
    fetchFTP
        .out
        .mix(local_genomes)
        .set {
            all_genomes
        }

    // Compute whole-genome similarity with mashtree
    mashtree(
        all_genomes.toSortedList()
    )

    // If the user has selected DIAMOND for alignment
    if (params.aligner == "diamond"){

        // Make a DIAMOND database for the input genes
        makedb_diamond(
            Channel
                .fromPath(
                    params.genes.split(",").toList()
                )
                .toSortedList()
        )

        // Align the query genes against the genomes
        align_diamond(
            makedb_diamond.out,
            all_genomes
        )

        // Channel with all alignment results
        alignments_output = align_diamond.out
    }

    // If the user has selected BLAST for alignment
    if (params.aligner == "blast"){

        // Make a BLAST database for the input genes
        makedb_blast(
            Channel
                .fromPath(
                    params.genes.split(",").toList()
                )
                .toSortedList()
        )

        // Align the query genes against the genomes
        align_blast(
            makedb_blast.out,
            all_genomes
        )

        // Channel with all alignment results
        alignments_output = align_blast.out
    }

}