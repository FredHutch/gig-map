#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.output_folder = false
params.output_prefix = false
params.genomes = false
params.genes_fasta = false
params.genes_dmnd = false
params.genome_tables = false
params.min_coverage = 50
params.min_identity = 50
params.aligner = 'diamond'
params.ftp_threads = 25
params.query_gencode = 11
params.max_evalue = 0.001
params.culling_limit = 5
params.max_target_seqs = 100000
params.annotate_geneshot = false
params.max_n_genes_train_pca = 10000
params.max_pcs_tsne = 50
params.sketch_size = 10000
params.genome_distances = false


// Import the processes to run in this workflow
include {
    parse_genome_csv;
    fetchFTP;
    extract_dmnd;
    mash_sketch;
    mash_join;
    mash_dist;
    aggregate_distances;
    makedb_blast;
    align_blast;
    align_diamond;
    makedb_diamond;
    add_genome_name;
    concatenate_alignments;
    concatenate_annotations;
    order_genes;
    generate_gene_map;
    annotate_genes;
    aggregate_results;
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
    max_n_genes_train_pca: params.max_n_genes_train_pca,
    max_pcs_tsne: params.max_pcs_tsne,
    sketch_size: params.sketch_size,
)

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/gig-map <ARGUMENTS>

    Required Arguments:
      --genomes             Genome sequences in FASTA format (see note below)
      --genome_tables       Tables of NCBI genomes to analyze (see note below)
      --genes_fasta         Amino acid sequences to search for (multi-FASTA format)
      --genes_dmnd          Amino acid sequences to search for (DIAMOND database format *.dmnd)
                            (either --genes_fasta or --genes_dmnd is required)
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
      --max_target_seqs     Maximum number of alignments to keep, per genome (default: 100000)
      --annotate_geneshot   Optionally format annotations from the geneshot pipeline in a format
                            which can be easily loaded into the gig-map visualization app.
                            The expected file is the output of geneshot named *.results.hdf5.
      --max_n_genes_train_pca
                            The maximum number of genes used to train the PCA model used
                            for ordering genes based on the similarity of the genomes
                            which they align to (default: 10000)
      --max_pcs_tsne        The maximum number of dimensions from the PCA output to use
                            for ordering genes by 1-dimensional t-SNE (default: 50)
      --sketch_size         Sketch size (see mash documentation) (default: 10000)
      --genome_distances    If the pairwise genome distances have already been computed,
                            you can directly import them into the analysis instead of repeating
                            that time-consuming process. This flag should be used to indicate
                            the distances.csv.gz file produced for this exact set of genomes.

    
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


    Specifying Genes for Alignment:

    Genes may be provided for alignment either in FASTA format or as a pre-formatted DIAMOND
    database. There is no computational benefit derived from providing a pre-formatted database,
    and so that option is provided merely for convenience. Note that DIAMOND database formats
    depend on the version of the software being used, and so there may be compatibility issues
    when using databases compiled with versions of DIAMOND which are significantly different
    from the version used for alignment in this workflow (v2.0.6 at time of writing).

    To provide genes in FASTA format (gzip-optional), use the --genes_fasta flag.
    To provide genes in DIAMOND database format (v2.0.6 recommended), use the --genes_dmnd flag.

    Multiple files may be specified with either flag, which will be concatenated for alignment.

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
    if (!params.output_folder || !params.output_prefix){
        log.info"""

        -----------------------
        MISSING REQUIRED INPUTS
        -----------------------

        """.stripIndent()
        helpMessage()

        // Exit out and do not run anything else
        exit 0
    }

    // The user must specify genes as either FASTA or DMND
    if (!params.genes_dmnd && !params.genes_fasta){
        log.info"""

        -------------------
        MISSING INPUT GENES
        -------------------

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

    // If one or more local files have been provided
    if(params.genomes){

        // Populate a channel with those local genomes
        Channel
            .fromPath(
                params.genomes.split(',').toList()
            ).set {
                local_genomes
            }

    // If no --genomes flag was used
    }else{

        // Set up an empty channel
        Channel.empty().set{local_genomes}

    }

    // NCBI genomes
    if(params.genome_tables){
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
                .out[0]
                .splitText()
        )

        // Join the genome annotations from the NCBI table
        concatenate_annotations(
            parse_genome_csv
                .out[1]
                .toSortedList()
        )

        // Join together the genomes from both sources
        fetchFTP
            .out
            .mix(local_genomes)
            .set {
                all_genomes
            }
    }else{

        // The only genomes to process will come from the
        // set of local genome
        local_genomes.set{all_genomes}

    }

    // If a file with all pairwise distances were provided by the user
    if (params.genome_distances){

        // Point to that file
        Channel
            .fromPath(params.genome_distances)
            .set { genome_distances_csv }

    // Otherwise
    }else{

        // Compute compressed genome representations (sketches) with mash
        mash_sketch(
            all_genomes
        )

        // Combine all of the sketches
        mash_join(
            mash_sketch.out.toSortedList()
        )

        // Calculate the distance of each genome against the total
        mash_dist(
            mash_sketch.out.combine(mash_join.out)
        )

        // Join together those distances
        aggregate_distances(
            mash_dist.out.toSortedList()
        )

        // Point to the file with all aggregated distances
        genome_distances_csv = aggregate_distances.out

    }

    // If the --genes_fasta flag was specified
    if (params.genes_fasta){

        // Make a channel containing the files from a comma-delimited list
        Channel
            .fromPath(
                params.genes_fasta.split(",").toList()
            )
            .set {
                genes_fasta_ch
            }

    } else {

        // Set up an empty channel
        Channel.empty().set{genes_fasta_ch}

    }

    // If the --genes_dmnd flag was specified
    if (params.genes_dmnd){

        // Make a channel containing the files from a comma-delimited list
        Channel
            .fromPath(
                params.genes_dmnd.split(",").toList()
            )
            .set {
                genes_dmnd_ch
            }

        // Reformat those DMND database files as FASTA
        extract_dmnd(
            genes_dmnd_ch
        )

        // Set up a channel containing the extracted FASTA files
        genes_dmnd_extracted = extract_dmnd.out

    } else {

        // Set up an empty channel
        Channel.empty().set{genes_dmnd_extracted}

    }

    // Combine the genes from FASTA and DMND input files
    genes_fasta_ch
        .mix(genes_dmnd_extracted)
        .set {
            combined_genes_ch
        }

    // If the user has selected DIAMOND for alignment
    if (params.aligner == "diamond"){

        // Make a DIAMOND database for the input genes
        makedb_diamond(
            combined_genes_ch.toSortedList()
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
            combined_genes_ch.toSortedList()
        )

        // Align the query genes against the genomes
        align_blast(
            makedb_blast.out,
            all_genomes
        )

        // Channel with all alignment results
        alignments_output = align_blast.out
    }

    // Add the name of the query genome to the alignments file
    add_genome_name(
        alignments_output
    )

    // Concatenate the results
    concatenate_alignments(
        add_genome_name.out.toSortedList()
    )

    // Order the genes based on the genomes they align to
    order_genes(
        concatenate_alignments.out
    )

    // Generate 2-dimensional t-SNE coordinates for genes based on their alignment to genomes
    generate_gene_map(
        concatenate_alignments.out
    )

    // If a set of geneshot results were provided
    if (params.annotate_geneshot){

        // Format those annotations as a CSV
        annotate_genes(
            Channel
                .fromPath(
                    params.annotate_geneshot
                )
        )
    }

    // Group together all results into a single HDF5 file object
    aggregate_results(
        concatenate_alignments.out,
        order_genes.out,
        genome_distances_csv,
        generate_gene_map.out,
    )

}