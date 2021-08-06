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
params.min_coverage = 90
params.min_identity = 90
params.aligner = 'diamond'
params.ftp_threads = 25
params.query_gencode = 11
params.max_evalue = 0.001
params.culling_limit = 5
params.max_target_seqs = 100000
params.annotate_geneshot = false
params.abundances_geneshot = false
params.max_n_genes_train_pca = 100000
params.max_pcs_tsne = 50
params.sketch_size = 10000
params.genome_distances = false
params.ani_thresholds = "99,95,90,80,70,60,50"
params.cluster_similarity = 0.9
params.cluster_coverage = 0.9
params.marker_genes = false
params.min_marker_coverage = 50


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
    makedb_blast as makedb_markers;
    align_blast;
    align_blast as align_markers;
    align_diamond;
    makedb_diamond;
    add_genome_name;
    concatenate_alignments;
    concatenate_annotations;
    order_genes;
    generate_gene_map;
    annotate_genes;
    annotate_genes_with_abundances;
    extract_markers;
    reorganize_markers;
    combine_markers;
    calc_marker_distances;
    cluster_genomes;
    aggregate_results;
} from './modules' params(
    output_folder: params.output_folder,
    ftp_output_folder: "${params.output_folder}/genomes",
    output_prefix: params.output_prefix,
    min_identity: params.min_identity,
    min_coverage: params.min_coverage,
    min_marker_coverage: params.min_marker_coverage,
    ftp_threads: params.ftp_threads,
    query_gencode: params.query_gencode,
    max_evalue: params.max_evalue,
    culling_limit: params.culling_limit,
    max_target_seqs: params.max_target_seqs,
    max_n_genes_train_pca: params.max_n_genes_train_pca,
    max_pcs_tsne: params.max_pcs_tsne,
    sketch_size: params.sketch_size,
    ani_thresholds: params.ani_thresholds,
    parse_genome_csv_suffix: "_genomic.fna.gz",
    publishFTP: 'false'
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
      --marker_genes       Optionally provide marker genes in FASTA format which will be used
                            in addition to ANI to estimate a phylogeny for the provided genomes.
                            More than one FASTA file can be specified with comma delimiters.
                            Genes should be provided in amino acid sequence.
      --min_marker_coverage Minimum percent coverage required to use the aligned marker sequence
                            from a particular genome (default: 50)
      --aligner             Alignment algorithm to use (options: diamond, blast; default: diamond)
      --min_identity        Percent identity threshold used for alignment (default: 90)
      --min_coverage        Percent coverage threshold used for alignment (default: 90)
      --ftp_threads         Number of FTP downloads to execute concurrently (default: 25)
      --query_gencode       Genetic code to use for conceptual translation (default: 11)
      --max_evalue          Maximum E-value for any alignment (default: 0.001)
      --culling_limit       If the query range of a hit is enveloped by that of at least
                            this many higher-scoring hits, delete the hit (default: 5, for BLAST)
      --max_target_seqs     Maximum number of alignments to keep, per genome (default: 100000)
      --annotate_geneshot   Optionally format annotations from the geneshot pipeline in a format
                            which can be easily loaded into the gig-map visualization app.
                            This flag does not require the use of --abundances_geneshot (below).
                            The expected file is the output of geneshot named *.results.hdf5.
      --abundances_geneshot Optionally annotate genes using their observed relative abundances
                            from a set of metagenomic datasets, using the output of the geneshot
                            analysis pipeline. This flag can only be used in combination with
                            --annotate_geneshot.
                            The expected file is the output of geneshot named *.details.hdf5.
      --max_n_genes_train_pca
                            The maximum number of genes used to train the PCA model used
                            for ordering genes based on the similarity of the genomes
                            which they align to (default: 100000)
                            NOTE: If this dataset has fewer genes than this, the ordering of
                            genes by similarity of genome alignments will be performed more 
                            precisely by linkage clustering.
      --max_pcs_tsne        The maximum number of dimensions from the PCA output to use
                            for ordering genes by 1-dimensional t-SNE (default: 50)
      --sketch_size         Sketch size (see mash documentation) (default: 10000)
      --genome_distances    If the pairwise genome distances have already been computed,
                            you can directly import them into the analysis instead of repeating
                            that time-consuming process. This flag should be used to indicate
                            the distances.csv.gz file produced for this exact set of genomes.
      --ani_thresholds      Comma-delimited list of percent identity thresholds at which
                            genomes will be clustered for more rapid visualization.
                            Default: 99,95,90,80,70,60,50

    
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


    Subcommands:

    To download a set of genomes from NCBI without performing any alignment, run the command
        nextflow run FredHutch/gig-map -entry download --help

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
                .map({it.trim()})
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

        // If the metagenomic abundances were also provided
        if (params.abundances_geneshot){

            // Format the gene annotations and abundances as a CSV
            annotate_genes_with_abundances(
                Channel
                    .fromPath(
                        params.annotate_geneshot
                    ),
                Channel
                    .fromPath(
                        params.abundances_geneshot
                    )
            )

        }else{

            // Format the gene annotations (only) as a CSV
            annotate_genes(
                Channel
                    .fromPath(
                        params.annotate_geneshot
                    )
            )

        }

    }

    // Cluster genomes by ANI
    cluster_genomes(
        concatenate_alignments.out.combine(
            genome_distances_csv
        ).combine(
            Channel.of(
                params.ani_thresholds.split(/,/)
            )
        )
    )

    // If the user specified any marker sequences
    if (params.marker_genes){

        // Build a BLAST database with the marker sequences
        makedb_markers(
            Channel.fromPath(
                params.marker_genes.split(",").toList()
            )
            .toSortedList()
        )

        // Align the genomes against those markers
        align_markers(
            makedb_markers.out,
            all_genomes
        )

        // Extract the aligned regions for each marker
        extract_markers(
            align_markers.out.join(
                all_genomes.map({
                    it -> [it.name, it]
                })
            )
        )

        // The output of extract_markers has one file per genome, with the
        // FASTA headers indicating the marker of origin

        // Next we will reformat the markers to have one file per marker,
        // with the FASTA headers indicating the genome of origin
        reorganize_markers(
            extract_markers.out.toSortedList()
        )

        // Run the MSA
        combine_markers(
            reorganize_markers.out.flatten()
        )

        // Generate the distance matrices
        calc_marker_distances(
            combine_markers.out
        )

    }

    // Group together all results into a single HDF5 file object
    aggregate_results(
        concatenate_alignments.out,
        order_genes.out,
        genome_distances_csv,
        generate_gene_map.out,
        cluster_genomes.out.toSortedList()
    )

}

// Workflow dedicated to download only

include {
    fetchFTP as saveFTP;
    concatenate_annotations as join_annotations;
} from './modules' params(
    output_prefix: 'downloaded',
    output_folder: params.output_folder,
    ftp_output_folder: "${params.output_folder}/genomes",
    publishFTP: 'true',
    ftp_threads: params.ftp_threads,
)

// Function which prints help message text
def downloadHelpMessage() {
    log.info"""
    Dedicated genome download utility.
    Downloads a set of genomes from NCBI and writes the files to a local folder.

    Genomes in FASTA format will be written to the genomes/ subdirectory,
    while genome metadata will be written to downloaded.genome.annotations.csv.gz.

    Usage:

    nextflow run FredHutch/gig-map -entry download <ARGUMENTS>

    Required Arguments:
      --genome_tables       Tables of NCBI genomes to analyze (see note below)
      --output_folder       Folder to write genome files

    Optional Arguments:
      --ftp_threads         Number of FTP downloads to execute concurrently (default: 25)


    Specifing Genomes for Download:

    Genomes are selected for download directly from the NCBI Prokaryotic Genome Browser
    found at (https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/). After selecting
    your genomes of interest, click on the "Download" button to save a CSV listing all
    of the genomes for alignment. That CSV file must be specified with the --genome_tables
    flag. More than one table of genomes may be specified using a comma delimiter.

    """.stripIndent()
}


workflow download {

    // Show help message if the user specifies the --help flag at runtime
    if (params.help){
        // Invoke the function above which prints the help message
        downloadHelpMessage()
        // Exit out and do not run anything else
        exit 0
    }
    

    // The user must specify an output folder
    if (!params.output_folder){
        log.info"""

        -----------------------
        MISSING --output_folder
        -----------------------

        """.stripIndent()
        downloadHelpMessage()

        // Exit out and do not run anything else
        exit 0
    }

    // The user must specify genomes from NCBI
    if (!params.genome_tables){
        log.info"""

        -----------------------
        MISSING --genome_tables
        -----------------------

        """.stripIndent()
        downloadHelpMessage()

        // Exit out and do not run anything else
        exit 0
    }

    // Parse the set of genomes specified by the user

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
    saveFTP(
        parse_genome_csv
            .out[0]
            .splitText()
            .map({it.trim()})
    )

    // Join the genome annotations from the NCBI table
    join_annotations(
        parse_genome_csv
            .out[1]
            .toSortedList()
    )

}


// Workflow dedicated to gene deduplication only

include {
    parse_genome_csv as parse_genes_csv;
    cdhit;
    fetchFTP as tryFetchFTP;
    annotate_centroids;
} from './modules' params(
    output_folder: params.output_folder,
    ftp_output_folder: "${params.output_folder}/ncbi_genes",
    output_prefix: params.output_prefix,
    cluster_similarity: params.cluster_similarity,
    cluster_coverage: params.cluster_coverage,
    ftp_threads: params.ftp_threads,
    parse_genome_csv_suffix: "_protein.faa.gz",
    skip_missing_ftp: "true",
    publishFTP: 'true'
)

// Function which prints help message text
def deduplicateHelpMessage() {
    log.info"""
    Dedicated gene deduplication utility
    
    Clusters a collection of gene sequences by amino acid similarity and outputs
    the centroids of each cluster.

    Output files:
        clustered.genes.fasta.gz - Amino acid sequences of the centroids of each cluster
        clustered.genes.csv.gz - Table with annotations for each cluster
        clustered.membership.csv.gz - Table linking inputs to each cluster

    Usage:

    nextflow run FredHutch/gig-map -entry deduplicate <ARGUMENTS>

    Required Arguments:
      --genome_tables       Tables of NCBI genomes to analyze (see note below)
      --genes_fasta         Amino acid sequences to search for (multi-FASTA format)
      --output_folder       Folder to write output

    Optional Arguments:
      --cluster_similarity  Similarity threshold for clustering, ranges 0-1 (default: 0.9)
      --cluster_coverage    Overlap (coverage) threshold for clustering, ranges 0-1 (default: 0.9)
      --ftp_threads         Number of FTP downloads to execute concurrently (default: 25)


    Specifing genes for deduplication:

    Gene sequences can be specified for input in two complementary ways, either using
    one or more tables describing a list of genomes to download directly from NCBI,
    as well as one or more FASTA files containing the amino acid sequences of the genes
    to cluster.

    Genomes can be selected for download directly from the NCBI Prokaryotic Genome Browser
    found at (https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/). After selecting
    your genomes of interest, click on the "Download" button to save a CSV listing all
    of the genomes for alignment. That CSV file must be specified with the --genome_tables
    flag. More than one table of genomes may be specified using a comma delimiter.

    NOTE: Not all genomes in NCBI have gene sequences provided. Any genomes in the table
    provided by the user which are missing gene sequences will be ignored.
    """.stripIndent()
}


workflow deduplicate {

    // Show help message if the user specifies the --help flag at runtime
    if (params.help){
        // Invoke the function above which prints the help message
        deduplicateHelpMessage()
        // Exit out and do not run anything else
        exit 0
    }
    

    // The user must specify an output folder
    if (!params.output_folder){
        log.info"""

        -----------------------
        MISSING --output_folder
        -----------------------

        """.stripIndent()
        downloadHelpMessage()

        // Exit out and do not run anything else
        exit 0
    }

    // The user must specify genomes from NCBI or genes in FASTA format
    if (!params.genome_tables && !params.genes_fasta){
        log.info"""

        -----------------------------------------
        MISSING --genome_tables and --genes_fasta
        -----------------------------------------

        """.stripIndent()
        downloadHelpMessage()

        // Exit out and do not run anything else
        exit 0
    }

    // If the user provided a table of genomes from NCBI
    if (params.genome_tables != false){

        // Set up a channel with those files
        Channel
            .fromPath(
                params.genome_tables.split(",").toList()
            )
            .set {
                genome_manifests
            }
        
        // Read the contents of each manifest file
        parse_genes_csv(
            genome_manifests
        )

        // Download each of the files
        tryFetchFTP(
            parse_genes_csv
                .out[0]
                .splitText()
                .map({it.trim()})
        )

        // Set the channel with those files to a channel
        genes_from_ncbi = tryFetchFTP.out
        
    // Otherwise, if --genome_tables was not set
    }else{

        // Make an empty channel with the same name
        Channel
            .empty()
            .set {
                genes_from_ncbi
            }

    }

    // If the user provided a list of FASTAs with genes in amino acid sequence
    if (params.genes_fasta != false){

        // Set up a channel with those files
        Channel
            .fromPath(
                params.genes_fasta.split(",").toList()
            )
            .set {
                genes_from_fasta
            }
        
    // Otherwise, if --genes_fasta was not set
    }else{

        // Make an empty channel with the same name
        Channel
            .empty()
            .set {
                genes_from_fasta
            }

    }

    // Run CD-HIT on all of the gene sequences
    cdhit(
        genes_from_ncbi
            .mix(genes_from_fasta)
            .toSortedList()
    )

    // Generate a simple annotation file for each centroid
    annotate_centroids(
        cdhit.out[0]
    )

}