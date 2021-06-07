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

// Parse the NCBI Genome Browser CSV 
process parse_genome_csv {
    container "${container__pandas}"
    label "io_limited"

    input:
        path "input.csv"
    
    output:
        path "url_list.txt"
    
"""
#!/usr/bin/env python3

import pandas as pd
import re

df = pd.read_csv("input.csv")

# Make a list of URLs to download
url_list = []

# Function to format the path to the genome FASTA
def format_ftp(ftp_prefix):

    assert isinstance(ftp_prefix, str)
    assert ftp_prefix.startswith("ftp://")
    assert "/" in ftp_prefix
    assert not ftp_prefix.endswith("/")

    # The ID of the assembly is the final directory name
    id_str = ftp_prefix.rsplit("/", 1)[-1]

    # Return the path to the genome FASTA
    return f"{ftp_prefix}/{id_str}_genomic.fna.gz"

# Iterate over each row in the table
for _, r in df.iterrows():
    
    # If there is no value in the 'GenBank FTP' column
    if pd.isnull(r['GenBank FTP']):
    
        # Skip it
        continue
    
    # If the 'GenBank FTP' column doesn't start with 'ftp://'
    elif not r['GenBank FTP'].startswith('ftp://'):
    
        # Skip it
        continue

    # Otherwise
    else:

        # Format the path to the genome in that folder
        url_list.append(
            format_ftp(r['GenBank FTP'])
        )

# Write the list to a file
with open("url_list.txt", "w") as handle:

    # Each URL on its own line
    handle.write("\\n".join(url_list))

"""
}


// Fetch a file via FTP
process fetchFTP {
    container "${container__wget}"
    label 'io_limited'

    maxForks params.ftp_threads

    input:
        val ftp_url
    
    output:
        file "*"
    
"""
#!/bin/bash
set -e

echo "Downloading from ${ftp_url}"

wget --quiet ${ftp_url}

# For any gzip-compressed files
for fp in *.gz; do

    # If any such file exists
    if [ -s \$fp ]; then
        # Make sure that it is appropriately compressed
        echo "Checking that \$fp is gzip-compressed"
        gzip -t \$fp
    fi

done

"""
}

// process align_genes {
//     container "${}"
// }


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

    all_genomes.view()

}