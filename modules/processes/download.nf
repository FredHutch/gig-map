// Fetch a file via FTP
process fetchFTP {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.ftp_output_folder}", mode: 'copy', overwrite: true, enabled: "${params.publishFTP}" == "true"

    maxForks params.ftp_threads

    input:
        val ftp_url
    
    output:
        path "*" optional params.skip_missing_ftp == "true"

    script:
    template "fetchFTP.py"
    
}


// Parse the NCBI Genome Browser CSV 
process parse_genome_csv {
    container "${params.container__pandas}"
    label "io_limited"

    input:
        path "input.csv"
    
    output:
        path "url_list.txt"
        path "genome_annotations.csv.gz"

    script:
    template "parse_genome_csv.py"

}


// Combine all of the genome annotations into a single file
process concatenate_annotations {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.project_folder}/downloaded_${output_prefix}", mode: 'copy', overwrite: true
   
    input:
    path "annotations/*.csv.gz"
    val output_prefix

    output:
    path "${output_prefix}.annot.csv.gz"

"""#!/usr/bin/env python3
import pandas as pd
import os

# Set the name of the output file based on the output_prefix value (genes/genomes)
output_fp = "${output_prefix}.annot.csv.gz"

# If there are any annotations available
if os.path.exists('annotations'):

    # Read in all of the inputs
    df = pd.concat([
        pd.read_csv(
            os.path.join('annotations', fp)
        )
        for fp in os.listdir('annotations')
        if fp.endswith('.csv.gz')
    ])

    # Write out the table
    df.to_csv(
        output_fp,
        index=None
    )

# If there are no annotations available
else:

    # Write out an empty file
    with open(output_fp, 'wt') as handle:
        handle.write("${output_prefix}_id\\n")

"""
}