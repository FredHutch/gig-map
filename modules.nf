// Docker containers reused across processes
container__wget = "quay.io/fhcrc-microbiome/wget:latest"
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
container__mashtree = "quay.io/hdc-workflows/mashtree:1.2.0"
container__blast = "quay.io/biocontainers/blast:2.11.0--pl526he19e7b1_0"
container__diamond = "quay.io/fhcrc-microbiome/docker-diamond:v2.0.6-biopython"

// Default values for parameters
params.output_folder = 'output'
params.output_prefix = 'output'
params.min_coverage = 50
params.min_identity = 50
params.ftp_threads = 25
params.query_gencode = 11
params.max_evalue = 0.001
params.culling_limit = 5
params.max_target_seqs = 100000
params.aln_fmt = "qseqid sseqid pident length qstart qend qlen sstart send slen"


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


process mashtree {
    container "${container__mashtree}"
    label 'mem_medium'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
    
    input:
        file "inputs/*"
    
    output:
        file "${params.output_prefix}.dnd"
        file "${params.output_prefix}.dists.tsv.gz"
    
"""
#!/bin/bash

set -Eeuo pipefail

# Decompress the genome FASTA files
for f in inputs/*.gz; do
    if [ ! -s \$f ]; then continue; fi
    gunzip -c \$f > \${f%.gz}
done

# Run mashtree
mashtree \
    --outmatrix ${params.output_prefix}.dists.tsv \
    --numcpus ${task.cpus} \
    inputs/* \
    > ${params.output_prefix}.dnd

# Compress the distance matrix
gzip ${params.output_prefix}.dists.tsv

"""
}

process makedb_blast {
    container "${container__blast}"
    label 'io_limited'

    input:
        file "inputs/*"

    output:
        file "database.tar"

    script:
    """#!/bin/bash

set -Eeuo pipefail

# Iterate over the inputs
for f in inputs/*; do

    # If the file exists
    if [ -s \$f ]; then

        # Try to unzip it, otherwise just open it
        gunzip -c \$f || cat \$f

    fi

# Write all of the files to a single file
done > database.fasta

# Make the database
makeblastdb \
    -in database.fasta \
    -dbtype prot

# Tar up the database
tar cvf database.tar database.fasta*
    """
}

process align_blast {
    container "${container__blast}"
    label "mem_medium"

    input:
        file database_tar
        file query_fasta

    output:
        tuple val("${query_fasta.name}"), file("alignments.gz")

"""#!/bin/bash

set -e

# Untar the database
tar xvf ${database_tar}

blastx \
    -query <(gunzip -c ${query_fasta}) \
    -db database.fasta \
    -query_gencode ${params.query_gencode} \
    -evalue ${params.max_evalue} \
    -outfmt 6 ${params.aln_fmt} \
    -culling_limit ${params.culling_limit} \
    -max_target_seqs ${params.max_target_seqs} \
    -num_threads ${task.cpus} \
    | gzip -c \
    > alignments.gz

"""
}


process align_diamond {
    container "${container__diamond}"
    label 'mem_medium'
    
    input:
    file refdb
    file query_fasta
    
    output:
    tuple val("${query_fasta.name}"), file("alignments.gz")

"""#!/bin/bash

set -e

diamond \
    blastx \
    --threads ${task.cpus} \
    --db ${refdb} \
    --out alignments.gz \
    --outfmt 6 ${params.aln_fmt} \
    --query ${query_fasta} \
    --unal 0 \
    --max-target-seqs ${params.max_target_seqs} \
    --evalue ${params.max_evalue} \
    --id ${params.min_identity} \
    --subject-cover ${params.min_coverage} \
    --compress 1 \
    --more-sensitive \
    --query-gencode ${params.query_gencode} \
    --range-culling \
    -F 1 \
    --block-size ${task.memory.toMega() / (1024 * 6 * task.attempt)}

    """

}

process extract_dmnd {
    container "${container__diamond}"
    label 'io_limited'
    
    input:
        file dmnd
    
    output:
        file "${dmnd}.fasta.gz"

"""#!/bin/bash

set -e

diamond \
    getseq \
    --db ${dmnd} \
    --out "${dmnd}.fasta.gz" \
    --compress 1

    """

}

// Align each sample against the reference database of genes using DIAMOND
process makedb_diamond {
    container "${container__diamond}"
    label 'mem_medium'
    
    input:
    file fasta

    output:
    file "database.dmnd"

    """
    set -e
    diamond \
      makedb \
      --in ${fasta} \
      --db database.dmnd \
      --threads ${task.cpus}
    """

}

// Add the query genome file name as the last column to the alignments
process add_genome_name {
    container "${container__pandas}"
    label 'io_limited'
    
    input:
    tuple val(genome_name), file(alignments_gz)

    output:
    file "alignments.named.gz"

"""#!/usr/bin/env python3
import pandas as pd

# Read the table
df = pd.read_csv(
    "${alignments_gz}",
    header=None,
    sep="\\t"
)

# Add the name of the genome FASTA
df = df.assign(genome_name="${genome_name}")

# Write out the table
df.to_csv(
    "alignments.named.gz",
    header=None,
    index=None,
    sep="\\t"
)

"""

}

// Combine all of the outputs into a single file
process concatenate_results {
    container "${container__pandas}"
    label 'io_limited'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    file "inputs/*.tsv.gz"

    output:
    file "${params.output_prefix}.csv.gz"

"""#!/usr/bin/env python3
import pandas as pd
import os

# Define the names for the columns
header = "${params.aln_fmt}".split(" ")
# Add the column for the genome name
header.append("genome")

# Read in all of the inputs
df = pd.concat([
    pd.read_csv(
        os.path.join('inputs', fp),
        sep='\\t',
        header=None,
        names=header
    )
    for fp in os.listdir('inputs')
    if fp.endswith('.tsv.gz')
])

# Write out the table
df.to_csv(
    "${params.output_prefix}.csv.gz",
    index=None
)

"""

}