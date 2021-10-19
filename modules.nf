// Docker containers reused across processes
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:4a6179f"
container__mashtree = "quay.io/hdc-workflows/mashtree:1.2.0"
container__blast = "quay.io/biocontainers/blast:2.11.0--pl526he19e7b1_0"
container__diamond = "quay.io/biocontainers/diamond:2.0.11--hdcc8f71_0"
container__cdhit = "quay.io/biocontainers/cd-hit:4.8.1--h2e03b76_5"
container__clustal = "quay.io/hdc-workflows/clustalo:main"
container__emboss = "biocontainers/emboss:v6.6.0dfsg-7b1-deb_cv1"

// Default values for parameters
params.output_folder = 'output'
params.output_prefix = 'output'
params.min_coverage = 50
params.min_marker_coverage = 50
params.pick_marker_genes = 10
params.min_identity = 50
params.ftp_threads = 25
params.query_gencode = 11
params.max_evalue = 0.001
params.max_overlap = 50
params.aln_fmt = "qseqid sseqid pident length qstart qend qlen sstart send slen"
params.max_n_genes_train_pca = 10000
params.max_pcs_tsne = 50
params.sketch_size = 10000
params.parse_genome_csv_suffix = "_genomic.fna.gz"
params.cluster_similarity = 0.9
params.cluster_coverage = 0.9
params.skip_missing_ftp = "false"

// Parse the NCBI Genome Browser CSV 
process parse_genome_csv {
    container "${container__pandas}"
    label "io_limited"

    input:
        path "input.csv"
    
    output:
        path "url_list.txt"
        path "genome_annotations.csv.gz"
    
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
    return f"{ftp_prefix}/{id_str}${params.parse_genome_csv_suffix}"

# Populate a list with formatted annotations for each file
# which will be downloaded from NCBI
annotation_list = []

# Function to format the annotation for each row
def format_annotation(r):

    # Get the name of the file which will be downloaded
    annots = dict(
        genome_id=r["GenBank FTP"].rsplit("/", 1)[-1] + "${params.parse_genome_csv_suffix}"
    )

    # Rename fields
    for k, n in [
        ("#Organism Name", "Organism"),
        ("Strain", "Strain"),
        ("BioSample", "BioSample"),
        ("BioProject", "BioProject"),
        ("Assembly", "Assembly"),
        ("Size(Mb)", "Size_Mb"),
        ("GC%", "GC_percent"),
        ("CDS", "CDS")
    ]:
        annots[n] = r[k]

    # Format a combined name
    combined_name = annots["Organism"]
    if isinstance(annots["Strain"], str):
        combined_name = combined_name + "(" + annots["Strain"] + ")"
    annots["Formatted Name"] = combined_name + "[" + annots["Assembly"] + "]"

    return annots

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

        # Format the annotations
        annotation_list.append(
            format_annotation(r)
        )

# Write the list to a file
with open("url_list.txt", "w") as handle:

    # Each URL on its own line
    handle.write("\\n".join(url_list))

# Write the annotations to a CSV
pd.DataFrame(annotation_list).set_index("genome_id").to_csv(
    "genome_annotations.csv.gz"
)

"""
}


// Remove any trailing whitespaces from genome FASTAs
process clean_genomes {
    container "${container__pandas}"
    label "io_limited"

    input:
        tuple val(output_file_name), path("INPUT.${output_file_name}")
    
    output:
        path "${output_file_name}"
    
"""#!/bin/bash

set -Eeuo pipefail

# Function to clean input genomes
clean_genome(){
    tr -d '\\r' | sed 's/[ \\t]*\$//'
}

# If the genome is gzip-compressed
if gzip -t "INPUT.${output_file_name}"; then
    gunzip -c "INPUT.${output_file_name}" | clean_genome | gzip -c > "${output_file_name}"
else
    # Otherwise, the file is not compressed
    cat "INPUT.${output_file_name}" | clean_genome > "${output_file_name}"
fi

"""
}


// Fetch a file via FTP
process fetchFTP {
    container "${container__pandas}"
    label 'io_limited'
    publishDir "${params.ftp_output_folder}", mode: 'copy', overwrite: true, enabled: "${params.publishFTP}" == "true"

    maxForks params.ftp_threads

    input:
        val ftp_url
    
    output:
        file "*" optional params.skip_missing_ftp == "true"
    
"""#!/usr/bin/env python3

import shutil
import urllib.request as request
from urllib.error import URLError
from contextlib import closing


remote_path = "${ftp_url}"
local_path = remote_path.rsplit("/", 1)[-1]
print(f"Downloading from {remote_path}")
print(f"Local path is {local_path}")

try:
    with closing(
        request.urlopen(
            remote_path,
            timeout=60
        )
    ) as r:
        with open(local_path, 'wb') as f:
            shutil.copyfileobj(r, f)
except URLError as e:
    print("The URL appears to not be valid")
    if "${params.skip_missing_ftp}" == "true":
        print("Missing FTP paths will be ignored")
    else:
        print("Raising error")
        raise e

"""
}


process mash_sketch {
    container "${container__mashtree}"
    label 'io_limited'
    
    input:
        file fasta
    
    output:
        file "${fasta.name}.msh"
    
"""
#!/bin/bash

set -Eeuo pipefail

mash \
    sketch \
    -p ${task.cpus} \
    -s ${params.sketch_size} \
    "${fasta}"
"""
}

process mash_join {
    container "${container__mashtree}"
    label 'io_limited'
    
    input:
        file "inputs/*"
    
    output:
        file "combined.msh"
    
"""
#!/bin/bash

set -Eeuo pipefail

mash \
    paste \
    combined \
    inputs/*

"""
}

process mash_dist {
    container "${container__mashtree}"
    label 'io_limited'
    
    input:
        tuple file(query_msh), file(combined_msh)
    
    output:
        file "${query_msh.name.replaceAll(/.msh/, '.tsv.gz')}"
    
"""
#!/bin/bash

set -Eeuo pipefail

mash \
    dist \
    -p ${task.cpus} \
    ${combined_msh} \
    ${query_msh} \
| gzip -c \
> "${query_msh.name.replaceAll(/.msh/, '.tsv.gz')}"

"""
}

process filter_genes {
    container "${container__pandas}"
    label 'io_limited'
    
    input:
        path input_fasta
    
    output:
        path "${input_fasta.name}.filtered.fasta"
    
"""#!/usr/bin/env python3

import gzip
import os

input_fp = "${input_fasta}"
assert os.path.exists(input_fp)
output_fp = "${input_fasta.name}.filtered.fasta"
print(f"Reading in {input_fp}, writing to {output_fp}")

min_gene_len = ${params.min_gene_length}
assert isinstance(min_gene_len, int), min_gene_len
print(f"Removing all genes shorter than {min_gene_len}aa")

def parse_fasta(handle):

    header = None
    seq = []
    counter = 0
    
    for line in handle:
        if line[0] == ">":
            if header is not None and len(seq) > 0:
                yield header, "".join(seq)
                counter += 1
            header = line[1:].split(" ")[0].rstrip("\\n")
            seq = []
        else:
            if len(line) > 1:
                seq.append(line.rstrip("\\n"))
                    
    if header is not None and len(seq) > 0:
        yield header, "".join(seq)
        counter += 1

    print(f"Read in a total of {counter:,} FASTA records")


def is_gzip_compressed(fp):
    with gzip.open(fp, 'rt') as handle:
        try:
            handle.read(1)
            return True
        except:
            return False


if is_gzip_compressed(input_fp):
    input_handle = gzip.open(input_fp, "rt")
else:
    input_handle = open(input_fp, "r")

counter = 0
with open(output_fp, 'w') as output_handle:

    for header, seq in parse_fasta(input_handle):

        if len(seq) >= min_gene_len:

            output_handle.write(f">{header}\\n{seq}\\n")
            counter += 1

print(f"Wrote out {counter:,} sequences")

input_handle.close()
print("DONE")

"""
}

process aggregate_distances {
    container "${container__pandas}"
    label 'io_limited'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true

    
    input:
        file "inputs/*"
    
    output:
        file "distances.csv.gz"
    
"""
#!/usr/bin/env python3

import pandas as pd
import os

# Read in all of the distances
df = pd.concat(
    [
        pd.read_csv(
            os.path.join('inputs', fp),
            sep="\\t",
            header=None,
            names=['query', 'ref', 'dist', 'n', 'ratio']
        ).reindex(
            columns=['query', 'ref', 'dist']
        )
        for fp in os.listdir('inputs')
    ]
).pivot_table(
    index='query',
    columns='ref',
    values='dist'
).to_csv(
    "distances.csv.gz"
)

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
    -outfmt '6 delim=, ${params.aln_fmt}' \
    -num_threads ${task.cpus} \
    | tr ',' '\\t' \
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
    --max-target-seqs 100000 \
    --evalue ${params.max_evalue} \
    --id ${params.min_identity} \
    --subject-cover ${params.min_coverage} \
    --compress 1 \
    --more-sensitive \
    --query-gencode ${params.query_gencode} \
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


// Extract the marker sequence from a BLAST alignment
process extract_markers {
    container "${container__pandas}"
    label "io_limited"

    input:
        tuple val(genome_name), file(alignments), file(genome)

    output:
        file "*.markers.fasta.gz" optional true

"""#!/bin/bash

set -e

extract_markers.py \
    "${alignments}" \
    "${genome}" \
    "${genome_name}" \
    "${params.aln_fmt}" \
    "${params.min_marker_coverage}"
"""
}


// Filter the alignments
process filter_alignments {
    container "${container__pandas}"
    label "io_limited"

    input:
        tuple val(query_name), file("unfiltered.alignments.gz")

    output:
        tuple val("${query_name}"), file("alignments.gz") optional true

"""#!/bin/bash

set -e

filter_alignments.py \
    --input unfiltered.alignments.gz \
    --output alignments.gz \
    --max-overlap "${params.max_overlap}" \
    --aln-fmt "${params.aln_fmt}"
"""
}


// Select a set of marker genes
process select_markers {
    container "${container__pandas}"
    label "io_limited"

    input:
        file "alignments.gz"

    output:
        file "markers.txt"

"""#!/bin/bash

set -e

select_markers.py \
    --input alignments.gz \
    --output markers.txt \
    --max-n "${params.pick_marker_genes}" \
    --min-coverage "${params.min_marker_coverage}" \
    --aln-fmt "${params.aln_fmt}"
"""
}


// Subset a set of alignments to a particular set of genes
process subset_alignments_by_genes {
    container "${container__pandas}"
    label "io_limited"

    input:
        tuple val(query_name), file(alignments_tsv_gz), file(gene_list_txt)

    output:
        tuple val("${query_name}"), file("${alignments_tsv_gz}.subset.tsv.gz")

"""#!/bin/bash

set -e

subset_alignments_by_genes.py \
    --input "${alignments_tsv_gz}" \
    --query-list "${gene_list_txt}" \
    --output "${alignments_tsv_gz}.subset.tsv.gz" \
    --aln-fmt "${params.aln_fmt}"

"""
}


// Go from nucleotide sequences to amino acid
// CURRENTLY UNUSED
process translate_markers {
    container "${container__emboss}"
    label "io_limited"

    input:
        file input_fasta

    output:
        file "*.gz"

"""#!/bin/bash

set -e

OUTPUT_PATH=\$(echo "${input_fasta}" | sed 's/.fasta.gz/.fastp/')
echo "Input: ${input_fasta}"
echo "Output: \$OUTPUT_PATH"

transeq \
    -sequence <(gunzip -c ${input_fasta}) \
    -outseq \$OUTPUT_PATH \
    -table ${params.query_gencode}

echo Compressing output
gzip \$OUTPUT_PATH

"""
}


// Reorganize the marker sequence
process reorganize_markers {
    container "${container__pandas}"
    label "io_limited"

    input:
        file "fastas_by_genome/*.markers.fasta.gz"

    output:
        file "fastas_by_marker/*"

"""#!/bin/bash

set -e

mkdir fastas_by_marker
reorganize_markers.py
"""
}


// Make the multiple sequence alignment
process combine_markers {
    container "${container__clustal}"
    label "mem_medium"
    publishDir "${params.output_folder}/markers/", mode: 'copy', overwrite: true

    input:
        file unaligned_fasta

    output:
        file "${unaligned_fasta}.distmat"
        file "${unaligned_fasta}.msa"

"""#!/bin/bash

set -e

# Make a local copy of the decompressed FASTA
gunzip -c ${unaligned_fasta} > input.fasta

# Run Clustal-Omega
clustalo \
    --in input.fasta \
    -t DNA \
    --full \
    --distmat-out="${unaligned_fasta.name}.distmat" \
    --out="${unaligned_fasta.name}.msa" \
    --threads ${task.cpus} \
    --verbose \
    --outfmt=clustal

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

"""#!/bin/bash

set -Eeuo pipefail

gunzip -c "${alignments_gz}" \
    | while read line; do echo -e "\$line\\t${genome_name}"; done \
    | gzip -c \
    > alignments.named.gz

"""

}

// Combine all of the outputs into a single file
process concatenate_alignments {
    container "${container__pandas}"
    label 'io_limited'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    file "inputs/*.tsv.gz"

    output:
    file "${params.output_prefix}.csv.gz"

"""#!/bin/bash

set -Eeuo pipefail

echo "${params.aln_fmt} genome" | tr ' ' ',' > "${params.output_prefix}.csv"

gunzip -c inputs/*.tsv.gz \
    | tr '\t' ',' \
    >> "${params.output_prefix}.csv"

gzip "${params.output_prefix}.csv"

"""

}

// Combine all of the genome annotations into a single file
process concatenate_annotations {
    container "${container__pandas}"
    label 'io_limited'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    file "annotations/*.csv.gz"

    output:
    file "${params.output_prefix}.genome.annotations.csv.gz"

"""#!/usr/bin/env python3
import pandas as pd
import os

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
    "${params.output_prefix}.genome.annotations.csv.gz",
    index=None
)

"""

}

// Order genes based on their alignment to genomes
process order_genes {
    container "${container__pandas}"
    label 'mem_medium'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    file alignments_csv_gz

    output:
    file "${params.output_prefix}.gene_order.txt.gz"

"""#!/bin/bash

set -Eeuo pipefail

order_genes.py \
    "${alignments_csv_gz}" \
    ${params.max_n_genes_train_pca} \
    ${params.max_pcs_tsne} \
    "${params.output_prefix}.gene_order.txt.gz" \
    ${task.cpus}

"""

}


// Generate 2-dimensional t-SNE coordinates for genes based on their alignment to genomes
process generate_gene_map {
    container "${container__pandas}"
    label 'mem_medium'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    file alignments_feather

    output:
    file "${params.output_prefix}.tsne.coords.csv.gz"

"""#!/bin/bash

set -Eeuo pipefail

generate_gene_map.py \
    "${alignments_feather}" \
    ${params.max_n_genes_train_pca} \
    ${params.max_pcs_tsne} \
    "${params.output_prefix}.tsne.coords.csv.gz" \
    ${task.cpus}

"""

}


// Format a set of gene annotations from the geneshot pipeline as a flat CSV
process annotate_genes {
    container "${container__pandas}"
    label 'mem_medium'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    file geneshot_results_hdf

    output:
    file "${params.output_prefix}.gene_annotations.csv.gz"

"""#!/bin/bash

set -Eeuo pipefail

format_geneshot_annotations.py \
    --input "${geneshot_results_hdf}" \
    --output "${params.output_prefix}.gene_annotations.csv.gz"

"""

}

// Format a set of gene annotations from the geneshot pipeline as a flat CSV
process annotate_genes_with_abundances {
    container "${container__pandas}"
    label 'mem_medium'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    file geneshot_results_hdf
    file geneshot_details_hdf

    output:
    file "${params.output_prefix}.gene_annotations.csv.gz"

"""#!/bin/bash

set -Eeuo pipefail

format_geneshot_annotations.py \
    --input "${geneshot_results_hdf}" \
    --details "${geneshot_details_hdf}" \
    --output "${params.output_prefix}.gene_annotations.csv.gz"

"""

}


// Cluster genomes by ANI
process cluster_genomes {
    container "${container__pandas}"
    label 'mem_medium'
   
    input:
    tuple file(alignments_csv_gz), file(dists_csv_gz), val(ani_threshold)

    output:
    file "*.hdf5"

"""#!/bin/bash

set -Eeuo pipefail

cluster_genomes.py \
    --alignments "${alignments_csv_gz}" \
    --dists "${dists_csv_gz}" \
    --ani-threshold ${ani_threshold}

"""

}

// Group together all results into a single HDF5 file object
process aggregate_results {
    container "${container__pandas}"
    label 'mem_medium'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    file alignments_csv_gz
    file gene_order_txt_gz
    file dists_csv_gz
    file tsne_coords_csv_gz
    file "genome_clusters/*"
    file "marker_clusters/*"
    file "marker_distances/*"

    output:
    file "${params.output_prefix}.rdb"

"""#!/bin/bash

set -Eeuo pipefail

# Start a redis server in the background
redis-server \
    --port 6379 \
    --bind 127.0.0.1 \
    --rdbcompression yes \
    --dbfilename "${params.output_prefix}.rdb" \
    --save "" \
    --dir \$PWD &

aggregate_results.py \
    --alignments "${alignments_csv_gz}" \
    --gene-order "${gene_order_txt_gz}" \
    --dists "${dists_csv_gz}" \
    --tnse-coords "${tsne_coords_csv_gz}" \
    --port 6379 \
    --host 127.0.0.1 || \
    redis-cli shutdown  # In case of failure

# Save the redis store
echo "Saving the redis store"
redis-cli save

# Shutdown the redis server
echo "Shutting down the redis server"
redis-cli shutdown

echo "Done"

"""

}


// Cluster genomes by ANI
process cdhit {
    container "${container__cdhit}"
    label 'mem_medium'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    file "input.genes.*.fasta.gz"
    
    output:
    file "clustered.genes.fasta.gz"
    file "clustered.membership.csv.gz"
    
"""
#!/bin/bash

set -e

# Combine all of the files

# Iterate over each of the input files
for f in input.genes.*.fasta.gz; do
    
    # Make sure the input file exits
    if [[ -s \$f ]]; then

        # If the file is compressed
        if gzip -t \$f; then

            # Decompress it to a stream
            gunzip -c \$f

        else

            # Cat to a stream
            cat \$f

        fi

    fi

# Write the stream to a file
done \
    > input.genes.fasta


# Cluster the inputs
cd-hit \
    -i input.genes.fasta \
    -o clustered.genes.fasta \
    -c ${params.cluster_similarity} \
    -aS ${params.cluster_coverage} \
    -T ${task.cpus} \
    -M ${task.memory.toMega()} \
    -p 1 \
    -d 0 \

# Compress the outputs
gzip clustered.genes.fasta
gzip clustered.genes.fasta.clstr
mv clustered.genes.fasta.clstr.gz clustered.membership.csv.gz
"""
}

// Generate a simple annotation file for each centroid
process annotate_centroids {
    container "${container__pandas}"
    label 'io_limited'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    file "clustered.genes.fasta.gz"
    
    output:
    file "clustered.genes.csv.gz"
    
"""#!/usr/bin/env python3

import gzip

fpi = "clustered.genes.fasta.gz"
fpo = "clustered.genes.csv.gz"

# Open both file paths, input and output
with gzip.open(fpi, "rt") as i, gzip.open(fpo, "wt") as o:

    # Write a header line
    o.write("gene_id,combined_name\\n")

    # Write to the output
    o.write(
        # A newline delimited list
        "\\n".join(
            [
                # Where each value is formatted from a line of the input
                # The formatting will 
                    # - remove the leading '>' character
                    # - replace the first ' ' with ','
                    # - and remove the trailing newline, if any
                line[1:].rstrip("\\n").replace(",", " ").replace(" ", ",", 1)
                # Parsed from each line of the input
                for line in i
                # For that subset of lines which start with '>' and which contain a space
                if line[0] == '>' and ' ' in line
            ]
        )
    )
"""
}
