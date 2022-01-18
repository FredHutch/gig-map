// Docker containers reused across processes
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:4a6179f"
container__mashtree = "quay.io/hdc-workflows/mashtree:1.2.0"
container__blast = "quay.io/biocontainers/blast:2.11.0--pl526he19e7b1_0"
container__diamond = "quay.io/biocontainers/diamond:2.0.11--hdcc8f71_0"
container__cdhit = "quay.io/biocontainers/cd-hit:4.8.1--h2e03b76_5"
container__clustal = "quay.io/hdc-workflows/clustalo:main"
container__emboss = "biocontainers/emboss:v6.6.0dfsg-7b1-deb_cv1"
container__raxml = "quay.io/biocontainers/raxml-ng:1.0.3--h32fcf60_0"

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

    script:
    template "parse_genome_csv.py"

}


// Remove any trailing whitespaces from genome FASTAs
process clean_genomes {
    container "${container__pandas}"
    label "io_limited"

    input:
        tuple val(output_file_name), path("INPUT.${output_file_name}")
    
    output:
        path "${output_file_name}"

    script:
    template "clean_genomes.sh"
    
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
        path "*" optional params.skip_missing_ftp == "true"

    script:
    template "fetchFTP.py"
    
}


process mash_sketch {
    container "${container__mashtree}"
    label 'io_limited'
    
    input:
        path fasta
    
    output:
        path "${fasta.name}.msh"
    
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
        path "inputs/*"
    
    output:
        path "combined.msh"
    
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
        tuple path(query_msh), path(combined_msh)
    
    output:
        path "${query_msh.name.replaceAll(/.msh/, '.tsv.gz')}"
    
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
    
    script:
    template "filter_genes.py"

}

process aggregate_distances {
    container "${container__pandas}"
    label 'io_limited'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true

    
    input:
        path "inputs/*"
    
    output:
        path "distances.csv.gz"
    
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
        path "inputs/*"

    output:
        path "database.tar"

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
        path database_tar
        path query_fasta

    output:
        tuple val("${query_fasta.name}"), path("alignments.gz")

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
    path refdb
    path query_fasta
    
    output:
    tuple val("${query_fasta.name}"), path("alignments.gz")

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
        path dmnd
    
    output:
        path "${dmnd}.fasta.gz"

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
    path fasta

    output:
    path "database.dmnd"

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
        tuple val(genome_name), path(alignments), path(genome)

    output:
        path "*.markers.fasta.gz" optional true

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
        tuple val(query_name), path("unfiltered.alignments.gz")

    output:
        tuple val("${query_name}"), path("alignments.gz") optional true

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
        path "alignments.gz"

    output:
        path "markers.txt"

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
        tuple val(query_name), path(alignments_tsv_gz), path(gene_list_txt)

    output:
        tuple val("${query_name}"), path("${alignments_tsv_gz}.subset.tsv.gz")

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
        path input_fasta

    output:
        path "*.gz"

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
    publishDir "${params.output_folder}/genes/", mode: 'copy', overwrite: true, enabled: "${params.publishGenes}" == "true"

    input:
        path "fastas_by_genome/*.markers.fasta.gz"

    output:
        path "fastas_by_marker/*.markers.fasta.gz"

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
        path unaligned_fasta

    output:
        path "${unaligned_fasta}.distmat", emit: distmat
        path "${unaligned_fasta}.msa", emit: msa

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
    --outfmt=fasta

"""
}

// Build an ML tree from the MSA
process raxml {
    container "${container__raxml}"
    label 'mem_veryhigh'
    publishDir "${params.output_folder}/markers/", mode: 'copy', overwrite: true

    input:
        path aln_fasta
    
    output:
        path "${aln_fasta.name}.raxml.bestTree"
        path "${aln_fasta.name}.raxml.bestModel"
        path "${aln_fasta.name}.raxml.log"
    
    script:
    template 'raxml.sh'
}


// Add the query genome file name as the last column to the alignments
process add_genome_name {
    container "${container__pandas}"
    label 'io_limited'
    
    input:
    tuple val(genome_name), path(alignments_gz)

    output:
    path "alignments.named.gz"

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
    path "inputs/*.tsv.gz"

    output:
    path "alignments.csv.gz"

"""#!/bin/bash

set -Eeuo pipefail

echo "${params.aln_fmt} genome" | tr ' ' ',' > "alignments.csv"

gunzip -c inputs/*.tsv.gz \
    | tr '\t' ',' \
    >> "alignments.csv"

gzip "alignments.csv"

"""

}

// Combine all of the genome annotations into a single file
process concatenate_annotations {
    container "${container__pandas}"
    label 'io_limited'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    path "annotations/*.csv.gz"

    output:
    path "genome.annotations.csv.gz"

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
    "genome.annotations.csv.gz",
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
    path alignments_csv_gz

    output:
    path "gene_order.txt.gz"

"""#!/bin/bash

set -Eeuo pipefail

order_genes.py \
    "${alignments_csv_gz}" \
    "gene_order.txt.gz"

"""

}


// Generate 2-dimensional t-SNE coordinates for genes based on their alignment to genomes
process generate_gene_map {
    container "${container__pandas}"
    label 'mem_medium'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    path alignments_feather

    output:
    path "tsne.coords.csv.gz"

"""#!/bin/bash

set -Eeuo pipefail

generate_gene_map.py \
    "${alignments_feather}" \
    ${params.max_n_genes_train_pca} \
    ${params.max_pcs_tsne} \
    "tsne.coords.csv.gz" \
    ${task.cpus}

"""

}


// Format a set of gene annotations from the geneshot pipeline as a flat CSV
process annotate_genes {
    container "${container__pandas}"
    label 'mem_medium'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    path geneshot_results_hdf

    output:
    path "gene_annotations.csv.gz"

"""#!/bin/bash

set -Eeuo pipefail

format_geneshot_annotations.py \
    --input "${geneshot_results_hdf}" \
    --output "gene_annotations.csv.gz"

"""

}

// Format a set of gene annotations from the geneshot pipeline as a flat CSV
process annotate_genes_with_abundances {
    container "${container__pandas}"
    label 'mem_medium'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    path geneshot_results_hdf
    path geneshot_details_hdf

    output:
    path "gene_annotations.csv.gz"

"""#!/bin/bash

set -Eeuo pipefail

format_geneshot_annotations.py \
    --input "${geneshot_results_hdf}" \
    --details "${geneshot_details_hdf}" \
    --output "gene_annotations.csv.gz"

"""

}


// Cluster genomes by ANI
process cluster_genomes {
    container "${container__pandas}"
    label 'mem_medium'
   
    input:
    tuple path(alignments_csv_gz), path(dists_csv_gz), val(ani_threshold)

    output:
    path "*.hdf5"

"""#!/bin/bash

set -Eeuo pipefail

cluster_genomes.py \
    --alignments "${alignments_csv_gz}" \
    --dists "${dists_csv_gz}" \
    --ani-threshold ${ani_threshold}

"""

}

// Group together all results into a single HDF5 path object
process aggregate_results {
    container "${container__pandas}"
    label 'mem_medium'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    path alignments_csv_gz
    path gene_order_txt_gz
    path dists_csv_gz
    path tsne_coords_csv_gz
    path "genome_clusters/*"
    path "marker_clusters/*"
    path "marker_distances/*"

    output:
    path "gigmap_output.rdb"

    script:
    template "aggregate_results.sh"

}


// Cluster genomes by ANI
process cdhit {
    container "${container__cdhit}"
    label 'mem_medium'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    path "input.genes.*.fasta.gz"
    
    output:
    path "clustered.genes.fasta.gz"
    path "clustered.membership.csv.gz"
    
    script:
    template "cdhit.sh"
}

// Generate a simple annotation path for each centroid
process annotate_centroids {
    container "${container__pandas}"
    label 'io_limited'
    publishDir "${params.output_folder}", mode: 'copy', overwrite: true
   
    input:
    path "clustered.genes.fasta.gz"
    
    output:
    path "clustered.genes.csv.gz"
    
    script:
    template "annotate_centroids.py"
}
