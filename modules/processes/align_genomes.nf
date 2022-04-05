// Align each sample against the reference database of genes using DIAMOND
process makedb_diamond {
    container "${params.container__diamond}"
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

// Remove any trailing whitespaces from genome FASTAs
process clean_genomes {
    container "${params.container__pandas}"
    label "io_limited"

    input:
        tuple val(output_file_name), path("INPUT.${output_file_name}")
    
    output:
        path "${output_file_name}"

    script:
    template "clean_genomes.sh"
    
}

// Align a genome to a set of genes with the DIAMOND algorithm
process align_diamond {
    container "${params.container__diamond}"
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

// Filter the alignments
process filter_alignments {
    container "${params.container__pandas}"
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

// Make a BLAST database
process makedb_blast {
    container "${params.container__blast}"
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

// Align a genome to a set of genes with the BLAST algorithm
process align_blast {
    container "${params.container__blast}"
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

// Add the query genome file name as the last column to the alignments
process add_genome_name {
    container "${params.container__pandas}"
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
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}", mode: 'copy', overwrite: true
   
    input:
    path "inputs/*.tsv.gz"

    output:
    path "genomes.aln.csv.gz"

"""#!/bin/bash

set -Eeuo pipefail

echo "${params.aln_fmt} genome" | tr ' ' ',' > "genomes.aln.csv"

gunzip -c inputs/*.tsv.gz \
    | tr '\t' ',' \
    >> "genomes.aln.csv"

gzip "genomes.aln.csv"

"""

}

// Select a set of marker genes
process select_markers {
    container "${params.container__pandas}"
    label "io_limited"

    input:
        path "alignments.gz"
        path "centroids.faa.gz"

    output:
        path "markers.fasta.gz"

"""#!/bin/bash

set -e

select_markers.py \
    --input alignments.gz \
    --all-genes centroids.faa.gz \
    --output markers.fasta.gz \
    --max-n "${params.pick_marker_genes}" \
    --min-coverage "${params.min_marker_coverage}" \
    --aln-fmt "${params.aln_fmt}"
"""
}


// Subset a set of alignments to a particular set of genes
process subset_alignments_by_genes {
    container "${params.container__pandas}"
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

// Order genes based on their alignment to genomes
process order_genes {
    container "${params.container__pandas}"
    label 'mem_medium'
    publishDir "${params.output}", mode: 'copy', overwrite: true
   
    input:
    path alignments_csv_gz

    output:
    path "genomes.gene.order.txt.gz"

"""#!/bin/bash

set -euo pipefail

order_genes.py \
    "${alignments_csv_gz}" \
    "genomes.gene.order.txt.gz"

"""

}

// Extract the marker sequence from a genome based on the alignment information
process extract_genes {
    container "${params.container__pandas}"
    label "io_limited"

    input:
        tuple val(genome_name), path(alignments), path(genome)

    output:
        path "*.fasta.gz" optional true

"""#!/bin/bash

set -e

extract_genes.py \
    "${alignments}" \
    "${genome}" \
    "${genome_name}" \
    "${params.aln_fmt}" \
    "${params.min_marker_coverage}"
"""
}

// Reorganize a set of FASTA sequences across files to be grouped by header
process reorganize_fastas {
    container "${params.container__pandas}"
    label "io_limited"
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:
        path "inputs/*.fasta.gz"
        val subfolder

    output:
        path "*.fasta.gz"

"""#!/bin/bash

set -e

reorganize_fastas.py
"""
}
