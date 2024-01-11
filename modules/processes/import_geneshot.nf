// Format a set of gene annotations from the geneshot pipeline as a flat CSV
process annotate_genes {
    container "${params.container__pandas}"
    label 'mem_medium'
    publishDir "${params.output}", mode: 'copy', overwrite: true
   
    input:
    path geneshot_results_hdf

    output:
    path "gene_annotations.csv.gz"

"""#!/bin/bash

set -euo pipefail

format_geneshot_annotations.py \
    --input "${geneshot_results_hdf}" \
    --output "gene_annotations.csv.gz"

"""

}

// Format a set of gene annotations from the geneshot pipeline as a flat CSV
process annotate_genes_with_abundances {
    container "${params.container__pandas}"
    label 'mem_medium'
    publishDir "${params.output}", mode: 'copy', overwrite: true
   
    input:
    path geneshot_results_hdf
    path geneshot_details_hdf

    output:
    path "gene_annotations.csv.gz"

"""#!/bin/bash

set -euo pipefail

format_geneshot_annotations.py \
    --input "${geneshot_results_hdf}" \
    --details "${geneshot_details_hdf}" \
    --output "gene_annotations.csv.gz"

"""
}