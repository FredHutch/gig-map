// Group together all results into a single HDF5 path object
process aggregate_results {
    container "${params.container__pandas}"
    label 'mem_medium'
    publishDir "${params.project_folder}/package/", mode: 'copy', overwrite: true
   
    input:
    path alignments_csv_gz
    path gene_order_txt_gz
    path dists_csv_gz
    path "genome_clusters/*"
    path "marker_clusters/*"
    path "marker_distances/*"

    output:
    path "gigmap.rdb"

    script:
    template "aggregate_results.sh"

}

// Cluster genomes by ANI
process cluster_genomes {
    container "${params.container__pandas}"
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