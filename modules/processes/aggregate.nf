// Group together all results into a single HDF5 path object
process aggregate_results {
    container "${params.container__pandas}"
    label 'mem_medium'
    publishDir "${params.output}", mode: 'copy', overwrite: true
   
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

set -euo pipefail

cluster_genomes.py \
    --alignments "${alignments_csv_gz}" \
    --dists "${dists_csv_gz}" \
    --ani-threshold ${ani_threshold}

"""

}

// Create a manifest of the genomes
process create_genome_manifest {
    container "${params.container__pandas}"
    publishDir "${params.output}", mode: 'copy', overwrite: true
    label 'io_limited'
   
    input:
    path "genome_aln.csv.gz"

    output:
    path "genome.manifest.csv"

"""#!/usr/bin/env python3

import os
import pandas as pd

# Read the table of genome alignments
pd.read_csv(
    "genome_aln.csv.gz"
).reindex(
    # Only keep the genome IDs
    columns=["genome"]
).drop_duplicates(
).rename(
    columns=dict(
        genome="genome_id"
    )
).sort_values(
    by="genome_id"
).assign(
    **{
        "Formatted Name": lambda d: d["genome_id"].apply(
            lambda s: s.replace(".gz", "").replace(".fna", "").replace(".fasta", "")
        )
    }
).to_csv(
    "genome.manifest.csv",
    index=None
)

"""

}

// Create a manifest of the genes
process create_gene_manifest {
    container "${params.container__pandas}"
    publishDir "${params.output}", mode: 'copy', overwrite: true
    label 'io_limited'
   
    input:
    path "genome_aln.csv.gz"

    output:
    path "gene.manifest.csv"

"""#!/usr/bin/env python3

import os
import pandas as pd

# Read the table of genome alignments
pd.read_csv(
    "genome_aln.csv.gz"
).reindex(
    # Only keep the gene IDs
    columns=["sseqid"]
).drop_duplicates(
).rename(
    columns=dict(
        sseqid="gene_id"
    )
).sort_values(
    by="gene_id"
).assign(
    **{
        "Formatted Name": lambda d: d["gene_id"].apply(
            lambda s: s.replace(".gz", "").replace(".faa", "").replace(".fasta", "")
        )
    }
).to_csv(
    "gene.manifest.csv",
    index=None
)

"""

}