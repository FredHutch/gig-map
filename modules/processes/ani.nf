process mash_sketch {
    container "${params.container__mashtree}"
    label 'io_limited'
    publishDir "${params.output}/msh/", mode: 'copy', overwrite: true
    
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
    container "${params.container__mashtree}"
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
    container "${params.container__mashtree}"
    label 'io_limited'
    publishDir "${params.output}/tsv/", mode: 'copy', overwrite: true
    
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

process aggregate_distances {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}", mode: 'copy', overwrite: true

    
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