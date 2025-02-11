#!/bin/bash

set -euo pipefail

# Decompress the input
gunzip -c deduplicated.genes.fasta.gz > input.genes.fasta

# Cluster the inputs
cd-hit \
    -i input.genes.fasta \
    -o centroids.faa \
    -c ${params.cluster_similarity} \
    -aS ${params.cluster_coverage} \
    -T ${task.cpus} \
    -M ${task.memory.toMega()} \
    -p 1 \
    -d 0 \

# Compress the outputs
gzip centroids.faa
gzip centroids.faa.clstr
mv centroids.faa.clstr.gz centroids.membership.csv.gz