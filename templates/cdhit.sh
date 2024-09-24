#!/bin/bash

set -euo pipefail

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
    | tr '|' '_' \
    | deduplicate_gene_names.sh \
    > input.genes.fasta


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