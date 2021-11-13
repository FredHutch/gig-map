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