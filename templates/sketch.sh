#!/bin/bash

set -e

# Set up the path to the output by removing .orfs.fasta from the input
OUTPUT_PREFIX="\$(echo "${orfs}" | sed 's/.orfs.fasta//')"

mash \
    sketch \
    -p ${task.cpus} \
    -o "\${OUTPUT_PREFIX}" \
    -s ${params.sketch_size} \
    -k ${params.kmer_size} \
    -a \
    "${orfs}"
