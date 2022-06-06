#!/bin/bash

set -eu

# Set up the path to the output
OUTPUT_TSV="${genome_uri.replaceAll(/.*\//, '')}.${query_filename}.tsv"

# Record the genome URI and query name in the TSV
echo -n "${genome_uri}\t" > "\${OUTPUT_TSV}"

# Append the MASH distance
mash \
    dist \
    -p ${task.cpus} \
    "${query_msh}" "${genome_msh}" >> \
    "\${OUTPUT_TSV}"
