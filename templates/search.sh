#!/bin/bash

set -eu

# Set up the path to the output
OUTPUT_TSV="${query_msh.name.replaceAll(/.msh$/, '')}.${ref_msh.name.replaceAll(/.msh$/, '')}.tsv"

# Append the MASH distance
mash \
    dist \
    -p ${task.cpus} \
    "${ref_msh}" "${query_msh}" >> \
    "\${OUTPUT_TSV}"
