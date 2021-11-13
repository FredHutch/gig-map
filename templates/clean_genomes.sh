#!/bin/bash

set -Eeuo pipefail

# Function to clean input genomes
clean_genome(){
    tr -d '\\r' | sed 's/[ \\t]*\$//'
}

# If the genome is gzip-compressed
if gzip -t "INPUT.${output_file_name}"; then
    gunzip -c "INPUT.${output_file_name}" | clean_genome | gzip -c > "${output_file_name}"
else
    # Otherwise, the file is not compressed
    cat "INPUT.${output_file_name}" | clean_genome > "${output_file_name}"
fi