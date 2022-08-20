#!/bin/bash

set -euo pipefail

# Run the FAMLI algorithm to resolve multi-mapping reads
famli \
    filter \
    --input ${input_aln} \
    --output ${sample_name}.json \
    --threads ${task.cpus} \
    --batchsize ${params.famli_batchsize} \
    --sd-mean-cutoff ${params.famli_sd_mean_cutoff} \
    2> ${sample_name}.log

# Compress the output
gzip ${sample_name}.json
