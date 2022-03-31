#!/bin/bash

set -euo pipefail

# Run the workflow
nextflow \
    run \
    FredHutch/gig-map/download_genes.nf \
    -r ${REVISION} \
    --genome_csv "${GENOME_CSV}" \
    --output "${PWD}/genes/" \
    -resume \
    -latest

# Delete the temporary files created during execution
rm -r work
