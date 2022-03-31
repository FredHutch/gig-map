#!/bin/bash

set -euo pipefail

# Run the workflow
nextflow \
    run \
    "${TOOL_REPO}/download_genes.nf" \
    --genome_csv "${GENOME_CSV}" \
    --output "${PWD}/genes/" \
    -resume

# Delete the temporary files created during execution
rm -r work
