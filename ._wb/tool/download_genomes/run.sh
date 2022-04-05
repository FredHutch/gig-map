#!/bin/bash

set -euo pipefail

# Run the workflow
echo Starting workflow
nextflow \
    run \
    "${TOOL_REPO}/download_genomes.nf" \
    --genome_csv "${GENOME_CSV}" \
    --output "${PWD}/genomes/" \
    -resume

# Delete the temporary files created during execution
echo Removing temporary files
rm -r work

echo Done