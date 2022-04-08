#!/bin/bash

set -euo pipefail

date
echo
echo "Running workflow from ${PWD}"
echo

# Run the workflow
echo Starting workflow
nextflow \
    run \
    "${TOOL_REPO}/download_genomes.nf" \
    --output "${PWD}" \
    -params-file ._wb/tool/params.json \
    -resume \
    -profile "${PROFILE}"

# Delete the temporary files created during execution
echo Removing temporary files
rm -r work

echo
date
echo Done
