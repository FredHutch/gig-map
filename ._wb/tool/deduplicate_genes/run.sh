#!/bin/bash

set -euo pipefail

echo "Running workflow from ${PWD}"

# Run the workflow
echo Starting workflow
nextflow \
    run \
    "${TOOL_REPO}/deduplicate.nf" \
    --output "${PWD}" \
    -params-file ._wb/tool/params.json \
    -resume \
    -profile "${PROFILE}"

# Delete the temporary files created during execution
echo Removing temporary files
rm -r work

echo Done
