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
    "${TOOL_REPO}/render.nf" \
    --output "${PWD}" \
    --render_options ._wb/tool/params.json \
    -params-file ._wb/tool/params.json \
    -profile "${PROFILE}" \
    -resume

# Delete the temporary files created during execution
echo Removing temporary files
rm -r work

echo
date
echo Done