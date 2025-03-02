#!/bin/bash

set -euo pipefail

date
echo
echo "Running workflow from ${PWD}"
echo

# Log the parameters being used
echo PARAMETERS
echo
cat ._wb/tool/params.json
echo

# Run the workflow
echo Starting workflow
nextflow \
    run \
    "${TOOL_REPO}/build_pangenome.nf" \
    --output "${PWD}" \
    -params-file ._wb/tool/params.json \
    -resume

# If temporary files were not placed in a separate location
if [ -d work ]; then
    # Delete the temporary files created during execution
    echo Removing temporary files
    rm -r work
fi


echo
date
echo Done
