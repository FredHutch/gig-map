#!/bin/bash

set -euo pipefail

date
echo
echo "Running workflow from ${PWD}"
echo

# If the TASK_LIMIT variable is set
if [ -z ${TASK_LIMIT} ]; then

    # Add the 'process.maxForks' parameter
    echo "process.maxForks = ${TASK_LIMIT}" >> nextflow.config

fi

# Log the parameters being used
echo PARAMETERS
echo
cat ._wb/tool/params.json
echo

# Run the workflow
echo Starting workflow
nextflow \
    run \
    "${TOOL_REPO}/test_reads.nf" \
    --output "${PWD}" \
    -params-file ._wb/tool/params.json \
    -resume \
    -profile "${PROFILE}"

# If temporary files were not placed in a separate location
if [ -d work ]; then
    # Delete the temporary files created during execution
    echo Removing temporary files
    rm -r work
fi


echo
date
echo Done
