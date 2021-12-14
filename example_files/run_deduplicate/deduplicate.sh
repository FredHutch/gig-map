#!/bin/bash

set -euo pipefail

# Generate a timestamp which will be used to label the report file
TIMESTAMP=$(date "+%Y-%m-%d-%H-%M-%S")

# Launch the workflow
nextflow \
    run \
    FredHutch/gig-map \
    -params-file deduplicate.params.json \
    -entry deduplicate \
    -r main \
    -with-report deduplicate.${TIMESTAMP}.report.html \
    -latest \
    -resume \
    "${@:1}"