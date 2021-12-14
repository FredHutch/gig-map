#!/bin/bash

set -euo pipefail

# Generate a timestamp which will be used to label the report file
TIMESTAMP=$(date "+%Y-%m-%d-%H-%M-%S")

# Launch the workflow
nextflow \
    run \
    FredHutch/gig-map \
    -params-file download.params.json \
    -entry download \
    -r main \
    -with-report download.${TIMESTAMP}.report.html \
    -latest \
    -resume \
    "${@:1}"