#!/bin/bash

set -euo pipefail

# Generate a timestamp which will be used to label the report file
TIMESTAMP=$(date "+%Y-%m-%d-%H-%M-%S")

# Launch the workflow
nextflow \
    run \
    FredHutch/gig-map/render.nf \
    -params-file render.params.json \
    -r main \
    -with-report render.${TIMESTAMP}.report.html \
    -latest \
    -resume \
    "${@:1}"