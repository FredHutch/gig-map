#!/bin/bash

set -euo pipefail

# Run the workflow
nextflow \
    run \
    "${TOOL_REPO}/deduplicate.nf" \
    --genes "${GENES}/*" \
    --output "${PWD}" \
    --cluster_similarity ${CLUSTER_SIMILARITY} \
    --cluster_coverage ${CLUSTER_COVERAGE} \
    -resume

# Delete the temporary files created during execution
rm -r work
