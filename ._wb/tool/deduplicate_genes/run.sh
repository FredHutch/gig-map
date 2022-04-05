#!/bin/bash

set -euo pipefail

# Run the workflow
nextflow \
    run \
    "${TOOL_REPO}/deduplicate.nf" \
    --genes "${GENES}/*.gz" \
    --output "${PWD}" \
    --cluster_similarity ${CLUSTER_SIMILARITY} \
    --cluster_coverage ${CLUSTER_COVERAGE} \
    --min_gene_length ${MIN_GENE_LENGTH} \
    -resume

# Delete the temporary files created during execution
rm -r work
