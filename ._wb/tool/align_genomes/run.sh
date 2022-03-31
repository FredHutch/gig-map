#!/bin/bash

set -euo pipefail

# Run the workflow
nextflow \
    run \
    "${TOOL_REPO}/align_genomes.nf" \
    --genes "${GENES}" \
    --genomes "${GENOMES}/*.gz" \
    --min_coverage "${MIN_COVERAGE}" \
    --min_identity "${MIN_IDENTITY}" \
    --max_evalue "${MAX_EVALUE}" \
    --aligner "${ALIGNER}" \
    --query_gencode "${QUERY_GENCODE}" \
    --max_overlap "${MAX_OVERLAP}" \
    --output "${PWD}" \
    -resume

# Delete the temporary files created during execution
rm -r work
