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

# Set up the input files that are expected in the PANGENOME directory
GENES="${PANGENOME}/gene_catalog/centroids.faa.gz"
GENE_BINS="${PANGENOME}/bin_pangenome/gene_bins.csv"
GROUP_PROFILE="${PANGENOME}/bin_pangenome/group_profile.csv"
GENOME_GROUPS="${PANGENOME}/bin_pangenome/genome_groups.csv"

echo "Pangenome Reference Files:"
for fp in "${GENES}" "${GENE_BINS}" "${GROUP_PROFILE}" "${GENOME_GROUPS}"; do
    if [ ! -f "${fp}" ]; then
        echo "ERROR: Missing file: ${fp}"
        exit 1
    else
        echo "${fp}"
    fi
done

# Run the workflow
echo Starting workflow
nextflow \
    run \
    "${TOOL_REPO}/align_pangenome.nf" \
    --genes "${GENES}" \
    --gene_bins "${GENE_BINS}" \
    --group_profile "${GROUP_PROFILE}" \
    --genome_groups "${GENOME_GROUPS}" \
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
