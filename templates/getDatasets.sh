#!/bin/bash
set -e

datasets --version

hasData(){
    # Check to see if an accession has data available
    echo "Checking for data in \$1" >&2
    datasets \
        download \
        genome \
        accession \
        \$1 \
        --no-progressbar \
        --include ${params.ncbi_datasets_type} \
        --preview > preview.json 2> error.txt

    # If there was a download error, raise an error to retry
    if grep -q "Download error" error.txt; then
        return 1
    fi

    # Parse the number of records
    parse_preview_json.py
}

INPUT_ACC=${dataset_acc}
for ACC in \${INPUT_ACC} \${INPUT_ACC/A_/F_} \${INPUT_ACC/F_/A_}; do
    hasData "\$ACC"
    if [ -s accession.has.data ]; then
        echo "Downloading accession \${ACC}"
        datasets download genome accession \${ACC} --no-progressbar --include ${params.ncbi_datasets_type}
        echo "Finished downloading \${ACC}"
        [ -s ncbi_dataset.zip ] || echo "ERROR: No data found for \${ACC}" >&2
        [ -s ncbi_dataset.zip ]
        break
    fi
done
