#!/bin/bash

datasets --version

hasData(){
    # Check to see if an accession has data available
    echo "Checking for data in \$1"
    datasets \
        download \
        genome \
        accession \
        \$1 \
        --no-progressbar \
        --include ${params.ncbi_datasets_type} \
        --preview > preview.json

    # Parse the number of records
    parse_preview_json.py
    if [ -s accession.has.data ]; then
        echo "Data found for \$1"
        return 0
    else
        echo "No data found for \$1"
        return 1
    fi
}

INPUT_ACC=${dataset_acc}
for ACC in \${INPUT_ACC} \${INPUT_ACC/A_/F_} \${INPUT_ACC/F_/A_}; do
    if hasData \$ACC; then
        echo "Downloading accession \${ACC}"
        datasets download genome accession \${ACC} --no-progressbar --include ${params.ncbi_datasets_type}
        break
    fi
done
