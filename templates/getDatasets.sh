#!/bin/bash

datasets --version

echo "Downloading accession ${dataset_acc}"
datasets download genome accession ${dataset_acc} --no-progressbar --include ${params.ncbi_datasets_type}
