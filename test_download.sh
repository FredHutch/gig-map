#!/bin/bash

set -e

NXF_VER=21.04.1 nextflow run download -c nextflow.config -profile testing -w work/ -with-docker ubuntu:latest --genome_tables test_data/NCBI_Escherichia_coli_genomes.test_subset.csv --output_folder download_genomes -resume