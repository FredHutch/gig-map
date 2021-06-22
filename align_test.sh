#!/bin/bash

set -e

NXF_VER=21.04.1 \
nextflow \
    run \
    main.nf \
    -c nextflow.config \
    -profile testing \
    -w work/ \
    -with-docker ubuntu:latest \
    --genome_tables test_data/NCBI_Escherichia_coli_genomes.test_subset.csv \
    --genes_fasta test_data/GCA_000005845.2_ASM584v2_protein.faa.gz \
    --output_folder output_diamond_fasta \
    --output_prefix output \
    -resume
