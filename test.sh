#!/bin/bash

python3 app.py \
    --alignments test_data/output.csv.gz \
    --distances test_data/output.dists.tsv.gz \
    --gene-annotations test_data/GCA_000005845.2_ASM584v2_protein.annotations.csv \
    --genome-annotations test_data/NCBI_Escherichia_coli_genomes.test_subset.annotations.csv
