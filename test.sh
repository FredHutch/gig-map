#!/bin/bash

# Load the data in the background with redis
redis-server \
    --port 6379 \
    --bind 127.0.0.1 \
    --dbfilename output.rdb \
    --dir $PWD/test_data &

# Start the app
python3 app.py \
    --port 6379 \
    --host 127.0.0.1 \
    --gene-annotations test_data/GCA_000005845.2_ASM584v2_protein.annotations.csv \
    --genome-annotations test_data/NCBI_Escherichia_coli_genomes.test_subset.annotations.csv

# Close the redis server
redis-cli shutdown
