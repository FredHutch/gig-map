#!/bin/bash

set -Eeuo pipefail

# Start a redis server in the background
redis-server \
    --port 6379 \
    --bind 127.0.0.1 \
    --rdbcompression yes \
    --dbfilename "gigmap.rdb" \
    --save "" \
    --dir \$PWD &

aggregate_results.py \
    --alignments "${alignments_csv_gz}" \
    --gene-order "${gene_order_txt_gz}" \
    --dists "${dists_csv_gz}" \
    --tnse-coords "${tsne_coords_csv_gz}" \
    --port 6379 \
    --host 127.0.0.1 || \
    redis-cli shutdown  # In case of failure

# Save the redis store
echo "Saving the redis store"
redis-cli save

# Shutdown the redis server
echo "Shutting down the redis server"
redis-cli shutdown

echo "Done"