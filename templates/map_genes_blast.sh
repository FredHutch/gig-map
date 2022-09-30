#!/bin/bash

set -e

blastp \
    -query <(gunzip -c queries.fasta.gz) \
    -subject <(gunzip -c refs.fasta.gz) \
    -evalue ${params.max_evalue} \
    -outfmt '6 delim=, ${params.aln_fmt}' \
    -num_threads ${task.cpus} \
    | tr ',' '\\t' \
    | gzip -c \
    > unfiltered_gene_mapping.csv.gz