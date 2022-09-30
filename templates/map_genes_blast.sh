#!/bin/bash

set -e

blastp \
    -query <(gunzip -c ${query_fasta}) \
    -subject <(gunzip -c ${ref_fasta}) \
    -evalue ${params.max_evalue} \
    -outfmt '6 delim=, ${params.aln_fmt}' \
    -num_threads ${task.cpus} \
    | tr ',' '\\t' \
    | gzip -c \
    > unfiltered_gene_mapping.csv.gz