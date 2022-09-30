#!/bin/bash

set -e

diamond \
    makedb \
    --in refs.fasta.gz \
    --db database.dmnd \
    --threads ${task.cpus}

diamond \
    blastp \
    --threads ${task.cpus} \
    --db database.dmnd \
    --out unfiltered_gene_mapping.csv.gz \
    --outfmt 6 ${params.aln_fmt} \
    --query queries.fasta.gz \
    --unal 0 \
    --max-target-seqs 100000 \
    --evalue ${params.max_evalue} \
    --id ${params.min_identity} \
    --subject-cover ${params.min_coverage} \
    --compress 1 \
    --more-sensitive \
    --block-size ${task.memory.toMega() / (1024 * 6 * task.attempt)}
