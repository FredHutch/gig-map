#!/bin/bash

raxml-ng \
    --all \
    --msa ${aln_fasta} \
    --msa-format FASTA \
    --model ${params.raxml_model} \
    --tree pars{${params.raxml_starting_trees}} \
    --bs-trees ${params.raxml_bs_trees} \
    --threads ${task.cpus} \
    --force perf_threads
