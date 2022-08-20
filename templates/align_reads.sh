#!/bin/bash

set -euo pipefail

# Keep an integer value keeping track of which input file we are working on
index=0

# Loop over each file in the input folder
for fp in input/*.fastq.gz; do

    # Make sure that the file exists
    [[ -s ${fp} ]] || echo "File not found: $fp"

    # Increment the counter
    let "index+=1"

    # Rename each read header before aligning
    gunzip -c $fp | \
        awk "{if(NR % 4 == 1){print \$1 \"_$index\"  }else{if(NR % 4 == 3){print \"+\"}else{print}}}" | \
        gzip -c \
        > query.fastq.gz

    # Align the reads
    diamond \
        blastx \
        --query query.fastq.gz \
        --out $index.aln.gz \
        --threads !{task.cpus} \
        --db !{refdb} \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
        --min-score !{params.min_score_reads} \
        --query-cover !{params.min_coverage_reads} \
        --id !{params.min_pctiden_reads} \
        --top !{params.top_pct_reads} \
        --block-size !{task.memory.toMega() / (1024 * 6 * task.attempt)} \
        --query-gencode !{params.query_gencode} \
        --compress 1 \
        --unal 0 \
        2> !{sample_name}.$index.log

    # Add the output file to the aggregate, and then delete the shard
    cat $index.aln.gz >> !{sample_name}.aln.gz
    rm $index.aln.gz
done
