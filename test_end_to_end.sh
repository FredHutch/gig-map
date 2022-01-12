#!/bin/bash

set -e

OUTPUT_FOLDER=output_end_to_end
OUTPUT_PREFIX=gigmap_alignment
GENOME_TABLES=test_data/NCBI_Escherichia_coli_genomes.test_subset.csv

NXF_CONFIG="-c nextflow.config -profile testing -w work/ -with-docker ubuntu:latest -resume"

# Download genomes
NXF_VER=21.04.1 \
nextflow \
    run \
    ${NXF_CONFIG} \
    main.nf \
    -entry download \
    --genome_tables ${GENOME_TABLES} \
    --output_folder ${OUTPUT_FOLDER}

# Download genes
NXF_VER=21.04.1 \
nextflow \
    run \
    ${NXF_CONFIG} \
    main.nf \
    -entry deduplicate \
    --genome_tables ${GENOME_TABLES} \
    --output_folder ${OUTPUT_FOLDER}

# Align genes against genomes
NXF_VER=21.04.1 \
nextflow \
    run \
    ${NXF_CONFIG} \
    main.nf \
    --genome_tables ${GENOME_TABLES} \
    --genes_fasta ${OUTPUT_FOLDER}/clustered.genes.fasta.gz \
    --output_folder ${OUTPUT_FOLDER} \
    --output_prefix ${OUTPUT_PREFIX}

# Render a plot
NXF_VER=21.04.1 \
nextflow \
    run \
    ${NXF_CONFIG} \
    render.nf \
    --rdb ${OUTPUT_FOLDER}/${OUTPUT_PREFIX}.rdb \
    --gene_annotations ${OUTPUT_FOLDER}/clustered.genes.csv.gz \
    --genome_annotations ${OUTPUT_FOLDER}/downloaded.genome.annotations.csv.gz \
    --output_folder ${OUTPUT_FOLDER} \
    --output_prefix ${OUTPUT_PREFIX} \
    --mem_gbs 4
