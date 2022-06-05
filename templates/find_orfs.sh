#!/bin/bash

 set -e

 getorf \
    -sequence <(zcat "${genome}") \
    -outseq "${genome}.orfs.fasta" \
    -table "${params.gencode}" \
    -minsize "${params.minsize}"