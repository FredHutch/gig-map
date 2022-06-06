#!/bin/bash

set -e

# Combine multiple sketches into a single file
mash \
    paste \
    combined_genomes \
    inputs/*.msh
