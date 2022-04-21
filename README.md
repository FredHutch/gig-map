# GiG-map (Genes in Genomes - Map)
Build a map of genes present across a set of microbial genomes

[![Test gig-map](https://github.com/FredHutch/gig-map/actions/workflows/test.yaml/badge.svg?branch=main&event=push)](https://github.com/FredHutch/gig-map/actions/workflows/test.yaml)

[![Docker Repository on Quay](https://quay.io/repository/hdc-workflows/gig-map/status "Docker Repository on Quay")](https://quay.io/repository/hdc-workflows/gig-map)

[**Documentation Wiki**](https://github.com/FredHutch/gig-map/wiki)

> Note: The documentation below refers to `gig-map` versions starting with `0.2.0`.
> To run the initial version with a single monolithic workflow entrypoint, use the
> flag `0.1.0` to access that [pre-release](https://github.com/FredHutch/gig-map/releases/tag/0.1.0).

# Background

When trying to understand the biological context of a set of microbial protein-coding
genes, it can be very helpful to understand the distribution of organisms which contain
that coding sequence in their genome. One of the most useful characteristics used to
differentiate those microbial genomes which do or do not contain a set of genes is their
overall evolutionary history, for which taxonomic assignment is commonly used as a proxy.
With this type of analysis it can become more easily evident when a set of genes are,
for example, common to a cohesive taxonomic grouping (i.e. genus-, or species-specific),
found sporadically across representatives within a taxonomic grouping (i.e. strain-specific),
or found across disparate taxonomic groupings (i.e., horizontally transfered genes).

To implement this type of analysis, this repository contains code needed to align a set
of genes against a set of genomes (e.g. BLAST), calculate the similarity of those genomes
to each other (e.g. [mashtree](https://github.com/lskatz/mashtree)), and visualize the
results alongside annotations of the genes and/or genomes.

# Workflow Documentation

For a more complete set of documentation of the workflow commands which can be used to
organize and align microbial genes and genomes, please [visit the Wiki](https://github.com/FredHutch/gig-map/wiki).

# GiG-map Visualization

To visualize the results of the gig-map alignment workflow, the output files may be
rendered using the interactive tool provided as `app.py`. In order to use this tool:

1. Clone this repository locally
2. Create a python virtual environment (tested with Python 3.8 and higher)
3. Install the python dependencies contained in `requirements.txt`
4. Run `python3 app.py --help` for detailed instructions on launching the visualization

# Development Guide

Any contributions to this repository are more than welcome.
Please start by opening an issue to discuss any potential contribution, bugfix, or improvement.

## Testing

To test the `gig-map` toolset, clone the repository locally and execute `run_tests.sh` from `test_data`.
You will need to have Python3, Nextflow, BATS, and Docker configured appropriately for tests to run.

Note: Due to unresolved permissions issues, the `act` utility does not appear to work for testing.
Any suggestions for resolving this issue would be greatly appreciated.
