# GiG-map (Genes in Genomes - Map)
Build a map of genes present across a set of microbial genomes

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

# Implementation

The execution of the analysis encoded by this repository is roughly grouped into two
activities: 1) alignment of genes against genomes and 2) visualization of those alignment
results. 

## Alignment of Genes Against Genomes

To align a set of genomes against a set of genomes, while also comparing those genomes
to each other, this repository may be run as a Nextflow workflow. For more details on
how to run Nextflow on your system, take a look at [their documentation](https://nextflow.io/).
For a detailed set of instructions on how to specify your genes and genomes for analysis,
run `nextflow run FredHutch/gig-map --help` to display the associated help message.

## Download Bacterial Genomes

It can be useful to decouple the process of downloading reference genomes from the
work it takes to align genes against those genomes. Similarly, we often want to download
a set of genomes and then align multiple distinct sets of genes against those genomes.
For those cases, it can be very helpfult to download the genomes in one step with a dedicated
utility. To view the instructions for running this utility,
run `nextflow run FredHutch/gig-map -entry download --help`

## Dedupliate Bacterial Genes

When aligning a set of genes, it can be helpful to first combine any genes which have a
very similar sequence. We refer to that particular process as "deduplicating" a set of genes.
To deduplicate a set of genes based on a threshold of amino acid similarity and overlap,
there is a standalone utility provided along with `gig-map`. To view the instructions for
running this utility, run `nextflow run FredHutch/gig-map -entry deduplicate --help`.
This utility uses the CD-HIT algorithm under the hood to cluster coding sequences
(with support to download directly from NCBI via FTP) and save the combined gene
catalog for downstream alignment.

## GiG-map Visualization

To visualize the results of the gig-map alignment workflow, the output files may be
rendered using the interactive tool provided as `app.py`. In order to use this tool:

1. Clone this repository locally
2. Create a python virtual environment (tested with Python 3.8 and higher)
3. Install the python dependencies contained in `requirements.txt`
4. Run `python3 app.py --help` for detailed instructions on launching the visualization

