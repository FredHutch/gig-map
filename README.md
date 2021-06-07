# GiG-map (Genes in Genomes - Map)
Build a map of genes present across a set of microbial genomes

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
