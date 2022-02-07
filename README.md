# GiG-map (Genes in Genomes - Map)
Build a map of genes present across a set of microbial genomes

[![Docker Repository on Quay](https://quay.io/repository/hdc-workflows/gig-map/status "Docker Repository on Quay")](https://quay.io/repository/hdc-workflows/gig-map)

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

# Documentation

For a more complete set of documentation, please [visit the Wiki](https://github.com/FredHutch/gig-map/wiki)

# Implementation

The execution of the analysis encoded by this repository is roughly grouped into two
activities: 1) alignment of genes against genomes and 2) visualization of those alignment
results. 

## Aligning Genes To Genomes

The analysis workflow used to align genes to genomes (and metagenomes) is organized around
the concept of a "project folder" which contains the inputs and outputs for a particular
analysis. To provide the user with a good amount of flexibility, the inputs for an analysis
can be sourced from any location, but the _default_ location of the inputs will be from
the project folder.

### Quick Reference

```
gigmap-project/
│
│   # INPUT FILES
│
├── genome_tables/
│   └── *.csv              # All genomes from *.csv files in this folder will be downloaded;
│                          # Source: https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/
├── genomes/
│   └── *.gz               # All files in this folder must be FASTA-formatted (gzip-compressed)
│                          # genomes to include in this analysis;
├── gene_tables/
│   └── *.csv              # All genes from *.csv files in this folder will be downloaded;
│                          # Source: https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/
├── genes/
│   └── *.gz               # All files in this folder must be FASTA-formatted (gzip-compressed)
│                          # genes to include in this analysis;
├── markers/
│   └── *.gz               # (optional) Any genes included in this folder will be aligned across
│                          # all genomes and used to estimate phylogenies;
│
├── metagenomes/
│   └── single_end/*.fastq.gz
│                          # (optional) Short reads in this folder will be aligned to the gene catalog
│
│   └── paired_end/*_R{1,2}*.fastq.gz
│                          # (optional) Short reads in this folder will be aligned to the gene catalog
│
│   # OUTPUT FILES
│
├── downloaded_genomes/    # (module: download)
│   ├── genomes.annot.csv.gz
│   │                      # Annotations provided for the downloaded genomes
│   └── *_genomic.fna.gz   # Genomes downloaded from genome_tables/
│
├── downloaded_genes/      # (module: download)
│   └── *_protein.faa.gz   # Genes downloaded from gene_tables/
│
├── deduplicated_genes/    # (module: deduplicate)
│   ├── centroids.faa.gz   # Non-redundant set of genes 
│   ├── centroids.annot.csv.gz 
│   │                      # Annotations applied to deduplicated genes
│   └── centroids.membership.csv.gz
│                          # Table indicating which centroids represent which input genes
│
├── genome_alignment/      # (module: align_genomes)
│   ├── genomes.aln.csv.gz # Tabular results of the alignment of genes to genomes
│   │                      
│   ├── genomes.gene.order.txt.gz
│   │                      # Ordination of genes based on genome membership
│   └── genes/
│       └── *.gz           # Nucleotide sequences of each gene from all genomes
│ 
├── ani/                   # (module: ani)
│   ├── distances.csv.gz   # Square distance matrix with pairwise ANI for all genomes
│   │
│   ├── msh/
│   │   └── *.msh          # Mash sketch generated from all genomes
│   │
│   └── tsv/
│       └── *.tsv.gz       # Long-format table with distances for each genome
│
├── metagenomes/           # (module: align_reads)
│   ├── alignments/
│   │   └── *.json.gz      # Alignment information for each individual metagenome
│   │
│   └── alignments.csv.gz  # Summary of all alignment information per gene across samples
│
└── package                # (module: aggregate)
    └── gigmap.rdb         # Collection of gig-map output data in RDB format for
                           # interactive visualization with app/gig-map
```

## GiG-map Visualization

To visualize the results of the gig-map alignment workflow, the output files may be
rendered using the interactive tool provided as `app.py`. In order to use this tool:

1. Clone this repository locally
2. Create a python virtual environment (tested with Python 3.8 and higher)
3. Install the python dependencies contained in `requirements.txt`
4. Run `python3 app.py --help` for detailed instructions on launching the visualization

