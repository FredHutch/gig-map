#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.container__gigmap = "quay.io/hdc-workflows/gig-map:73f25a1"
params.mem_gbs = 4
params.rdb = false
params.output_folder = "output"
params.output_prefix = "gig-map"
params.options = ""
params.gene_annotations = false
params.genome_annotations = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/gig-map/render.nf <ARGUMENTS>

    Required Arguments:
      --rdb                 File contining output of the gig-map pipeline, file with the .rdb extension
      --output_prefix       Name of the HTML output (default: ${params.output_prefix})
      --output_folder       Folder for output file (default: ${params.output_folder})
      --gene_annotations    Optional CSV containing gene annotations
      --genome_annotations  Optional CSV containing genome annotations
      --mem_gbs             Amount of memory (in GBs) to use for rendering plots (default: ${params.mem_gbs})
      --container__gigmap   Docker container used for rendering plots (default: ${params.container__gigmap})

      --options             Any additional options which should be provided when rendering the static gig-map output
                            The format of this parameter will seem very odd, because it is actually a string
                            which will be added in its entirety to the gig-map-cli invokation.

                            The complete list of options available can be seen in the following file:
                            https://github.com/FredHutch/gig-map/blob/main/app/gig-map-cli

                            Note that the --rdb, --gene-annotations, and --genome-annotations options can be
                            omitted, since they are taken care of by the larger workflow.

                            As an example, to run the script with a more stringent threshold for alignment
                            identity and coverage, you would set the --options flag with the following string,
                            making sure to enclose the entire string in quotes as shown here:

                            "--min-pctid 95 --min-cov 95"

                            Additional options, in combination with the gene and genome annotation file options
                            can be extremely helpful in making an informative display.

    """.stripIndent()
}

// Define the process for rendering an HTML file from an RDB input
process render {
    container "${params.container__gigmap}"
    memory "${params.mem_gbs}.GB"
    publishDir "${params.output_folder}"

    input:
    path RDB

    output:
    path "${params.output_prefix}.html"

    script:
    """#!/bin/bash

set -e

gig-map-cli \
    --rdb "${RDB}" \
    --output-prefix ${params.output_prefix} \
    --output-folder ./
    """
}

process render_genes_genomes {
    container "${params.container__gigmap}"
    memory "${params.mem_gbs}.GB"
    publishDir "${params.output_folder}"

    input:
    path RDB
    path GENE_ANNOTATIONS
    path GENOME_ANNOTATIONS

    output:
    path "${params.output_prefix}.html"

    script:
    """#!/bin/bash

set -e

gig-map-cli \
    --rdb "${RDB}" \
    --gene-annotations "${GENE_ANNOTATIONS}" \
    --genome-annotations "${GENOME_ANNOTATIONS}" \
    --output-prefix ${params.output_prefix} \
    --output-folder ./
    """
}

process render_genes {
    container "${params.container__gigmap}"
    memory "${params.mem_gbs}.GB"
    publishDir "${params.output_folder}"

    input:
    path RDB
    path GENE_ANNOTATIONS
    path GENOME_ANNOTATIONS

    output:
    path "${params.output_prefix}.html"

    script:
    """#!/bin/bash

set -e

gig-map-cli \
    --rdb "${RDB}" \
    --gene-annotations "${GENE_ANNOTATIONS}" \
    --output-prefix ${params.output_prefix} \
    --output-folder ./
    """
}

process render_genomes {
    container "${params.container__gigmap}"
    memory "${params.mem_gbs}.GB"
    publishDir "${params.output_folder}"

    input:
    path RDB
    path GENOME_ANNOTATIONS

    output:
    path "${params.output_prefix}.html"

    script:
    """#!/bin/bash

set -e

gig-map-cli \
    --rdb "${RDB}" \
    --genome-annotations "${GENOME_ANNOTATIONS}" \
    --output-prefix ${params.output_prefix} \
    --output-folder ./
    """
}

workflow {

    // Show help message if the user specifies the --help flag at runtime
    if (params.help){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }

    // If an RDB file is not provided
    if (!params.rdb){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }

    // If both a gene and genome annotation file is provided
    if (params.gene_annotations && params.genome_annotations){
        render_genes_genomes(
            file(params.rdb),
            file(params.gene_annotations),
            file(params.genome_annotations)
        )
    } else {
        if (params.gene_annotations){
            render_genes(
                file(params.rdb),
                file(params.gene_annotations)
            )
        }
        if (params.genome_annotations){
            render_genomes(
                file(params.rdb),
                file(params.genome_annotations)
            )
        }
        if (!params.gene_annotations && !params.genome_annotations){
            render(
                file(params.rdb)
            )
        }
    }
}