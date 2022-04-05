#!/usr/bin/env bats

@test "Download genomes" {
    TOOL=download_genomes
    rm -rf ${TOOL}
    mkdir ${TOOL}
    cd ${TOOL}
    
    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher nextflow_docker

    # Specify which table of genomes to download and
    # run the dataset, waiting until it finishes
    wb run_dataset --genome_csv ../Escherichia_virus_T4.csv --wait

    # Print the process error and output to the screen
    cat ._wb/error.txt
    cat ._wb/output.txt
    ls -lahtr

    # Make sure that the genomes were downloaded
    (( $(ls genomes/*.fna.gz | wc -l) == 8 ))
}

@test "Download genes" {
    TOOL=download_genes
    rm -rf ${TOOL}
    mkdir ${TOOL}
    cd ${TOOL}

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher nextflow_docker
    
    # Specify which table of genomes to download and
    # run the dataset, waiting until it finishes
    wb run_dataset --genome_csv ../Escherichia_virus_T4.csv --wait

    # Print the process error and output to the screen
    cat ._wb/error.txt
    cat ._wb/output.txt
    ls -lahtr

    # Make sure that the genes were downloaded
    (( $(ls genes/*.faa.gz | wc -l) == 8 ))
}

@test "Deduplicate genes" {
    TOOL=deduplicate_genes
    rm -rf ${TOOL}
    mkdir ${TOOL}
    cd ${TOOL}

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher nextflow_docker
    
    # Specify the folder which contains the set of genes to deduplicate
    wb run_dataset --genes ../download_genes/genes/ --wait

    # Print the process error and output to the screen
    cat ._wb/error.txt
    cat ._wb/output.txt
    ls -lahtr

    # Make sure that the genes were deduplicated
    [ -s centroids.faa.gz ]
    [ -s centroids.membership.csv.gz ]
    [ -s centroids.annot.csv.gz ]
}
