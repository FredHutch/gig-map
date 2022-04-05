#!/usr/bin/env bats

@test "Download genomes" {
    rm -rf download_genomes
    mkdir download_genomes
    cd download_genomes
    
    # Specify the tool and launcher to use
    wb setup_dataset \
    --tool gig-map/download_genomes \
    --launcher FredHutch_bash-workbench-tools/nextflow_docker
    
    # Specify which table of genomes to download
    wb set_tool_params --genome_csv ../Escherichia_virus_T4.csv

    # Set the params for the nextflow docker launcher (there are none)
    wb set_launcher_params

    # Run the dataset and wait until it finishes
    wb run_dataset --wait

    # Make sure that the genomes were downloaded
    (( $(ls genomes/*.fna.gz | wc -l) == 8 ))
}

@test "Download genes" {
    rm -rf download_genes
    mkdir download_genes
    cd download_genes

    # Specify the tool and launcher to use
    wb setup_dataset \
    --tool gig-map/download_genes \
    --launcher FredHutch_bash-workbench-tools/nextflow_docker
    
    # Specify which table of genes to download
    wb set_tool_params --genome_csv ../Escherichia_virus_T4.csv

    # Set the params for the nextflow docker launcher (there are none)
    wb set_launcher_params

    # Run the dataset and wait until it finishes
    wb run_dataset --wait

    # Make sure that the genes were downloaded
    (( $(ls genes/*.faa.gz | wc -l) == 8 ))
}
