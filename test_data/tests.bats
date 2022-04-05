#!/usr/bin/env bats

@test "Download genomes" {
    rm -rf download_genomes
    mkdir download_genomes
    cd download_genomes
    
    # Specify the tool and launcher to use
    wb setup_dataset \
    --tool download_genomes \
    --launcher nextflow_docker

    # Specify which table of genomes to download and
    # run the dataset, waiting until it finishes
    wb run_dataset --genome_csv ../Escherichia_virus_T4.csv --wait

    # Make sure that the genomes were downloaded
    (( $(ls genomes/*.fna.gz | wc -l) == 8 ))
}

@test "Download genes" {
    rm -rf download_genes
    mkdir download_genes
    cd download_genes

    # Specify the tool and launcher to use
    wb setup_dataset \
    --tool download_genes \
    --launcher nextflow_docker
    
    # Specify which table of genomes to download and
    # run the dataset, waiting until it finishes
    wb run_dataset --genome_csv ../Escherichia_virus_T4.csv --wait

    # Make sure that the genes were downloaded
    (( $(ls genes/*.faa.gz | wc -l) == 8 ))
}
