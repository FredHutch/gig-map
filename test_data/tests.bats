#!/usr/bin/env bats

@test "Download genomes" {
    tree -lah
    TOOL=download_genomes
    rm -rf ${TOOL}
    mkdir ${TOOL}
    cd ${TOOL}
    
    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker

    # Specify which table of genomes to download and
    # run the dataset, waiting until it finishes
    wb run_dataset --genome_csv ../Escherichia_virus_T4.csv --wait

    # Print the process error and output to the screen
    cat ._wb/error.txt
    cat ._wb/output.txt
    ls -lahtr

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the genomes were downloaded
    (( $(ls genomes/*.fna.gz | wc -l) == 4 ))
}

@test "Download genes" {
    tree -lah
    TOOL=download_genes
    rm -rf ${TOOL}
    mkdir ${TOOL}
    cd ${TOOL}

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker
    
    # Specify which table of genomes to download and
    # run the dataset, waiting until it finishes
    wb run_dataset --genome_csv ../Escherichia_virus_T4.csv --wait

    # Print the process error and output to the screen
    cat ._wb/error.txt
    cat ._wb/output.txt
    ls -lahtr

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the genes were downloaded
    (( $(ls genes/*.faa.gz | wc -l) == 4 ))
}

@test "Deduplicate genes" {
    tree -lah
    TOOL=deduplicate_genes
    rm -rf ${TOOL}
    mkdir ${TOOL}
    cd ${TOOL}

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker
    
    # Specify the folder which contains the set of genes to deduplicate
    wb run_dataset --genes ../download_genes/genes/ --nxf_profile testing --wait

    # Print the process error and output to the screen
    cat ._wb/error.txt
    cat ._wb/output.txt
    cat ._wb/tool/run.sh
    cat ._wb/tool/env
    cat .nextflow.log
    ls -lahtr

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the genes were deduplicated
    [ -s centroids.faa.gz ]
    [ -s centroids.membership.csv.gz ]
    [ -s centroids.annot.csv.gz ]
}

@test "Align genomes" {
    tree -lah
    TOOL=align_genomes
    rm -rf ${TOOL}
    mkdir ${TOOL}
    cd ${TOOL}

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker

    # Specify the genes and genomes to align
    wb run_dataset \
        --genes ../deduplicate_genes/centroids.faa.gz \
        --genomes ../download_genomes/genomes \
        --pick_marker_genes 1 \
        --nxf_profile testing \
        --wait

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the outputs were created
    [ -s genomes.aln.csv.gz ]
    [ -s markers.fasta.gz ]
    [ -s gigmap.rdb ]

}

@test "Genome ANI" {
    tree -lah
    TOOL=ani
    rm -rf ${TOOL}
    mkdir ${TOOL}
    cd ${TOOL}

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker

    # Specify the genes and genomes to align
    wb run_dataset \
        --genomes ../download_genomes/genomes \
        --nxf_profile testing \
        --wait

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the outputs were created
    [ -s distances.csv.gz ]

}

@test "Collect results" {
    tree -lah

    TOOL=collect
    rm -rf ${TOOL}
    mkdir ${TOOL}
    cd ${TOOL}

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker

    # Specify the genes and genomes to align
    wb run_dataset \
        --genomes ../download_genomes/genomes \
        --genome_aln ../align_genomes/genomes.aln.csv.gz \
        --marker_genes ../align_genomes/markers.fasta.gz \
        --nxf_profile testing \
        --wait

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the outputs were created
    [ -s gigmap.rdb ]
    (( $(find ./ -name "gigmap*.feather" | wc -l) > 0 ))

}

@test "Render HTML" {
    tree -lah

    TOOL=render
    rm -rf ${TOOL}
    mkdir ${TOOL}
    cd ${TOOL}

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker

    # Specify the genes and genomes to align
    wb run_dataset \
        --genome_aln ../align_genomes/genomes.aln.csv.gz \
        --genome_annot ../align_genomes/genome.manifest.csv \
        --gene_annot ../align_genomes/gene.manifest.csv \
        --genome_distmat ../align_genomes/distances.csv.gz \
        --gene_order ../align_genomes/genomes.gene.order.txt.gz \
        --genomeAnnot_label_col "Formatted Name" \
        --geneAnnot_label_col "Formatted Name" \
        --nxf_profile testing \
        --wait

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the outputs were created
    (( $(find ./ -name "gigmap*.html" | wc -l) > 0 ))
    (( $(find ./ -name "gigmap*.feather" | wc -l) > 0 ))

}

@test "Align Reads" {
    tree -lah

    TOOL=align_reads
    rm -rf ${TOOL}_single
    mkdir ${TOOL}_single
    cd ${TOOL}_single

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker

    # Specify the genes and reads to align
    wb run_dataset \
        --genes ../GCA_000005845.2_ASM584v2_protein.faa.gz \
        --reads ../fastq \
        --nxf_profile testing \
        --wait

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the outputs were created
    (( $(find ./ -name "read_alignments.csv.gz" | wc -l) > 0 ))
    [ -s alignment_logs.csv ]
    [ -s famli_logs.csv ]

}

@test "Align Paired Reads" {
    tree -lah

    TOOL=align_reads
    rm -rf ${TOOL}_paired
    mkdir ${TOOL}_paired
    cd ${TOOL}_paired

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker

    # Specify the genes and reads to align
    wb run_dataset \
        --genes ../GCA_000005845.2_ASM584v2_protein.faa.gz \
        --reads ../fastq \
        --paired 1 \
        --nxf_profile testing \
        --wait

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the outputs were created
    (( $(find ./ -name "read_alignments.csv.gz" | wc -l) > 0 ))

}

@test "Sketch genomes" {
    tree -lah

    TOOL=sketch_genomes
    rm -rf ${TOOL}
    mkdir ${TOOL}
    cd ${TOOL}

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker

    # Specify the genes and reads to align
    wb run_dataset \
        --genomes ../download_genomes/genomes/ \
        --sketch_folder genome_sketches/ \
        --task_limit 2 \
        --nxf_profile testing \
        --wait

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the outputs were created
    (( $(find genome_sketches/ -name "combined_genomes.msh" | wc -l) == 1 ))

}

@test "Search sketches" {
    tree -lah

    TOOL=search_sketches
    rm -rf ${TOOL}
    mkdir ${TOOL}
    cd ${TOOL}

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker

    # Specify the genes and reads to align
    wb run_dataset \
        --query '../download_genes/genes/*.faa.gz' \
        --genome_sketches ../sketch_genomes/genome_sketches/combined_genomes.msh \
        --search_results search_results \
        --task_limit 2 \
        --nxf_profile testing \
        --wait

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the outputs were created
    (( $(find search_results -name "*.csv" | wc -l) > 0 ))

}

@test "Test Reads" {
    tree -lah

    TOOL=test_reads
    rm -rf ${TOOL}
    mkdir ${TOOL}
    cd ${TOOL}

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker

    # Specify the gene abundances, manifest, and columns to test
    wb run_dataset \
        --gene_abund ../align_reads_single/read_alignments.csv.gz \
        --manifest ../read_manifest.csv \
        --formula GROUP \
        --nxf_profile testing \
        --wait

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the outputs were created
    (( $(find ./ -name "mu.GROUP.results.csv.gz" | wc -l) > 0 ))

}

@test "Map Genes - blast" {
    tree -lah

    TOOL=map_genes
    rm -rf ${TOOL}-blast
    mkdir ${TOOL}-blast
    cd ${TOOL}-blast

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker

    # Specify the gene abundances, manifest, and columns to test
    wb run_dataset \
        --references ../deduplicate_genes/centroids.faa.gz \
        --queries ../download_genes/genes/GCA_002599625.1_ASM259962v1_protein.faa.gz \
        --map_batchsize 50 \
        --aligner blast \
        --nxf_profile testing \
        --wait

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the outputs were created
    [ -s gene_mapping.csv.gz ]

}

@test "Map Genes - diamond" {
    tree -lah

    TOOL=map_genes
    rm -rf ${TOOL}-diamond
    mkdir ${TOOL}-diamond
    cd ${TOOL}-diamond

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker

    # Specify the gene abundances, manifest, and columns to test
    wb run_dataset \
        --references ../deduplicate_genes/centroids.faa.gz \
        --queries ../download_genes/genes/GCA_002599625.1_ASM259962v1_protein.faa.gz \
        --map_batchsize 50 \
        --aligner diamond \
        --nxf_profile testing \
        --wait

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the outputs were created
    [ -s gene_mapping.csv.gz ]

}

@test "Pairwise Similarity" {
    tree -lah

    TOOL=pdist
    rm -rf ${TOOL}-diamond
    mkdir ${TOOL}-diamond
    cd ${TOOL}-diamond

    # Specify the tool and launcher to use
    wb setup_dataset --tool ${TOOL} --launcher gigmap_docker

    # Specify the gene abundances, manifest, and columns to test
    wb run_dataset \
        --genes ../deduplicate_genes/centroids.faa.gz \
        --map_batchsize 50 \
        --min_coverage 10 \
        --aligner diamond \
        --nxf_profile testing \
        --wait

    # Print the logs
    cat ._wb/output.txt
    cat ._wb/error.txt

    # Make sure that the outputs were created
    [ -s gene_pdist.json.gz ]

}
