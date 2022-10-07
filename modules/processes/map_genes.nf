// Align two gene collections against each other with BLAST
process shard {
    container "${params.container__pandas}"
    label 'io_limited'
    
    input:
    path "queries.fasta.gz"

    output:
    path "queries.shard.*.fasta.gz"

    script:
    template "map_genes_shard.py"
}

// Align two gene collections against each other with BLAST
process map_genes_blast {
    container "${params.container__blast}"
    label 'mem_veryhigh'
    
    input:
    path "queries.fasta.gz"
    path "refs.fasta.gz"

    output:
    path "unfiltered_gene_mapping.csv.gz"

    script:
    template "map_genes_blast.sh"
}

// Align two gene collections against each other with DIAMOND
process map_genes_diamond {
    container "${params.container__diamond}"
    label 'mem_veryhigh'
    
    input:
    path "queries.fasta.gz"
    path "refs.fasta.gz"

    output:
    path "unfiltered_gene_mapping.csv.gz"

    script:
    template "map_genes_diamond.sh"
}

// Only keep the top hit per gene
process filter {
    container "${params.container__pandas}"
    label 'io_limited'
   
    input:
    path "unfiltered_gene_mapping.csv.gz"

    output:
    path "gene_mapping.csv.gz"

    script:
    template "map_genes_filter.py"
}

// Only keep the top hit per gene
process join {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}", mode: 'copy', overwrite: true
   
    input:
    path "gene_mapping.shard.*.csv.gz"

    output:
    path "gene_mapping.csv.gz"

    script:
    template "map_genes_join.py"
}
