// Import the processes to run in this workflow
include {
    annotate_genes;
    annotate_genes_with_abundances;
} from './processes/import_geneshot'

workflow {
 
    // If a set of geneshot results were provided
    if (params.annotate_geneshot){

        // If the metagenomic abundances were also provided
        if (params.abundances_geneshot){

            // Format the gene annotations and abundances as a CSV
            annotate_genes_with_abundances(
                Channel
                    .fromPath(
                        params.annotate_geneshot
                    ),
                Channel
                    .fromPath(
                        params.abundances_geneshot
                    )
            )

        }else{

            // Format the gene annotations (only) as a CSV
            annotate_genes(
                Channel
                    .fromPath(
                        params.annotate_geneshot
                    )
            )

        }

    }
   
}