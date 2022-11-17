#!/usr/bin/env python3.8

import click
import gzip
import json
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial import distance

@click.command()
@click.option("--distmat")
@click.option("--alignments")
@click.option("--genome_annot")
@click.option("--gene_annot")
@click.option("--gene_order")
@click.option("--output_folder")
@click.option("--output_prefix")
@click.option("--options_json")
def serialize(
    distmat,
    alignments,
    genome_annot,
    gene_annot,
    gene_order,
    output_folder,
    output_prefix,
    options_json
):

    data = read_data(
        distmat,
        alignments,
        genome_annot,
        gene_annot,
        gene_order,
        options_json
    )

    # Order the genomes based on linkage clustering
    print("Ordering genomes")
    genome_order = order_genomes(data)
    
    # Make a wide table with the genome alignments
    wide_df = data["alignments"].pivot_table(
        index="genome",
        columns="sseqid",
        values="pident",
        aggfunc=max
    ).reindex(
        index=genome_order,
        columns=data["gene_order"]
    ).fillna(
        0
    )
    print(f"Constructed a gene-genome table with {wide_df.shape[0]:,} genomes and {wide_df.shape[1]:,} genes")

    # If there are genome annotations
    genomeAnnot_label_col = data["options"].get("genomeAnnot_label_col")
    if genomeAnnot_label_col is not None and genomeAnnot_label_col in data["genome_annot"].columns.values:
        print(f'Renaming genomes by {genomeAnnot_label_col}')
        wide_df = wide_df.rename(
            index=data["genome_annot"][genomeAnnot_label_col].get
        )

        data['distmat'] = data['distmat'].rename(
            index=data["genome_annot"][genomeAnnot_label_col].get,
            columns=data["genome_annot"][genomeAnnot_label_col].get
        ).reindex(
            index=wide_df.index.values,
            columns=wide_df.index.values,
        )

    # If there are gene annotations
    geneAnnot_label_col = data["options"].get("geneAnnot_label_col")
    if geneAnnot_label_col is not None and geneAnnot_label_col in data["gene_annot"].columns.values:
        print(f'Renaming genes by {geneAnnot_label_col}')
        wide_df = wide_df.rename(
            columns=data["gene_annot"][geneAnnot_label_col].get
        )

    # Write out two tables, the distance matrix and the gene table
    # Only columns will be labeled
    # The order of columns on the distance matrix will match the order of rows on the gene table
    dm_fp = f"{output_folder.rstrip('/')}/{output_prefix}.dm.feather"
    write_feather(data['distmat'], dm_fp)

    genes_fp = f"{output_folder.rstrip('/')}/{output_prefix}.genes.feather"
    write_feather(wide_df, genes_fp)

    print("Done")


def write_feather(df, fp):
    if not fp.endswith(".feather"):
        fp = fp + ".feather"

    # Get the column names
    cnames = df.columns.values
    
    # Write out the column names
    cnames_fp = fp.replace(".feather", ".columns.txt.gz")
    print(f"Writing to {cnames_fp}")
    with gzip.open(cnames_fp, "wt") as handle:
        handle.write("\n".join(cnames))

    # Override the DataFrame columns
    df.columns = list(map(str, range(df.shape[1])))

    # Write the feather file
    print(f"Writing to {fp}")
    df.reset_index(drop=True).to_feather(fp)


def order_genomes(data, method="average"):

    return data['distmat'].index.values[
        hierarchy.leaves_list(
            hierarchy.linkage(
                distance.squareform(data['distmat'].values),
                method=method
            )
        )
    ]


def read_data(
        distmat,
        alignments,
        genome_annot,
        gene_annot,
        gene_order,
        options_json
    ) -> dict:

    data = dict()

    # Read in the options file
    data["options"] = json.load(open(options_json))

    # Read in the distance matrix
    data["distmat"] = read_csv(distmat, index_col=0)

    # Read in the alignments
    data["alignments"] = read_csv(alignments, usecols=["genome", "sseqid", "pident"])
    
    # Read in the genome annotations
    data["genome_annot"] = read_csv(genome_annot).set_index("genome_id")

    # Read in the gene annotations
    data["gene_annot"] = read_csv(gene_annot).set_index("gene_id")

    # Read in the gene order
    data["gene_order"] = read_csv(gene_order, header=None)[0].to_list()

    return data


def read_csv(fp, **kwargs) -> pd.DataFrame:
    print(f"Reading in {fp}")
    df = pd.read_csv(fp, **kwargs)
    print(f"Read in {df.shape[0]:,} rows and {df.shape[1]:,} columns")
    return df


if __name__ == "__main__":
    serialize()
