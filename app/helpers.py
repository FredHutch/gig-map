#!/usr/bin/env python3

import pandas as pd
from skbio import DistanceMatrix
from skbio.tree import nj

def remove_genome_file_ext(fp):
    for ext in ['.gz', '.fna', '.fasta', '.fa']:
        if fp.endswith(ext):
            fp = fp[:-len(ext)]
    return fp

def read_data(alignments_csv, dists_tsv, gene_annotations):
    """Read in the data needed to render the heatmap."""

    # Read in the alignments
    alignments = pd.read_csv(alignments_csv)

    # Calculate the alignment coverage of each gene
    alignments = alignments.assign(
        coverage = alignments.apply(
            lambda r: 100 * (r['send'] - r['sstart'] + 1) / r['slen'],
            axis=1
        )
    )

    # Remove the file endings from the genome names
    alignments = alignments.apply(
        lambda c: c.apply(remove_genome_file_ext) if c.name == "genome" else c
    )

    # The default annotations available for each gene are the
    # alignment identity and coverage
    available_gene_annotations = [
        dict(
            label="Alignment Identity",
            value="pident"
        ),
        dict(
            label="Alignment Coverage",
            value="coverage"
        )
    ]

    # Read in the pairwise genome distances
    dists = pd.read_csv(
        dists_tsv,
        sep="\t",
        index_col=0
    )

    # Read in the gene annotations, if any
    if gene_annotations is None:
        gene_annotations = None

    # If a table was provided
    else:
        # Read in the table
        gene_annotations = pd.read_csv(gene_annotations)

        # Make sure that there is a column named "gene_id"
        msg = "Gene annotation CSV must contain a column named 'gene_id'"
        msg = f"{msg} - found {'; '.join(gene_annotations.columns.values)}"
        assert "gene_id" in gene_annotations.columns.values, msg

        # Set the index of the table as 'gene_id'
        gene_annotations.set_index('gene_id', inplace=True)

        # Add the other columns to the available gene annotations
        for col_name in gene_annotations.columns.values:
            available_gene_annotations.append(
                dict(
                    label=col_name,
                    value=col_name,
                )
            )

    # Return data formatted as a dict
    return dict(
        alignments=alignments,
        dists=dists,
        gene_annotations=gene_annotations,
        available_gene_annotations=available_gene_annotations,
    )


# Generate a neighbor-joining tree from a subset of genomes
def make_nj_tree(genome_list, dists_df):

    # Make sure that we have distances for every genome
    for genome_id in genome_list:
        msg = f"Could not find a distance for genome {genome_id}"
        msg = f"{msg} -- best performance for FASTA files ending with .fasta[.gz] or .fna[.gz]"
        assert genome_id in dists_df.index.values, msg

    # Subset the distance matrix to the genomes in the list
    subset_dists = dists_df.reindex(
        columns=genome_list,
        index=genome_list,
    )

    # Format the distances as expected by skbio
    dm = DistanceMatrix(subset_dists.values, genome_list)

    # Make a neighbor-joining tree
    tree = nj(dm)

    # Root at midpoint
    tree = tree.root_at_midpoint()

    # Return the tree
    return tree


def plot_tree(tree):
    """Return a Plotly trace rendered from a skbio tree."""
    
    # The first challenge is to assign an X-Y coordinate for each node in the tree

    # The total height of the tree is going to be the number of tips
    min_y = 0
    max_y = len(list(tree.tips()))

    # Assign x/y to create a DataFrame
    node_pos_df = CartesianTree(
        tree,
        max_y=max_y,
        min_y=min_y,
    ).df()

    print(node_pos_df)

class CartesianTree:

    def __init__(self, tree, max_y=1, min_y=0, x=0):

        # Set up each position as a dict in a list
        # keys will be name, x, y, and parent
        self.positions = []

        # Assign the root, and then recurse down to the tips
        self.add_clade(tree, max_y=max_y, min_y=min_y, x=0)

    def df(self):

        # Return a DataFrame
        return pd.DataFrame(self.positions)

    def add_clade(self, clade, max_y=1, min_y=0, x=0, parent=-1):
        """Add the node at the base of a clade, then add its children (if any)."""

        # Set up a numeric ID for this clade
        clade_id = len(self.positions)

        # The position of this clade is in the middle of max_y and min_y
        self.positions.append(
            dict(
                name=clade.name,
                x=x,
                y=(max_y - min_y) / 2.,
                parent=parent
            )
        )

        # Calculate the total y span
        total_y_span = max_y - min_y

        # Calculate the range provided for each tip
        tip_y_span = total_y_span / len(list(clade.tips(include_self=True)))

        # For each child
        for child in clade.children:

            # Set the minimum y value for the clade based on the number of tips
            child_min_y = max_y - (tip_y_span * len(list(child.tips(include_self=True))))

            # Calculate the x position of the clade by adding
            # the distance to its parent
            child_x = x + clade.distance(child)

            # Add the child
            self.add_clade(
                child,
                max_y = max_y,
                min_y = child_min_y,
                x=child_x,
                parent=clade_id,
            )

            # The next child will be positioned below
            max_y = child_min_y

