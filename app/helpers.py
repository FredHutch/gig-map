#!/usr/bin/env python3

from copy import copy
import logging
from numpy import select
import pandas as pd
import plotly.graph_objects as go
from skbio import DistanceMatrix
from skbio.tree import nj

def remove_genome_file_ext(fp):
    for ext in ['.gz', '.fna', '.fasta', '.fa']:
        if fp.endswith(ext):
            fp = fp[:-len(ext)]
    return fp

def read_data(args):
    """Read in the data needed to render the heatmap."""

    # Get the logger
    logger = logging.getLogger('gig-map')

    # Read in the alignments
    logger.info(f"Reading from {args['alignments']}")
    alignments = pd.read_csv(args['alignments'])

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

    # By default, there are no additional labels to apply to the genes
    available_gene_labels = []

    # Read in the pairwise genome distances
    logger.info(f"Reading from {args['distances']}")
    dists = pd.read_csv(
        args['distances'],
        sep="\t",
        index_col=0
    )

    # Read in the gene annotations, if any
    if args['gene_annotations'] is None:
        gene_annotations = None

    # If a table was provided
    else:
        # Read in the table
        logger.info(f"Reading from {args['gene_annotations']}")
        gene_annotations = pd.read_csv(args['gene_annotations'])

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

            # Also add those columns to the options available for labeling genes
            available_gene_labels.append(
                dict(
                    label=col_name,
                    value=col_name,
                )
            )

    # By default, there are no additional labels to apply to the genomes
    available_genome_labels = []

    # Read in the genome annotations, if any
    if args['genome_annotations'] is None:
        genome_annotations = None

    # If a table was provided
    else:
        # Read in the table
        logger.info(f"Reading from {args['genome_annotations']}")
        genome_annotations = pd.read_csv(args['genome_annotations'])

        # Make sure that there is a column named "genome_id"
        msg = "Genome annotation CSV must contain a column named 'genome_id'"
        msg = f"{msg} - found {'; '.join(genome_annotations.columns.values)}"
        assert "genome_id" in genome_annotations.columns.values, msg

        # Set the index of the table as 'genome_id'
        genome_annotations.set_index('genome_id', inplace=True)

        # Add the other columns to the available genome annotations
        for col_name in genome_annotations.columns.values:

            # Add those columns to the options available for labeling genomes
            available_genome_labels.append(
                dict(
                    label=col_name,
                    value=col_name,
                )
            )

    logger.info("Done reading all data")

    # Return data formatted as a dict
    return dict(
        alignments=alignments,
        dists=dists,
        gene_annotations=gene_annotations,
        available_gene_annotations=available_gene_annotations,
        available_gene_labels=available_gene_labels,
        genome_annotations=genome_annotations,
        available_genome_labels=available_genome_labels,
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

    # Assign x/y to create a DataFrame
    node_positions = CartesianTree(
        tree,
        y_offset=-0.5,
    )

    # Return the layout of the tree
    return node_positions

def plot_heatmap(plot_df, node_positions, selections, data, xaxis='x', yaxis='y'):
    """Return a heatmap rendered from the gene DataFrame and the genome positions."""

    plot_df = plot_df.rename(
        index=lambda s: s.replace("_", " ")
    ).reindex(
        index=node_positions.genome_order
    )

    # If the user elected to label the genes by something other than their ID
    if selections["label-genes-by"] != "":

        # Rename the columns of the DataFrame
        plot_df = plot_df.rename(
            columns=data["gene_annotations"][selections["label-genes-by"]].get
        )

    return go.Heatmap(
        x=list(plot_df.columns.values),
        z=plot_df.values,
        xaxis=xaxis,
        yaxis=yaxis,
        colorscale=selections["heatmap-colorscale"],
        showscale=selections["show-on-right"] == "colorscale",
        # Title of the colorbar
        colorbar=dict(
            title=dict(
                pident="Alignment Identity",
                coverage="Alignment Coverage",
            ).get(
                selections["color-genes-by"],
                selections["color-genes-by"]
            ),
        ),
    )


def plot_tree(node_positions, selections, data, xaxis='x', yaxis='y'):
    """Return a Plotly trace rendered from a skbio tree."""

    # If the user decided to label the genomes
    if selections["label-genomes-by"] != "":

        # Replace the values in the 'name' column of the DataFrame used for plotting
        node_positions.df = node_positions.df.replace(
            to_replace=dict(
                name=data[
                    "genome_annotations"
                ][
                    selections["label-genomes-by"]
                ].to_dict()
            )
        )

    # Return a ScatterGL
    return go.Scattergl(
        name="Neighbor Joining Tree",
        showlegend=False,
        mode="lines",
        x=node_positions.x_coords(),
        y=node_positions.y_coords(),
        text=node_positions.text(),
        hoverinfo="text",
        xaxis=xaxis,
        yaxis=yaxis,
    )

class CartesianTree:

    def __init__(self, tree, y_offset=0, x=0):

        # Set up each position as a dict in a list
        # keys will be name, x, y, and parent
        self.positions = []

        # Assign the root, and then recurse down to the tips
        self.add_clade(tree, y_offset=y_offset, x=x)

        # Set up a DataFrame with the coordinates
        self.df = pd.DataFrame(self.positions)

        # Save the list of genome positions from the tree layout
        self.genome_order = self.df.query("is_leaf").set_index(
            "name"
        )[
            "y"
        ].sort_values().index.values

    def x_coords(self):
        """Return a list of x-coordinates to use for plotting the tree."""

        # Format is a list with each node and its parent, separated by NaN values
        return self._list_link_to_parents(col_name="x")

    def y_coords(self):
        """Return a list of y-coordinates to use for plotting the tree."""

        # Format is a list with each node and its parent, separated by NaN values
        return self._list_link_to_parents(col_name="y")

    def text(self):
        """Return a list of leaf labels to use for plotting the tree."""

        # Format is a list with each node and its parent, separated by NaN values
        return self._list_link_to_parents(col_name="name")

    def _list_link_to_parents(self, col_name="x"):
        """Internal method for generating a list of values from self.df"""

        # col_name may only be x, y, or name
        assert col_name in ['x', 'y', 'name'], f"Not recognized: {col_name}"

        # Populate a list which will be output
        output_list = []

        # Iterate over each row
        for _, r in self.df.iterrows():

            # If there is no parent for this row
            if pd.isnull(r['parent']) or r['parent'] < 0:

                # Skip it
                continue

            # If there is a parent for this row
            else:

                # Add the item to the list
                output_list.append(r[col_name])

                # The value placed in between the node and its parent
                # depends on whether it is X or Y (or name)

                # Get the value for the parent
                parent_val = self.df.loc[r['parent'], col_name]

                # Moving along the x axis
                if col_name == "x":

                    # The intermediate node has the X coordinate of the parent
                    output_list.append(parent_val)

                # Moving along the y axis
                elif col_name == "y":

                    # The intermediate node has the Y coordinate of the child
                    output_list.append(r[col_name])

                # For the 'name' values
                elif col_name == "name":

                    # The intermediate node has no name
                    output_list.append(None)

                # Add its parent
                output_list.append(parent_val)

                # Add a NaN to separate it
                output_list.append(None)

        # Return the list
        return output_list

    def add_clade(self, clade, y_offset=0, x=0, parent=-1):
        """Add the node at the base of a clade, then add its children (if any)."""

        # Set up a numeric ID for this clade
        clade_id = len(self.positions)

        # Calculate the number of tips for this clade
        clade_n_tips = len(list(clade.tips(include_self=True)))

        # The Y position is based on the total number of tips, and the offset
        clade_y = y_offset + (clade_n_tips / 2.)

        # The position of this clade is in the middle of max_y and min_y
        self.positions.append(
            dict(
                name=clade.name,
                x=x,
                y=clade_y,
                parent=parent,
                is_leaf=clade_n_tips==1
            )
        )

        # For each child
        for child in clade.children:

            # Calculate the number of tips for this child
            child_n_tips = len(list(child.tips(include_self=True)))

            # Calculate the x position of the clade by adding
            # the distance to its parent
            child_x = x + clade.distance(child)

            # Add the child
            self.add_clade(
                child,
                y_offset=y_offset,
                x=child_x,
                parent=clade_id,
            )

            # The next child will be positioned above, based on the number of tips
            y_offset = y_offset + child_n_tips
