#!/usr/bin/env python3

# Import the helper functions

from app.helpers import plot_colorbar, read_data, make_nj_tree, plot_tree, plot_heatmap, plot_colorbar

# Import the menu-driven-figure library
from menu_driven_figure.app import MenuDrivenFigure

import argparse
from direct_redis import DirectRedis
from functools import lru_cache
import logging
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

##################
# SET UP LOGGING #
##################

# Set the level of the logger to INFO
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [gig-map] %(message)s'
)
logger = logging.getLogger('gig-map')
logger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)


###################
# PARSE ARGUMENTS #
###################

# Create the parser
parser = argparse.ArgumentParser(
    description='Display the distribution of annotated genes across microbial genomes'
)

# Add the arguments
parser.add_argument(
    '--host',
    type=str,
    default="localhost",
    help='Redis host used for reading alignment data'
)
parser.add_argument(
    '--port',
    type=int,
    default=6379,
    help='Redis port used for reading alignment data'
)
parser.add_argument(
    '--gene-annotations',
    type=str,
    default=None,
    help='(optional) Annotations for genes to use for plotting (must contain a column named `gene_id`)'
)
parser.add_argument(
    '--genome-annotations',
    type=str,
    default=None,
    help='(optional) Annotations for genomes to use for plotting (must contain a column named `genome_id`)'
)

# Parse the arguments
args = parser.parse_args()

# Read and format the data
data = read_data(args.__dict__)

# Define the menu items to be presented to the user
# The `menus` object is a list, which organizes the menus into tabs
menus = [
    # Second level is a dict, which defines the content of each menu tab
    # This example only has a single tab, but each additional dict
    # will add another tab to the menu display
    dict(
        # Label to be displayed at the top of the tab
        label="Customize Display",
        params=[
            # Filter alignments by a minimum percent identity (similarity)
            dict(
                # ID used to access the value of this menu item
                elem_id="minimum-pctid",
                # Define the type of menu item
                type="input",
                # Define the value type
                input_type="number",
                # Label displayed along this menu item
                label="Minimum Alignment Identity",
                # Default value
                value=90.,
            ),
            # Filter alignments by a minimum coverage percentage
            dict(
                elem_id="minimum-coverage",
                # Define the type of menu item
                type="input",
                # Define the value type
                input_type="number",
                # Label displayed along this menu item
                label="Minimum Alignment Coverage",
                # Default value
                value=90.,
            ),
            # Filter genomes by a minimum number of genes aligned
            dict(
                elem_id="minimum-genes-per-genome",
                type="input",
                input_type="number",
                label="Minimum Number of Genes to Display Genome",
                value=1,
            ),
            # Allow the user to color the genes by the
            # coverage or identity of the alignment, as well
            # as any additional user-provided metadata
            dict(
                elem_id="color-genes-by",
                label="Color Genes By",
                type="dropdown",
                options=data['available_gene_annotations'],
                value="pident",
            ),
            # Set up the labels for each gene
            dict(
                elem_id="label-genes-by",
                label="Label Genes By",
                type="dropdown",
                options=data['available_gene_labels'],
                value="",
            ),
            # Limit the length of each gene label
            dict(
                elem_id="max-gene-label-len",
                label="Maximum Gene Label Length",
                type="input",
                input_type="numeric",
                value=60,
                # Make sure that this option stays in the same column as the previous
                keep_with_previous=True
            ),
            # Set up the labels for each genome
            dict(
                elem_id="label-genomes-by",
                label="Label Genomes By",
                type="dropdown",
                options=data['available_genome_labels'],
                value="",
            ),
            # Limit the length of each genome label
            dict(
                elem_id="max-genome-label-len",
                label="Maximum Genome Label Length",
                type="input",
                input_type="numeric",
                value=60,
                # Make sure that this option stays in the same column as the previous
                keep_with_previous=True
            ),
            # Set the colorscale used for the heatmap
            dict(
                elem_id="heatmap-colorscale",
                label="Heatmap Color Scale",
                type="dropdown",
                options=[
                    dict(label=v, value=v)
                    for v in px.colors.named_colorscales()
                ],
                value="blues",
            ),
            # Set the width of the tree
            dict(
                elem_id="tree-width",
                label="Tree Width",
                type="slider",
                min_val=0.1,
                max_val=0.9,
                value=0.4,
                step=0.01,
            ),
            # Set the width of the figure
            dict(
                elem_id="figure-width",
                label="Figure Width",
                type="slider",
                min_val=200,
                max_val=2400,
                value=800,
                step=20,
            ),
            dict(
                elem_id="figure-height",
                label="Figure Height",
                type="slider",
                min_val=200,
                max_val=2400,
                value=650,
                step=20,
            ),
        ]
    ),
]

# Keep a list of params that were considered, but aren't currently being used
deprecated_menu_items = [
    # Show a heatmap or a t-SNE map
    dict(
        # ID used to access the value of this menu item
        elem_id="display-type",
        # Label displayed along this menu item
        label="Display Type",
        # Dropdown
        type="dropdown",
        # Available options
        options=[
            dict(
                label="Heatmap + Tree",
                value="heatmap"
            ),
            dict(
                label="t-SNE Map",
                value="tsne"
            ),
        ],
        # Default value
        value="heatmap"
    ),
    # Allow the user to set a title to the plot
    dict(
        # ID used to access the value of this menu item
        elem_id="plot-title",
        # Label displayed along this menu item
        label="Plot Title",
        # Free-form input box
        type="input",
        # Input must be a string
        input_type="string",
        # Default value
        value="",
    ),
]


@lru_cache(maxsize=1)
def get_mask(min_pctid, min_cov):

    # Compute the mask for which genes/genomes pass the filter
    return (data["alignments_pident"] >= min_pctid) & \
           (data["alignments_coverage"] >= min_cov)


# Generate a wide table of alignments based on minimum thresholds
def format_alignments_wide(min_pctid, min_cov, display_value):

    # The options for `display_value` are:
    #   pident
    #   coverage
    #   description
    #   mask (returns a bool for each cell reflecting the filter)
    assert display_value in ['pident', 'coverage', 'description', 'mask']

    # Compute the mask for which genes/genomes pass the filter
    # Using a subfunction allows us to cache the value to help with
    # multiple calls to the format_alignments_wide() function
    mask = get_mask(min_pctid, min_cov)

    # If that is all we need to return
    if display_value == "mask":

        # Return the DataFrame of bools
        return mask

    # Otherwise:
    else:

        # Format the key for the source data
        data_key = f"alignments_{display_value}"

        # Make sure that the key is valid
        assert data_key in list(data.keys())

        # Return the alignments which pass the filter
        return data[data_key].where(mask)


def plot_gig_map(_, selections):
    """Render the gig-map display based on the data and the user's menu selections."""

    # Render that figure
    return plot_gig_map_heatmap(selections)

    # DEPRECATED BEHAVIOR BELOW WAS TO SUPPORT T-SNE
    # # The two options for display are 'heatmap' and 'tnse'
    # assert selections["display-type"] in ['heatmap', 'tsne']

    # # If the user selected the option to display a heatmap + tree
    # if selections["display-type"] == "heatmap":

    #     # Render that figure
    #     return plot_gig_map_heatmap(selections)

    # # Otherwise
    # else:

    #     assert selections["display-type"] == 'tsne'

    #     # Render that figure
    #     return plot_gig_map_tsne(selections)

def get_gene_annot_values(selections):
    """Return a dict with the specified annotation for 'color-genes-by'."""

    # Make sure that the specified value is in the annotation table
    assert selections["color-genes-by"] in data["gene_annotations"].columns.values

    # Make a dict with numeric values
    value_dict = data[
        "gene_annotations"
    ][
        selections["color-genes-by"]
    ].apply(
        lambda v: pd.to_numeric(v, errors="coerce")
    ).dropna(
    ).to_dict()

    # Make sure that there are at least some numeric values
    msg = f"Column {selections['color-genes-by']} contains no numeric values"
    assert len(value_dict) > 0, msg

    return value_dict

def plot_gig_map_tsne(selections):
    """Render the t-SNE map display."""

    # Read the t-SNE coordinates
    tsne = data["tsne"]

    # If there are gene annotations
    if "gene_annotations" in data:

        # For all of the possible annotations
        for col_name, col_values in data["gene_annotations"].items():

            # Add the values as a column to the table to display
            tsne = tsne.assign(
                **{
                    col_name: col_values.apply(
                        # Limit the length of each label, if the label is a string
                        lambda v: v[:selections["max-gene-label-len"]] if isinstance(v, str) else v
                    )
                }
            )

    # Genes cannot be colored by genome-specific alignment data
    # If the user has selected another piece of metadata to color by
    if selections["color-genes-by"] not in [None, "pident", "coverage"]:

        # Set the name of the column to use for colors
        color_by_column=selections['color-genes-by']

    else:

        # Set a null value for the color_by_column
        color_by_column = None

    fig = px.scatter(
        data_frame=tsne.reset_index(),
        x='t-SNE 1',
        y='t-SNE 2',
        hover_name='index',
        color=color_by_column,
        hover_data=tsne.columns.values
    )

    # Set up the layout
    fig.update_layout(
        # White background
        paper_bgcolor='white',
        plot_bgcolor='white',
        # Figure height and width
        height=selections['figure-height'],
        width=selections['figure-width'],
    )

    return fig


def plot_gig_map_heatmap(selections):
    """Render the heatmap + tree display."""

    # Format a wide table with gene alignments
    # If the user wants to display alignment stats
    if selections["color-genes-by"] in ["pident", "coverage"]:

        # Make a table with the alignment stats
        plot_df = format_alignments_wide(
            selections["minimum-pctid"],
            selections["minimum-coverage"],
            selections["color-genes-by"]
        )

    # Otherwise, the user wants to color by gene annotation
    else:

        # Get the values for each gene as a dict
        value_dict = get_gene_annot_values(selections)

        # Start with a table filtered by alignment characteristics, and filled in with True/False
        plot_df = format_alignments_wide(
            selections["minimum-pctid"],
            selections["minimum-coverage"],
            "mask"
        # Now replace any non-null value (for which an alignment is present) with the annotation
        ).replace(
            to_replace={
                gene_name: {
                    True: gene_value,
                    False: None
                }
                for gene_name, gene_value in value_dict.items()
            }
        )

    # Drop any genomes which don't have the minimum number of alignments
    plot_df = plot_df.loc[
        plot_df.notnull().sum(axis=1) >= selections["minimum-genes-per-genome"]
    ]

    # Drop any genes which don't have any alignments
    plot_df = plot_df.loc[
        :,
        plot_df.notnull().sum() > 0
    ]

    # Make sure that there are at least two genomes with alignments
    assert plot_df.shape[0] > 1, "<=1 genome with alignments passing filter"

    # Make a table with the text descriptions of each alignment
    text_df = format_alignments_wide(
        selections["minimum-pctid"],
        selections["minimum-coverage"],
        "description"
    ).reindex(
        index=plot_df.index.values,
        columns=plot_df.columns.values,
    ).fillna(
        "No alignments found"
    )

    # Create a tree using the set of genomes which contain alignments
    node_positions = make_nj_tree(plot_df.index.values, data['distances'])

    # The figure will render with a dendrogram on the left and a heatmap on the right

    # Set up a base level figure
    fig = go.Figure()

    # Render the tree with multiple traces
    for trace in plot_tree(
        node_positions,
        selections,
        data,
        xaxis="x",  # Primary X axis
        yaxis="y",  # Primary Y axis
    ):
        fig.add_trace(
            trace
        )

    # Render the heatmap
    fig.add_trace(
        plot_heatmap(
            dict(
                values=plot_df,
                text=text_df
            ),
            node_positions,
            selections,
            data,
            xaxis="x2",  # Secondary X axis
            yaxis="y",   # Primary Y axis
        )
    )

    # Render the colorbar
    fig.add_trace(
        plot_colorbar(
            min_val=plot_df.fillna(100).min().min(),
            max_val=plot_df.fillna(0).max().max(),
            color_genes_by=selections["color-genes-by"],
            label={
                i['value']: i['label']
                for i in data["available_gene_annotations"]
            }[
                selections["color-genes-by"]
            ],
            colorscale=selections["heatmap-colorscale"],
            xaxis="x3",  # Tertiary X axis
            yaxis="y2",  # Secondary Y axis
        )
    )

    # Set up the labels for the genomes to add to the axis
    genome_labels = list(node_positions.genome_order)

    # If the user elected to label the genomes by something other than their ID
    if selections["label-genomes-by"] != "":

        # Map the labels to the values from the annotation table
        genome_labels = list(map(
            data["genome_annotations"][selections["label-genomes-by"]].get,
            genome_labels
        ))

    # Set up the layout
    fig.update_layout(
        # Set up the primary x-axis (with the tree)
        xaxis=dict(
            title_text="Genome Distance (ANI)",
            domain=[0, selections["tree-width"]],
            range=[
                node_positions.df['x'].max() * -0.01,
                node_positions.df['x'].max() * 1.01,
            ]
        ),
        # Primary y-axis (with the tree)
        yaxis=dict(
            tickmode="array",
            tickvals=list(range(len(genome_labels))),
            ticktext=[
                l[:selections["max-genome-label-len"]]
                for l in genome_labels
            ],
            side="right",
            anchor="x3",
            showticklabels=True,
            domain=[0, 0.9]
        ),
        # Secondary y-axis (with the colorbar)
        yaxis2=dict(
            anchor="x2",
            showticklabels=True,
            domain=[0.91, 1.0]
        ),
        # Secondary x-axis (with the heatmap),
        xaxis2=dict(
            domain=[selections["tree-width"], 1.0],
        ),
        # Tertiary x-axis (with the colorbar),
        xaxis3=dict(
            domain=[selections["tree-width"], 1.0],
            anchor="y2",
            side="top",
        ),
        paper_bgcolor='white',
        plot_bgcolor='white',
        # Figure height and width
        height=selections['figure-height'],
        width=selections['figure-width'],
    )

    return fig

# Instantiate the MenuDrivenFigure object
mdf = MenuDrivenFigure(
    data=data,
    menus=menus,
    function=plot_gig_map,
    title="Genes in Genomes Map"
)

# Launch the Dash/Flask app
mdf.run_server(
    host='0.0.0.0',
    port=8080,
    debug=True,
)