#!/usr/bin/env python3

# Import the helper functions
from os import read
from app.helpers import read_data, make_nj_tree, plot_tree

# Import the menu-driven-figure library
from menu_driven_figure.app import MenuDrivenFigure
import argparse
from functools import lru_cache
import logging
import pandas as pd
from plotly.subplots import make_subplots
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
    '--alignments',
    type=str,
    required=True,
    help='Table of gene alignments (suffix: .csv.gz)'
)
parser.add_argument(
    '--distances',
    type=str,
    required=True,
    help='Pairwise genome ANI values (suffix: .dists.tsv.gz)'
)
parser.add_argument(
    '--gene-annotations',
    type=str,
    default=None,
    help='(optional) Annotations for genes to use for plotting (must contain a column named `gene_id`)'
)

# Parse the arguments
args = parser.parse_args()

# Read and format the data
data = read_data(
    args.alignments,
    args.distances,
    args.gene_annotations
)

# Define the menu items to be presented to the user
# The `menus` object is a list, which organizes the menus into tabs
menus = [
    # Second level is a dict, which defines the content of each menu tab
    # This example only has a single tab, but each additional dict
    # will add another tab to the menu display
    dict(
        # Label to be displayed at the top of the tab
        label="Filter Alignments",
        params=[
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
            )
        ]
    ),
    dict(
        # Label to be displayed at the top of the tab
        label="Customize Display",
        params=[
            # Allow the user to color the genes by the
            # coverage or identity of the alignment, as well
            # as any additional user-provided metadata
            dict(
                # ID used to access the value of this menu item
                elem_id="color-genes-by",
                # Label displayed along this menu item
                label="Color Genes By",
                # Dropdown
                type="dropdown",
                # Available options
                options=data['available_gene_annotations'],
                # Default value
                value="pident",
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
            )
        ]
    ),
]

# Filter alignments based on minimum thresholds
@lru_cache(maxsize=128)
def filter_alignments(min_pctid, min_cov):

    return data['alignments'].query(
        f"pident >= {min_pctid}"
    ).query(
        f"coverage >= {min_cov}"
    )


# Generate a wide table of alignments based on minimum thresholds
@lru_cache(maxsize=128)
def format_alignments_wide(min_pctid, min_cov, display_value):

    return filter_alignments(
        min_pctid,
        min_cov
    ).pivot_table(
        index="genome",
        columns="sseqid",
        values=display_value
    )

def plot_gig_map(_, selections):
    """Render the gig-map display based on the data and the user's menu selections."""

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

        # Start with a table filtered by alignment characteristics, and filled in with pident
        plot_df = format_alignments_wide(
            selections["minimum-pctid"],
            selections["minimum-coverage"],
            "pident"
        # Now replace any non-null value (for which an alignment is present) with the annotation
        ).apply(
            lambda c: c.apply(lambda v: None if pd.isnull(v) else value_dict.get(c.name))
        )

    # Make sure that >1 genome contains alignments which pass the threshold
    assert plot_df.shape[0] > 1, "Not enough genomes have alignments passing the filter"

    # Create a tree using the set of genomes which contain alignments
    tree = make_nj_tree(plot_df.index.values, data['dists'])

    # The figure will render with a dendrogram on the left and a heatmap on the right
    fig = make_subplots(
        rows=1,
        cols=2,
        shared_yaxes=True
    )

    print(plot_tree(tree))

    # Render the tree
    fig.add_trace(
        plot_tree(tree)
    )

    # Render the heatmap


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