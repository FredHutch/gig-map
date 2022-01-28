from cartesian_tree import get_node_positions
from cartesian_tree import make_nj_tree
import gzip
import pandas as pd
import plotly.graph_objects as go


def plot_genome_tree(fig):
    """Plot a genome tree on the FigureBuilder."""

    # Define the position at which this panel will be rendered
    fig.log("Adding subplot for the genome tree")
    fig.subplots.add(
        # ID for the subplot
        id='genome_tree',
        # Ordinal position on the horizontal axis
        x_index=0,
        # Ordinal position on the vertial axis
        y_index=0,
        # Share the y-coordinates with all other subplots
        share_y=True,
    )
    
    # Make an object to map the neighbor joining tree on a cartesian plot
    fig.log("Formatting genome tree layout")
    node_positions = make_nj_tree(
        # Note that this data object was read in by read_genome_alignments()
        # in the 'genomeHeatmap' element of the FigureBuilder
        fig.data['genomeHeatmap']['genomes_to_plot'],
        fig.data['genomeTree']['dm']
    )

    # If the user has elected to rename the genomes in the plot
    if fig.params['global']['label_genomes_by'] is not None:

        # Reference the name of the column being used to label
        label_col = fig.params['global']['label_genomes_by']
        fig.log(f"Renaming genomes by {label_col}")

        # Reference the CSV used to annotate genomes
        genome_annotations = fig.data[
            "global"
        ][
            "genome_annotations"
        ]

        # A genome annotation table must have been provided
        msg = "Must provide --genome-annotations with --label-genomes-by"
        assert genome_annotations is not None, msg

        # The DataFrame must contain a column which matches the name
        msg = f"Cannot find a columns named {label_col} in {fig.params['genome-annotations']}"
        assert label_col in genome_annotations.columns.values, msg

        # Replace the values in the 'name' column of the DataFrame used for plotting
        label_genomes_dict = genome_annotations[label_col].to_dict()

        # If there is no annotation, keep the original label
        node_positions.genome_order = list(map(
            lambda l: label_genomes_dict.get(l, l),
            node_positions.genome_order
        ))

    # Add the traces to the plot

    # Scatterplot for the nodes
    fig.log("Plotting nodes for genome tree")
    fig.subplots.plot(
        id="genome_tree",
        trace=go.Scattergl(
            name="Neighbor Joining Tree",
            showlegend=False,
            mode="lines",
            x=node_positions.x_coords(),
            y=node_positions.y_coords(),
            text=node_positions.text(),
            hoverinfo="text",
        )
    )

    # Also show a set of dotted lines extending each tip
    fig.log("Plotting lines for genome tree")
    fig.subplots.plot(
        id="genome_tree",
        trace=go.Scattergl(
            showlegend=False,
            mode="lines",
            x=node_positions.extension_x_coords(),
            y=node_positions.extension_y_coords(),
            hoverinfo="skip",
            line=dict(
                color="black",
                dash="dot",
                width=1,
            )
        )
    )

    # Format the axes for this plot
    fig.log("Formatting axes for genome tree")
    fig.subplots.format_axis(
        id="genome_tree",
        ax="x",
        params=dict(
            range=[
                node_positions.df['x'].max() * -0.01,
                node_positions.df['x'].max() * 1.01,
            ]
        )
    )

    fig.subplots.format_axis(
        id="genome_tree",
        ax="y",
        params=dict(
            tickmode="array",
            tickvals=list(range(len(node_positions.genome_order))),
            ticktext=[
                l[:fig.params["global"]["max_genome_label_len"]]
                for l in node_positions.genome_order
            ],
            side="right",
            showticklabels=True,
            automargin=True,
        ),
        anchor="genome_tree"
    )

    fig.log("Done plotting genome tree")


def read_global(
    genome_annotations=None,
    gene_annotations=None,
    gene_order=None,
    **kwargs
):
    """Read in the data for gig-map in the global namespace."""

    # Format output as a dict
    data = dict()

    # Read in the genome annotations
    if genome_annotations is None:
        data['genome_annotations'] = None
    else:
        data['genome_annotations'] = pd.read_csv(genome_annotations)

        # The table must have a column named `genome_id`
        msg = "--genome-annotations must contain column 'genome_id'"
        assert 'genome_id' in data['genome_annotations'].columns.values, msg
        data['genome_annotations'].set_index("genome_id", inplace=True)

    # Read in the gene annotations
    if gene_annotations is None:
        data['gene_annotations'] = None
    else:
        data['gene_annotations'] = pd.read_csv(gene_annotations)

        # The table must have a column named `genome_id`
        msg = "--gene-annotations must contain column 'gene_id'"
        assert 'gene_id' in data['gene_annotations'].columns.values, msg
        data['gene_annotations'].set_index("gene_id", inplace=True)

    # Read in the gene order
    with gzip.open(gene_order, 'rt') as handle:
        data['gene_order'] = [
            line.rstrip("\n")
            for line in handle
        ]

    return data


def plot_global(fig):
    """Apply global plot configuration after all of the other params are set."""

    # Set the global layout elements
    fig.subplots.fig.update_layout(
        paper_bgcolor='white',
        plot_bgcolor='white'
    )


def read_genome_tree(distmat=None, **kwargs):
    """Read in data needed for the genome tree."""

    # Read in the distance matrix
    dm = pd.read_csv(distmat, index_col=0)

    # Format a CartesianTree object
    data = {
        "dm": dm
    }

    return data


def read_genome_alignments(csv=None):
    """Read in data needed for the genome alignments."""

    # Read in the table of alignments
    alignments = pd.read_csv(csv)

    # Get the unique list of genomes which have alignments
    genomes_to_plot = set(alignments["genome"].tolist())

    return {
        "alignments": alignments,
        "genomes_to_plot": genomes_to_plot
    }