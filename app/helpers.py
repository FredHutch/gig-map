from cartesian_tree import make_nj_tree
from direct_redis import DirectRedis
import logging
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from redis.exceptions import BusyLoadingError
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from time import sleep


####################
# HELPER FUNCTIONS #
####################

def remove_genome_file_ext(fp):
    for ext in ['.gz', '.fna', '.fasta', '.fa']:
        if fp.endswith(ext):
            fp = fp[:-len(ext)]
    return fp


def read_alignments(r):
    """Read in the alignment information from redis"""

    # Get the logger
    logger = logging.getLogger('gig-map')

    # Try to get the list of keys used to store slices of the alignments
    logger.info("Reading alignment shard keys")
    alignments_keys = r.get("alignments_keys")

    # If there are no keys
    if alignments_keys is None:

        # Just read in all of the alignments directly
        logger.info("Reading unsharded alignments")
        return r.get("alignments")

    # Otherwise, if there are keys
    else:

        logger.info("Reading alignments in shards")
        # Read in all of the chunks and combine them
        return pd.concat(map(r.get, alignments_keys))


def read_data(args):
    """Read in the data needed to render the heatmap."""

    # Get the logger
    logger = logging.getLogger('gig-map')

    # Return data formatted as a dict
    output = dict()

    # The default annotations available for each gene are the
    # alignment identity and coverage
    output["available_gene_annotations"] = [
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
    output["available_gene_labels"] = []

    # Open the redis store for reading
    logger.info(f"Connecting to redis at {args['host']}:{args['port']}")
    with DirectRedis(host=args['host'], port=args['port']) as r:

        # Wait for the database to load - up to 120 seconds
        for _ in range(120):

            # Try to load the keys
            try:
                _ = r.keys()
                break

            # Catch the BusyLoadingError
            except BusyLoadingError:
                logger.info("Redis is loading...")
                sleep(1)

        # Log the available keys, to help with debugging
        logger.info("Available keys in redis:")
        all_keys = r.keys()
        for k in all_keys:
            logger.info(k)

        # Read in the alignment information in long format
        logger.info("Reading alignments")
        output["alignments"] = read_alignments(r)
        assert output["alignments"] is not None

        # Get the mapping of genome_ix to filenames
        logger.info("Reading genome index")
        output["genome_ix"] = r.get("genome_ix")
        assert output["genome_ix"] is not None

        # Add the genome names to the alignments table
        output["alignments"] = output["alignments"].assign(
            genome_name = output["alignments"].genome_ix.apply(
                lambda i: output["genome_ix"][i]
            )
        )

        # Get the mapping of gene_ix to strings
        logger.info("Reading gene index")
        output["gene_ix"] = r.get("gene_ix")
        assert output["gene_ix"] is not None

        # Read the genome groups which have been pre-computed
        # Results will be output to keys in data formatted as:
        # "genome_groups {threshold}"
        # "genome_group_df {threshold}"
        # "group_distance_df {threshold}"
        # Keep track of all possible values for `threshold`
        output["genome_group_thresholds"] = set([])
        for key in all_keys:
            if key.startswith(("genome_groups ", "genome_group_df ", "group_distance_df ")):
                # Save the data
                output[key] = r.get(key)
                # Parse the threshold value
                threshold_value = int(key.rsplit(" ", 1)[-1])
                # Add the threshold value
                output["genome_group_thresholds"].add(threshold_value)

        # Format output["genome_group_thresholds"] as a list
        output["genome_group_thresholds"] = list(output["genome_group_thresholds"])

        # If marker genes were used, there will be a list in `marker_genes`
        output["marker_genes"] = r.get("marker_genes") if "marker_genes" in all_keys else []

        # If there were any marker genes used
        if len(output["marker_genes"]) > 0:

            # Set up a nested dict to use for the marker gene groups
            output["marker_gene_groups"] = dict()

            # Iterate over each marker gene
            for marker_gene in output["marker_genes"]:

                # Set up a dict
                output["marker_gene_groups"][marker_gene] = dict()

                # Iterate over each threshold value
                for threshold_value in output["genome_group_thresholds"]:

                    # Set up a dict
                    output["marker_gene_groups"][marker_gene][threshold_value] = dict()

                    # Iterate over each individual table
                    for table_name in ["genome_groups", "genome_group_df", "group_distance_df"]:

                        # Format the key based on the marker gene and the threshold
                        key = f"{marker_gene} {table_name} {threshold_value}"

                        # Make sure the key is present
                        assert key in all_keys, f"Expected to find key {key}"

                        # Read the table
                        output["marker_gene_groups"][marker_gene][threshold_value][table_name] = r.get(key)
        
        # Sort the keys numerically
        output["genome_group_thresholds"].sort()

        # Read in the t-SNE coordinates
        output["tsne"] = r.get("tsne")

        #####################
        # DISTANCE MATRICES #
        #####################

        # Set up a dict with all of the distance matrices
        output["distances"] = dict()

        # Read in the ANI pairwise genome distances, merging multiple shards
        output["distances"]["ani"] = pd.concat(
            [
                r.get(redis_key)
                for redis_key in r.get("distances_keys")
            ]
        )
        logger.info("Read in ANI distances")

        # Read in the distances for each marker gene (if any)
        for marker_gene in output["marker_genes"]:

            # Format the key which has the list of shards to read in
            marker_gene_distances_key = f"distances_keys {marker_gene}"

            # Make sure that the key is present
            msg = f"Could not find pairwise distances for {marker_gene}"
            assert marker_gene_distances_key in all_keys, msg

            # Use the name of the marker gene as the key of the dict
            output["distances"][marker_gene] = pd.concat(
                [
                    r.get(redis_key)
                    for redis_key in r.get(marker_gene_distances_key)
                ]
            )
            logger.info(f"Read in distances for {marker_gene}")
        
    # Read in the gene annotations, if any
    if args['gene_annotations'] is None:
        output["gene_annotations"] = None

    # If a table was provided
    else:
        # Read in the table
        logger.info(f"Reading from {args['gene_annotations']}")
        output["gene_annotations"] = pd.read_csv(args['gene_annotations'])

        # Make sure that there is a column named "gene_id"
        msg = "Gene annotation CSV must contain a column named 'gene_id'"
        msg = f"{msg} - found {'; '.join(output['gene_annotations'].columns.values)}"
        assert "gene_id" in output["gene_annotations"].columns.values, msg

        # Set the index of the table as 'gene_id'
        output["gene_annotations"].set_index('gene_id', inplace=True)

        # Add the other columns to the available gene annotations
        for col_name in output["gene_annotations"].columns.values:
            output["available_gene_annotations"].append(
                dict(
                    label=col_name,
                    value=col_name,
                )
            )

            # Also add those columns to the options available for labeling genes
            output["available_gene_labels"].append(
                dict(
                    label=col_name,
                    value=col_name,
                )
            )

    # By default, there are no additional labels to apply to the genomes
    output["available_genome_labels"] = []

    # Read in the genome annotations, if any
    if args['genome_annotations'] is None:
        output['genome_annotations'] = None

    # If a table was provided
    else:
        # Read in the table
        logger.info(f"Reading from {args['genome_annotations']}")
        output['genome_annotations'] = pd.read_csv(args['genome_annotations'])

        # Make sure that there is a column named "genome_id"
        msg = "Genome annotation CSV must contain a column named 'genome_id'"
        msg = f"{msg} - found {'; '.join(output['genome_annotations'] .columns.values)}"
        assert "genome_id" in output['genome_annotations'] .columns.values, msg

        # Set the index of the table as 'genome_id'
        output['genome_annotations'] .set_index('genome_id', inplace=True)

        # Add the other columns to the available genome annotations
        for col_name in output['genome_annotations'] .columns.values:

            # Add those columns to the options available for labeling genomes
            output["available_genome_labels"].append(
                dict(
                    label=col_name,
                    value=col_name,
                )
            )

    logger.info("Done reading all data")

    # Return data formatted as a dict
    return output

def filter_alignments(
        data=None,
        clustering_method=None,
        ani_threshold=None,
        selected_group=None,
        max_n_genomes=None,
        min_pctid=None,
        min_cov=None,
        min_genes_per_genome=None,
        min_genomes_per_gene=None,
        query=None,
):
    """
    Filter the alignments for genomes in the selected group
    based on pident, coverage, and the number of genes aligned per genome.
    """

    # Get the logger
    logger = logging.getLogger('gig-map')

    # If no ANI grouping has been selected
    if ani_threshold == "None":

        # There is no need to filter the alignments
        filtered_alignments = data["alignments"]

        # However, if a marker gene has been selected for clustering
        if clustering_method != "ani":

            # Then we will need to find the subset of genomes
            # which contain the marker gene

            # The genomes which have the marker are in this list
            genomes_with_marker = set(list(data[
                "distances"
            ][
                clustering_method
            ].index.values))

            # Filter the alignments
            filtered_alignments = filtered_alignments.loc[
                filtered_alignments.genome_name.isin(genomes_with_marker)
            ]

    # If an ANI threshold _has_ been selected
    else:

        # Get the genome groups based on this ANI threshold
        genome_groups = read_genome_groups(
            data,
            clustering_method,
            ani_threshold
        )

        # Filter by genome group
        genomes_in_selected_group = set([
            n
            for n, i in genome_groups.items()
            if i == selected_group
        ])
        filtered_alignments = data["alignments"].loc[
            data["alignments"]["genome_name"].isin(genomes_in_selected_group)
        ]
        logger.info(f"Filtered by genome group - {filtered_alignments.shape[0]:,} lines")

    # Filter by alignment quality
    filtered_alignments = filtered_alignments.query(
        f"pident >= {min_pctid}"
    ).query(
        f"coverage >= {min_cov}"
    ).assign(
        mask = True
    )
    logger.info(f"Filtered by alignment quality - {filtered_alignments.shape[0]:,} lines")

    # If the user has provided a `query` string to filter genes based on their annotations
    if query is not None:

        # Get the set of gene indices which pass the filter
        genes_passing_query_filter = filter_genes_by_query(data, query)

        # Subset the alignments to just include those genes
        filtered_alignments = filtered_alignments.loc[
            filtered_alignments.gene_ix.isin(genes_passing_query_filter)
        ]
        logger.info(f"Filtered alignments by query string - {filtered_alignments.shape[0]:,} lines")

    # If the user has decided to filter the genes displayed
    if min_genomes_per_gene > 1:

        # Count up the number of genomes per gene
        genomes_per_gene = filtered_alignments.gene_ix.value_counts()

        # Apply the minimum threshold of the number of genes per genome
        genes_passing_filter = set(genomes_per_gene.index.values[
            genomes_per_gene >= min_genomes_per_gene
        ])
        genes_passing_filter = set(list(genes_passing_filter))

        logger.info(f"{len(genes_passing_filter):,} genes")

        # Filter the alignments to only include those genes which pass the filter
        filtered_alignments = filtered_alignments.loc[
            filtered_alignments.gene_ix.isin(genes_passing_filter)
        ]
        logger.info(f"Filtered by final gene set - {filtered_alignments.shape[0]:,} lines")

    # Only keep a single alignment per gene / genome
    filtered_alignments = filtered_alignments.groupby(
        ["genome_ix", "gene_ix"]
    ).head(1)
    logger.info(f"Deduplicated alignments - {filtered_alignments.shape[0]:,} lines")

    # Count up the number of genes per genome
    genes_per_genome = filtered_alignments.genome_ix.value_counts()

    # If the user has specified a maximum number of genomes to display
    if max_n_genomes is not None:

        # Apply the maximum number of genomes to display
        genes_per_genome = genes_per_genome.head(max_n_genomes)

    # Apply the minimum threshold of the number of genes per genome
    genomes_passing_filter = set(genes_per_genome.index.values[
        genes_per_genome >= min_genes_per_genome
    ])
    genomes_passing_filter = set(list(genomes_passing_filter))

    logger.info(f"{len(genomes_passing_filter):,} genomes")

    # Filter the alignments to only include those genomes which pass the filter
    filtered_alignments = filtered_alignments.loc[
        filtered_alignments.genome_ix.isin(genomes_passing_filter)
    ]
    logger.info(f"Filtered by final genome set - {filtered_alignments.shape[0]:,} lines")

    return filtered_alignments


def filter_genes_by_query(data, query):
    """Return a set of `gene_ix` for those genes whose annotations satisfy the query string."""

    # Get the logger
    logger = logging.getLogger("gig-map")
    logger.info(f"Applying filter to genes: {query}")

    # There must be a gene annotation table
    msg = "A query string can only be applied with --gene-annotations"
    assert data["gene_annotations"] is not None, msg

    # Apply the query string to the gene annotation table
    filtered_genes = data["gene_annotations"].query(query)

    # At least one gene must pass the filter
    msg = f"Query string has 0 results: {query}"
    assert filtered_genes.shape[0] > 0, msg

    # Count up the total and filtered number of genes
    n_tot = data["gene_annotations"].shape[0]
    n_fil = filtered_genes.shape[0]

    logger.info(f"Number of annotated genes passing the filter: {n_fil:,} / {n_tot:,}")

    # Get the set of `gene_id` strings for those genes passing the query filter
    filtered_gene_ids = set(filtered_genes.index.values)

    # Return the set of `gene_ix` values for the filtered genes
    return {
        i
        for i, n in enumerate(data["gene_ix"])
        if n in filtered_gene_ids
    }


def get_gene_annot_values(data, selections):
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


def order_genes(data, plot_df, method="ward", metric="euclidean"):
    """Reorder the genes based on linkage clustering."""

    # Get the logger
    logger = logging.getLogger("gig-map")

    logger.info("Clustering genes by genome assignment")

    # Make a rectangular matrix of gene pident across genomes
    pident_df = plot_df.pivot(
        index="gene_ix",
        columns="genome_name",
        values="pident"
    ).fillna(
        0
    ).sort_index(
        # Sort on the existing `gene_ix`
        axis=0
    )

    logger.info(f"Sorting {pident_df.shape[0]:,} genes using data from {pident_df.shape[1]:,} genomes")

    # Calculate the pairwise distances between genes
    dists = pdist(pident_df.values, metric)

    # Perform linkage clustering
    Z = hierarchy.linkage(
        dists,
        method=method
    )

    # Order the tree so that adjacent leaves are more similar
    Z_ordered = hierarchy.optimal_leaf_ordering(
        Z,
        dists
    )

    # Get the ordered list of leaves
    leaves_list = hierarchy.leaves_list(Z_ordered)

    # Map the previous `gene_ix` to the new `gene_ix`
    gene_map = {
        pident_df.index.values[i]: i
        for i in leaves_list
    }

    # Replace the `gene_ix` values in `plot_df`
    plot_df = plot_df.assign(
        gene_ix=plot_df.gene_ix.apply(gene_map.get)
    )

    # Replace the `gene_ix` values in `data`
    data["gene_ix"] = [
        data["gene_ix"][pident_df.index.values[i]]
        for i in leaves_list
    ]

    return data, plot_df

def set_figure_height(plot_df, figure_height, height_frac=20, height_padding=200):
    """Dynamically set the figure height based on the number of genomes."""

    # If the user specified an 'auto' figure height
    if figure_height == "auto":

        # Get the total number of genomes
        n_genomes = plot_df["genome_ix"].unique().shape[0]

        # Apply the following formula to give each genome sufficient space
        return int((n_genomes * height_frac) + height_padding)

    # Otherwise, the user must have specified a number of pixels
    else:

        # Convert `figure_height` to an int
        return int(figure_height)

def plot_heatmap(data, plot_df, node_positions, selections, xaxis='x', yaxis='y'):
    """Return a heatmap rendered from the gene DataFrame and the genome positions."""

    # Format a set of wide tables
    tables = {
        table_key: plot_df.pivot(
            index="genome_name",
            columns="gene_ix",
            values=column_key
        ).sort_index(
            # Sorting on the column index preserves the t-SNE order
            axis=1
        ).reindex( 
            # Reorder the rows (genomes) to match the tree
            index=node_positions.genome_order,
        ).rename(
            # Transform the gene index to the gene ID string
            columns=lambda i: data["gene_ix"][i]
        )
        for table_key, column_key in [
            ("values", selections["color-genes-by"]),
            ("text", "description"),
        ]
    }

    # If the user elected to label the genes by something other than their ID
    if selections["label-genes-by"] != "":

        # Function to rename genes
        def format_gene_id(gene_id):

            # Get the value
            gene_label = data[
                "gene_annotations"
            ][
                selections["label-genes-by"]
            ].get(gene_id)

            # If there is no value
            if gene_label is None:

                # Just show the gene ID
                return gene_id

            # If there is a value
            else:

                # Join together the gene ID and the annotated label,
                #  up to a maximum length
                return f"{gene_label} ({gene_id})"[
                    :selections["max-gene-label-len"]
                ]

        # Rename the columns of the DataFrame
        for k in ["values", "text"]:
            tables[k] = tables[k].rename(
                columns=format_gene_id
            )

    return go.Heatmap(
        x=list(tables["values"].columns.values),
        z=tables["values"].values,
        text=tables["text"].values,
        xaxis=xaxis,
        yaxis=yaxis,
        colorscale=selections["heatmap-colorscale"],
        showscale=False,
        hovertemplate="%{text}<extra></extra>",
    )


def plot_tree(data, node_positions, selections, xaxis='x', yaxis='y', try_relabel=True):
    """Return a Plotly trace rendered from a skbio tree."""

    # If the user decided to label the genomes
    if try_relabel and selections["label-genomes-by"] != "":

        # Replace the values in the 'name' column of the DataFrame used for plotting
        label_genomes_dict = data[
            "genome_annotations"
        ][
            selections["label-genomes-by"]
        ].to_dict()

        # If there is no annotation, keep the original label
        node_positions.df = node_positions.df.replace(
            to_replace=dict(
                name=lambda l: label_genomes_dict.get(l, l)
            )
        )

    # Return a ScatterGL
    yield go.Scattergl(
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

    # Also show a set of dotted lines extending each tip
    yield go.Scattergl(
        showlegend=False,
        mode="lines",
        x=node_positions.extension_x_coords(),
        y=node_positions.extension_y_coords(),
        hoverinfo="skip",
        line=dict(
            color="black",
            dash="dot",
            width=1,
        ),
        xaxis=xaxis,
        yaxis=yaxis,
    )


def render_gig_map_heatmap(
    data,
    plot_df,
    selections,
    selected_group
):
    """Render a plotly figure showing a gig-map display based on plot selections."""

    # Get the logger
    logger = logging.getLogger('gig-map')

    # If the user wants to color by gene annotation
    if selections["color-genes-by"] not in ["pident", "coverage"]:

        # Get the values for each gene as a dict
        logger.info(f"Getting values for {selections['color-genes-by']}")
        value_dict = get_gene_annot_values(data, selections)

        # Add that column to the DataFrame
        logger.info(f"Adding column for {selections['color-genes-by']}")
        plot_df = plot_df.assign(
            **{
                selections["color-genes-by"]: plot_df["gene_ix"].apply(
                    lambda i: value_dict.get(data["gene_ix"][i])
                )
            }
        )

        # Filter down to just those genes the selected annotation
        plot_df = plot_df.loc[
            plot_df[selections["color-genes-by"]].apply(pd.notnull)
        ]

    # Count up the number of unique genomes and genes
    logger.info("Getting genome names")
    unique_genomes = list(map(
        lambda i: data["genome_ix"][i],
        plot_df.genome_ix.unique()
    ))
    unique_genes = plot_df.gene_ix.unique()

    logger.info("{:,} genomes, {:,} genes, {:,} alignments".format(
        len(unique_genomes),
        len(unique_genes),
        plot_df.shape[0]
    ))

    # Create a tree using the set of genomes which contain alignments
    node_positions = make_nj_tree(
        unique_genomes,
        data[
            'distances'
        ][
            selections["clustering-method"]
        ]
    )

    logger.info(f"Grouped genomes by {selections['clustering-method']}")

    # The figure will render with a dendrogram on the left and a heatmap on the right

    # Set up a base level figure
    fig = go.Figure()

    # Render the tree with multiple traces
    logger.info("Plotting tree")
    for trace in plot_tree(
        data,
        node_positions,
        selections,
        xaxis="x",  # Primary X axis
        yaxis="y",  # Primary Y axis
    ):
        fig.add_trace(
            trace
        )

    # Render the heatmap
    logger.info("Plotting heatmap")
    fig.add_traces(
        plot_heatmap(
            data,
            plot_df,
            node_positions,
            selections,
            xaxis="x2",  # Secondary X axis
            yaxis="y",   # Primary Y axis
        )
    )

    # Render the colorbar
    logger.info("Plotting colorbar")
    fig.add_trace(
        plot_colorbar(
            min_val=plot_df[
                selections["color-genes-by"]
            ].min(),
            max_val=plot_df[
                selections["color-genes-by"]
            ].max(),
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
            lambda l: data["genome_annotations"][selections["label-genomes-by"]].get(l, l),
            genome_labels
        ))

    # Set up the label which will be used underneath the tree
    if selections["clustering-method"] == "ani":
        tree_xaxis_label = "Genome Distance<br>(ANI)"
    else:
        tree_xaxis_label = f"Marker Distance<br>({selections['clustering-method']})"

    # Set up the plot title
    if selections["ani-threshold"] == "None":
        plot_title=""
    else:
        plot_title=f"Genome Group {selected_group}"

    # Set the height of the colorbar dynamically based on the total figure size
    # so that it always takes up about 50 pixels
    colorbar_height = 50 / selections["figure-height"]

    # Set up the layout
    logger.info("Updating layout")
    fig.update_layout(
        # Set up the title of the plot
        title=plot_title,
        # Set up the primary x-axis (with the tree)
        xaxis=dict(
            title_text=tree_xaxis_label,
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
            domain=[0, 0.99 - colorbar_height],
            automargin=True,
        ),
        # Secondary y-axis (with the colorbar)
        yaxis2=dict(
            anchor="x2",
            showticklabels=True,
            domain=[1.0 - colorbar_height, 1.0],
            automargin=True,
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
    logger.info("Rendering plot")

    return fig


def plot_colorbar(
    min_val=0.,
    max_val=100.,
    colorscale="blues",
    color_genes_by="pctid",
    label="Percent Identity",
    xaxis="x3",
    yaxis="y2",
):
    """Render a colorbar as a heatmap on a dedicated axis."""

    value_list = [
        v
        for v in np.linspace(
            min_val, max_val, num=100
        )
    ]

    return go.Heatmap(
        y=[label],
        x=value_list,
        z=[value_list],
        xaxis=xaxis,
        yaxis=yaxis,
        colorscale=colorscale,
        showscale=False,
        hovertemplate="%{z}<extra></extra>",
    )


def read_genome_groups(data, clustering_method, ani_threshold):
    """Read in the assignment of genomes to groups at a particular ANI threshold."""

    if clustering_method == "ani":
        return data[
            f"genome_groups {ani_threshold}"
        ].set_index(
            "genome"
        )[
            "group"
        ]
    else:
        return data[
            "marker_gene_groups"
        ][
            clustering_method
        ][
            ani_threshold
        ][
            "genome_groups"
        ].set_index(
            "genome"
        )[
            "group"
        ]


def empty_plot():
    return go.Figure(
        layout=dict(
            xaxis=dict(
                visible=False
            ),
            yaxis=dict(
                visible=False
            ),
            plot_bgcolor="rgba(0, 0, 0, 0)",
            paper_bgcolor="rgba(0, 0, 0, 0)",
            height=10,
            width=10,
        )
    )