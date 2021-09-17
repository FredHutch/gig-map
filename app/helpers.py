import logging
import pandas as pd
from direct_redis import DirectRedis
from redis.exceptions import BusyLoadingError
import plotly.graph_objects as go
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