#!/usr/bin/env python3

# Import libraries
import gzip
import logging
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from sklearn.cluster import KMeans
import sys
from time import time
from uuid import uuid4


class OrderRows:
    """
    Order the rows of a DataFrame with linkage clustering, breaking
    up the rows initially into k clusters based on the total counts.
    """
    
    def __init__(
        self,
        df:pd.DataFrame,
        method="average",
        metric="dice",
        logger=None,
        k:int=5
    ):

        # Attach the data to the object
        self.df = df
        self.row_n = df.shape[0]

        # Save the parameters
        self.method = method
        self.metric = metric
        self.start_time = None
        self.k = k

        # Attach the logger
        self.logger = logger

        # Break up the rows into clusters with similar counts
        self.group_rows_kmeans()

        # Assign an order within each of those clusters
        self.order_rows_in_groups()

        # Arrange each of the clusters in descending order of counts
        self.calc_cluster_order()

        # Concatenate the clusters to make the overall row order
        self.concatenate_clusters()

    def concatenate_clusters(self):
        """Concatenate the clusters to make the overall row order."""

        self.row_order = [
            row_name
            for cluster_ix in self.cluster_order
            for row_name in self.grouped_row_order[cluster_ix]
        ]

    def calc_cluster_order(self):
        """Arrange each of the clusters in descending order of counts."""

        self.set_timer()

        # Count up the average number of counts for each cluster
        mean_counts = {
            cluster_ix: cluster_df.mean().mean()
            for cluster_ix, cluster_df in self.df.groupby(self.kmeans_labels)
        }

        # Set the order of the clusters with the highest counts first
        self.cluster_order = pd.Series(
            mean_counts
        ).sort_values(
            ascending=False
        ).index.values

        self.stop_timer("Determined cluster order")

    def group_rows_kmeans(self):
        """
        Group rows together by k-means clustering, selecting
        the number of clusters which results in the lowest silhouette score.
        """

        self.set_timer()

        # Compute the number of genomes that each gene is found within
        gene_counts = self.df.sum(axis=1).values

        # Convert to a 2D array
        gene_counts = np.transpose(np.array([gene_counts]))

        assert gene_counts.shape[0] == self.row_n
        assert gene_counts.shape[1] == 1

        # Intialize the object
        kmeans = KMeans(n_clusters=self.k, random_state=0)

        # Fit the clusters
        kmeans.fit(gene_counts)

        # Predict the labels
        self.kmeans_labels = kmeans.predict(gene_counts)

        self.stop_timer("Performed k-means clustering")

    def order_rows_in_groups(self):
        """Assign an order for the rows within each cluster."""

        self.set_timer()
        self.grouped_row_order = {
            group_ix: self.calc_row_order(
                list(row_list.values)
            )
            for group_ix, row_list in pd.Series(
                self.df.index.values
            ).groupby(
                self.kmeans_labels
            )
        }
        self.stop_timer("Ordered rows within all groups")

    def set_timer(self):
        """Set the timer."""
        self.start_time = time()

    def stop_timer(self, msg, decimals=2):
        """Report the elapsed time."""

        # The timer must have been started
        assert self.start_time is not None, "Forgot to start the timer"

        rounded_seconds = round(
            time() - self.start_time,
            decimals
        )

        # Log the message
        self.log(f"{msg} in {rounded_seconds} seconds")

        # Unset the start time variable
        self.start_time = None

    def log(self, msg):
        """Print a message to the logger, if provided, otherwise print."""
        if self.logger is None:
            print(msg)
        else:
            self.logger.info(msg)

    def calc_row_order(self, row_list):
        """Set up the initial row order with linkage clustering."""

        # If there are < 3 rows
        if len(row_list) < 3:

            # No need to change things
            return row_list

        # Perform linkage clustering
        L = hierarchy.linkage(
            self.df.reindex(index=row_list).values,
            metric=self.metric,
            method=self.method,
            optimal_ordering=False
        )

        # Get the order of rows based on that list
        # and map those indices back to the input list
        return [row_list[i] for i in hierarchy.leaves_list(L)]

# If this is being run as a script
if __name__ == "__main__":

    ##################
    # SET UP LOGGING #
    ##################

    # Set the level of the logger to INFO
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [Order Genes] %(message)s'
    )
    logger = logging.getLogger(str(uuid4()))
    logger.setLevel(logging.INFO)

    # Write to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    # Get the user-supplied positional arguments
    input_fp = sys.argv[1]
    logger.info(f"Reading input from {input_fp}")
    output_fp = sys.argv[2]

    # Read in all of the alignment information
    df = pd.read_csv(
        input_fp
    # Pivot to wide format
    ).pivot_table(
        columns="genome",
        index="sseqid",
        values="pident"
    ).fillna(
        0
    # Convert to a presence / absence boolean
    ).applymap(
        lambda v: int(v > 0)
    )

    # Get the order of rows based on linkage clustering
    gene_order = OrderRows(
        df,
        logger=logger
    ).row_order

    logger.info(f"Writing output to {output_fp}")

    # Write out to a file
    with gzip.open(
        output_fp,
        "wt"
    ) as handle:
        handle.write("\n".join(gene_order))

    logger.info("DONE")