#!/usr/bin/env python3

import logging
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial import distance
from time import time
from uuid import uuid4


class Timer:
    """Keep track of an interval of time."""

    def __init__(self, logger=None):

        # Start the timer
        self.start_time = time()

        # Attach the logger, if any
        self.logger = logger

    def report(self, msg):

        # Track the amount of time since the timer was started
        elapsed_time = time() - self.start_time

        # If a logger was provided
        if self.logger is not None:

            # Print the message to the logger
            self.logger.info(msg.format(elapsed_time))

        # If no logger was provided
        else:

            # Print the message
            print(msg.format(elapsed_time))


class GreedyLinkageSorting:
    """
    Sort the rows of a DataFrame using a greedy linkage clustering approach.
    """

    def __init__(
        self,
        # Input data
        df:pd.DataFrame,
        # Distance metric used for comparing items
        metric='braycurtis',
        # Method used for linkage clustering
        method='average',
        # Maximum distance value used to group initial clusters
        threshold=0.25,
        # Number of seconds between status messages
        logging_interval:int=2,
        # Print logging message
        verbose:bool=True,
    ):
        """Save the input data to the object and order the rows."""
        
        ##################
        # SET UP LOGGING #
        ##################

        # Set the level of the logger to INFO
        logFormatter = logging.Formatter(
            '%(asctime)s %(levelname)-8s [GreedyLinkageSorting] %(message)s'
        )
        self.logger = logging.getLogger(str(uuid4()))
        self.logger.setLevel(logging.INFO)

        # Write to STDOUT
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        self.logger.addHandler(consoleHandler)

        ###############
        # ATTACH DATA #
        ###############

        # Save the values
        self.df = df

        # Save the parameters
        self.metric = metric
        self.logger.info(f"metric: {metric}")
        self.method = method
        self.logger.info(f"method: {method}")
        self.threshold = threshold
        self.logger.info(f"threshold: {threshold}")
        self.logging_interval = logging_interval
        self.verbose = verbose
        
        #############
        # SORT DATA #
        #############

        # Get the order of the rows
        self.row_order = self.calc_row_order()
        
    def calc_row_order(self):
        """Calculate the optimal order of rows."""
        
        msg = f"Sorting {self.df.shape[0]:,} features using data from {self.df.shape[1]:,} observations"
        self.logger.info(msg)

        # Set up the clusters that each row will be assigned to
        self.clusters = []
        
        # Track the amount of time needed to build clusters
        self.start_time = time()

        # Keep track of how long it has been since the last status
        # message was printed (if verbose=True).
        status_time = time()

        # Keep track of the time it took to join each set of pairs
        time_to_join = []

        # Keep track of the total number of rows which have been added
        n_added = 0

        # Iterate over each row in the table
        for row_name, row_vals in self.df.iterrows():
            
            # Keep track of how long it takes to join each new pair
            t = time()

            # Add the row to the first cluster which it matches
            self.add_row_to_cluster([row_name], row_vals.values)

            # Increment the counter
            n_added += 1

            # Keep track of how long it takes to join each new pair
            time_to_join.append(time() - t)
            
            # If logging was turned on
            if self.verbose:

                # And more than logging_interval has passed
                if (time() - status_time) > self.logging_interval:

                    # Calculate the mean amount of time it took to join each pair
                    mean_time_to_join = np.mean(time_to_join)

                    # Log the message
                    self.logger.info(f"Rows added: {len(time_to_join):,} ({n_added:,} total) - Iteration time (s): {mean_time_to_join:.2E}")
                    
                    # Reset the time
                    status_time = time()
                    
                    # Reset the list
                    time_to_join = []

        # Calculate the amount of time it took to build each cluster
        elapsed = time() - self.start_time
        self.logger.info(f"Added {n_added:,} rows to {len(self.clusters):,} clusters in {round(elapsed, 2):,} seconds")

        # Order the clusters by linkage clustering
        self.order_clusters()

        # Build the row order as a list
        row_order = []

        # Keep track of the amount of time needed to sort the rows within each cluster
        timer = Timer(logger=self.logger)
        
        # Iterate over the clusters, which have been ordered by linkage clustering
        for cluster in self.clusters:
            
            # Add the rows in that cluster to the overall list
            row_order.extend(self.order_rows_in_cluster(cluster))

        # Report the time elapsed
        timer.report("Rows sorted within each cluster - {:,} seconds")
            
        # Return the overall list of rows
        return row_order
        
    def add_row_to_cluster(self, row_name_list, row_vals):
        """Add a single row to the cluster that it is most similar to."""

        # If there are no existing clusters
        if len(self.clusters) == 0:

            # Make a new cluster containing this row
            self.clusters.append(
                RowCluster(
                    row_names=row_name_list,
                    row_values=row_vals
                )
            )

        # If there are already >0 clusters created
        else:

            # Calculate the distance to each cluster
            cluster_dist = self.distance_to_clusters(row_vals)

            # If any cluster is below the threshold
            if cluster_dist.min() <= self.threshold:

                # Reference the best cluster
                best_cluster_ix = cluster_dist.idxmin()
                best_cluster = self.clusters[best_cluster_ix]

                # Add the row to the cluster
                best_cluster.add_row(row_name_list, row_vals)

                # Since the cluster has been modified, we now need to
                # check and see if it should be merged with another cluster

                # Remove this cluster from the list of clusters
                self.clusters.pop(best_cluster_ix)

                # Recursively add this cluster back to the batch
                self.add_row_to_cluster(
                    best_cluster.row_names,
                    best_cluster.avg_values
                )

            # Otherwise, if no cluster is below the threshold
            else:
                
                # Make a new cluster containing this row
                self.clusters.append(
                    RowCluster(
                        row_names=row_name_list,
                        row_values=row_vals
                    )
                )

    def distance_to_clusters(self, row_vals):
        """Return a Series with the distance of a vector to all clusters."""

        cluster_dist = pd.Series(
            distance.cdist(
                np.array([row_vals]),
                np.array([
                    cluster.avg_values
                    for cluster in self.clusters
                ]),
                metric=self.metric
            )[0]
        )

        # The number of distances must equal the number of clusters
        assert cluster_dist.shape[0] == len(self.clusters)

        return cluster_dist

    def order_clusters(self):
        """Sort the clusters by linkage clustering."""

        # If there are fewer than 3 clusters
        if len(self.clusters) < 3:

            # There is no need to do any sorting
            self.logger.info(f"There is no need to sort {len(self.clusters)} clusters")

            # Take no further action
            return

        # Otherwise, if there are >=3 clusters

        # Keep track of the time needed to sort the clusters
        timer = Timer(logger=self.logger)

        # Make a matrix with the average values across all clusters
        cluster_mat = np.array(
            [
                cluster.avg_values
                for cluster in self.clusters
            ]
        )

        # Make sure that there are no NaN values
        assert not np.isnan(cluster_mat).any(), f"Error: matrix contains NaN values ({cluster_mat.shape})"

        # Get the list of sorted index positions using exhaustive linkage clustering
        ix_list = self.sort_by_linkage_clustering(cluster_mat)

        # Reorder the self.clusters object
        self.clusters = [self.clusters[ix] for ix in ix_list]

        # Report the time elapsed
        timer.report("Ordered all clusters in {:,} seconds")

    def order_rows_in_cluster(self, cluster):
        """Return the sorted list of row names for a cluster."""

        # If there are < 3 rows in this cluster
        if cluster.size < 3:

            # Sorting is not required
            return cluster.row_names

        # First, make a matrix with the values for each of the rows in this cluster
        row_mat = self.df.reindex(cluster.row_names).values

        # Get the sorted list of index positions
        ix_list = self.sort_by_linkage_clustering(row_mat)

        # Return the sorted list of row names
        return [
            cluster.row_names[i]
            for i in ix_list
        ]            

    def sort_by_linkage_clustering(self, mat):
        """Return the sorted list of index positions for a numpy array."""

        # Calculate the pairwise distances
        dists = distance.pdist(
            mat,
            metric=self.metric
        )

        # If there are any NaNs
        if np.isnan(dists).any():

            # Replace them with 0
            np.nan_to_num(
                dists,
                copy=False,
                nan=0.0,
                posinf=1.0,
                neginf=0.0
            )

        try:
            return hierarchy.leaves_list(
                hierarchy.linkage(
                    dists,
                    method=self.method,
                    optimal_ordering=True
                )
            )

        except Exception as e:
            self.logger.info(f"Error clustering matrix ({mat.shape})")
            self.logger.info(f"Distances ({dists})")
            self.logger.info(mat)
            raise e
    
            
class RowCluster:
    """Collection of rows which have been grouped together based on similarity."""
    
    def __init__(
        self,
        row_names:list=[],
        row_values:np.array=np.array([])
    ):
        """The RowCluster must be initialized with a group of rows."""
        
        assert len(row_names) > 0, "RowCluster must contain > 0 rows"
        assert row_values.shape[0] > 0, "RowCluster must contain > 0 values"
        
        # Save the row names
        self.row_names = row_names

        # Save the number of rows        
        self.size = len(row_names)
        
        # Save the values
        self.avg_values = row_values

        # Make sure that there are no NaN values
        assert not np.isnan(self.avg_values).any(), f"Error: matrix contains NaN values"

    def add_row(self, row_name_list, row_vals):
        """Add a row to this cluster."""
        
        # Make sure that the input has the right number of columns
        assert row_vals.shape[0] == self.avg_values.shape[0], (row_vals.shape[0], self.avg_values.shape[0])
        
        # Compute a new weighted average of the values
        self.avg_values = np.array(
            [
                self.avg_values * self.size,
                row_vals * len(row_name_list)
            ]
        ).sum(axis=0) / (self.size + len(row_name_list))
        
        # Make sure that there are no NaN values
        assert not np.isnan(self.avg_values).any(), f"Error: matrix contains NaN values"

        # Add the row name to the list
        self.row_names.extend(row_name_list)
        
        # Increment the size counter
        self.size += len(row_name_list)
