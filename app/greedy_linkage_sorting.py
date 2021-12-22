#!/usr/bin/env python3

import logging
import numpy as np
import pandas as pd
from random import choice
from scipy.spatial.distance import cdist, squareform
from time import time
from uuid import uuid4

class GreedyLinkageSorting:
    """
    Sort one axis of a DataFrame using a greedy linkage clustering approach.
    """

    def __init__(
        self,
        # Input data
        data_frame:pd.DataFrame,
        # Axis to sort
        axis:int=0,
        # Distance metric used for comparing items
        metric:str="euclidean",
        # Print logging message
        verbose:bool=True,
        # Number of seconds between status messages
        logging_interval:int=2,
    ):

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

        # Attach all provided data objects
        assert isinstance(data_frame, pd.DataFrame)
        if axis == 0:
            self.values = data_frame.values
        else:
            assert axis == 1, "axis must be 0 or 1"
            self.values = data_frame.T.values

        # Keep track of the number of elements which were input
        self.input_size = self.values.shape[0]

        self.metric = metric

        # After sorting, this object will be populated with a list of
        # sorted index positions
        self.row_order = None

        # Set up a cache for the best hits
        self.best_hit_cache = dict()

        self.verbose = verbose
        self.logging_interval = logging_interval

        #############
        # SORT DATA #
        #############
        self.sort()

    def sort(self):
        """Sort the rows of self.data_frame, populating row_order."""

        # Start the clock
        self.start_time = time()

        # Initialize the clusters with each row in its own cluster
        # self.clusters is a list where each item is a nested list
        # containing the ordered index positions from the original
        # table which are contained within it.
        # The average value of each feature for this cluster can be
        # found in self.values. Any changes to either must be kept
        # in sync with the other.
        self.clusters = [ [i] for i in range(self.values.shape[0]) ]

        self.logger.info(f"Sorting {self.values.shape[0]:,} features using data from {self.values.shape[1]:,} observations")
        
        # Keep track of how long it has been since the last status
        # message was printed (if verbose=True).
        status_time = time()

        # Keep track of the time it took to join each set of pairs
        time_to_join = []

        # Keep joining groups until there is only one cluster left
        while len(self.clusters) > 1:

            # Keep track of how long it takes to join each new pair
            t = time()

            # Join the first pair of clusters we can find which are mutual best hits
            self.join_mutual_best_hit()
            
            # Keep track of how long it takes to join each new pair
            time_to_join.append(time() - t)

            # If logging was turned on
            if self.verbose:

                # And more than logging_interval has passed
                if (time() - status_time) > self.logging_interval:

                    # Calculate the mean amount of time it took to join each pair
                    mean_time_to_join = np.mean(time_to_join)
                    # Log the message
                    self.logger.info(f"Clusters remaining: {len(self.clusters) - 1:,} - Iteration time (s): {mean_time_to_join:.2E}")
                    # Reset the time
                    status_time = time()
                    # Reset the list
                    time_to_join = []

        elapsed = time() - self.start_time
        self.logger.info(f"Sorted all rows in {round(elapsed, 2):,} seconds")

        # Save the list of ordered index positions
        self.row_order = self.clusters[0]

        # Make sure that the final ordered list contains all of the inputs
        assert len(self.row_order) == self.input_size, (len(self.row_order), self.input_size)

    def find_size(self, ix):
        """Return the number of items in a particular cluster."""
        return len(self.clusters[ix])

    def join_mutual_best_hit(self):
        """Find a pair of clusters which are mutual best hits and join them."""

        # Find a mutual best hit from self.values
        seed_ix, match_ix = self.find_mutual_best_hit()

        # Get the size of each cluster
        seed_ix_n = self.find_size(seed_ix)
        match_ix_n = self.find_size(match_ix)

        # Calculate the weights of each cluster based on their weights
        seed_ix_weight = seed_ix_n / (seed_ix_n + match_ix_n)
        match_ix_weight = match_ix_n / (seed_ix_n + match_ix_n)

        # Calculate the values which are the weighted average of the two clusters
        combined_row = np.mean(
            self.values[[seed_ix, match_ix]].T * np.array([seed_ix_weight, match_ix_weight]),
            axis=1
        )

        # Calculate the average value for each of the two clusters
        seed_ix_avg_value = np.mean(self.values[seed_ix])
        match_ix_avg_value = np.mean(self.values[match_ix])

        # Remove the two clusters and add the new cluster at the end
        # From the values
        self.values = np.insert(
            np.delete(
                self.values,
                [seed_ix, match_ix],
                axis=0
            ),
            self.values.shape[0] - 2,
            [combined_row],
            axis=0
        )

        # First add to the list of clusters

        # If the seed cluster has a higher average value,
        if seed_ix_avg_value > match_ix_avg_value:
            # place it first in the new cluster
            self.clusters.append(self.clusters[seed_ix] + self.clusters[match_ix])
        # Otherwise
        else:
            # Place the other cluster first
            self.clusters.append(self.clusters[match_ix] + self.clusters[seed_ix])

        # And then remove the clusters which were merged
        # Starting with the higher index value
        if seed_ix > match_ix:
            del self.clusters[seed_ix]
            del self.clusters[match_ix]
        else:
            del self.clusters[match_ix]
            del self.clusters[seed_ix]

        # Clear the cache of best hits
        self.best_hit_cache = dict()

    def check_cluster_size(self, checkpoint):
        """Make sure that the total number of clusters matches the input."""
        n_tot = sum(map(len, self.clusters))

        assert n_tot == self.input_size, f"{checkpoint}: Found {n_tot:,} clusters, expected {self.input_size:,}"

    def find_mutual_best_hit(self, seed_ix=None):
        """
        Find a pair of rows in self.values which are mutual best hits.
        Return the tuple of the index positions of each.
        """

        # If no seed was specified
        if seed_ix is None:

            # Start by picking a random row
            seed_ix = choice(range(self.values.shape[0]))

        # Find the list of most similar rows
        for match_ix in self.find_best_hit(seed_ix):

            # If seed_ix is also the closest match to match_ix
            if seed_ix in self.find_best_hit(match_ix):

                # Then return the pair
                return seed_ix, match_ix

        # If no match was found, start the search again by daisy-chaining from the best hit
        return self.find_mutual_best_hit(seed_ix=match_ix)

    def find_best_hit(self, ix):
        """Return the list of index positions for rows which are most similar to this one."""

        # If there is a cached value
        if self.best_hit_cache.get(ix) is not None:

            # Return it
            return self.best_hit_cache.get(ix)

        # Calculate the distnaces of all rows to this one
        dists = cdist(
            self.values[[ix]],
            self.values,
            metric=self.metric
        )[0]

        # Find the lowest distance to any other row
        lowest_value = np.min([d for other_ix, d in enumerate(dists) if other_ix != ix])

        # Find the index positions of the rows which are the closest match
        best_hits = [other_ix for other_ix, d in enumerate(dists) if other_ix != ix and d == lowest_value]

        assert len(best_hits) > 0

        # Cache the hit
        self.best_hit_cache[ix] = best_hits

        return best_hits
