#!/usr/bin/env python3

from collections import defaultdict
import logging
from operator import index
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from time import time
from uuid import uuid4

class NestedLinkageSorting:
    """
    Sort one axis of a DataFrame using a nested linkage clustering approach.
    """

    def __init__(
        self,
        # Input data
        data_frame:pd.DataFrame,
        # Linkage clustering method
        method:str="centroid",
        # Distance metric used for 
        metric:str="braycurtis",
        # The number of clusters to generate at each level
        k:int=100,
        # Only perform clustering for datasets with sufficient rows
        min_cluster_threshold:int=100,
        # Function to use when aggregating members of groups
        agg_func=np.mean
    ):

        ##################
        # SET UP LOGGING #
        ##################

        # Set the level of the logger to INFO
        logFormatter = logging.Formatter(
            '%(asctime)s %(levelname)-8s [NestedLinkageSorting] %(message)s'
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
        self.data_frame = data_frame
        self.method = method
        self.metric = metric
        assert isinstance(k, int)
        self.k = k

        assert min_cluster_threshold > 0
        self.min_cluster_threshold = min_cluster_threshold

        self.agg_func = agg_func
        self.n = self.data_frame.shape[0]

        #################
        # ORDINATE DATA #
        #################
        self.ordinate()

    def ordinate(self):
        """Run PCA to generate new coordinates in the embedded distance metric space."""

        # Calculate the pairwise distances between genes
        t = time()
        self.dists = pd.DataFrame(
            squareform(
                pdist(
                    self.data_frame,
                    self.metric
                )
            ),
            index=self.data_frame.index,
            columns=self.data_frame.index,
        )
        elapsed = time() - t
        self.logger.info(f"Calculated pairwise distances in {elapsed:.2E} seconds")

        # Initialize the PCA object
        self.pca = PCA(
            n_components=self.dists.shape[0]
        )

        # Fit with PCA
        t = time()
        self.pca_df = pd.DataFrame(
            self.pca.fit_transform(
                self.dists
            ),
            index=self.dists.index
        )
        elapsed = time() - t
        self.logger.info(f"Ran PCA in {elapsed:.2E} seconds")

    def sort(self):
        """
        Sort the index of the `data_frame` with the following recursive approach.
        1. Split up the rows into k groups;
        2. Perform linkage clustering on the aggregated data for each group;
        3. Optimize the leaf ordering of groups to minimize differences between adjacent branches;
        4. Reorganize `data_frame` based on the order of the nested groups.
        """

        # Start the clock
        self.start_time = time()

        # Perform k-means clustering
        self.find_groups()

        # Perform linkage clustering on the aggregated observations from each group
        self.order_groups()

        # Concatenate the values from each group to reorder `data_frame`
        self.concatenate_groups()

        elapsed = time() - self.start_time
        self.logger.info(f"Sorted {self.n:,} rows in {elapsed:.2E} seconds")

    def find_groups(self):
        """Perform k-means clustering on all of the rows in `data_frame`."""

        # If there aren't at least `min_cluster_threshold` rows in `data_frame`
        if self.n <= self.min_cluster_threshold:

            self.logger.info(f"Table with {self.n:,} rows -- skipping further clustering")

            # Set up each row as its own group
            self.groups = [
                [i]
                for i in self.data_frame.index.values
            ]

            # The `group_data_frame` is the same as `data_frame`
            self.group_data_frame = self.data_frame

        # Otherwise, there are enough rows for clustering
        else:

            # Perform k-means clustering
            self.cluster()

    def cluster(self):
        """Perform k-means clustering."""

        # Set up the object
        model = KMeans(
            n_clusters=self.k
        )

        # Fit to the data and predict labels for each row
        labels = model.fit_predict(
            self.pca_df,
            sample_weight=self.pca.explained_variance_ratio_
        )

        # Make a list of lists with the members of each cluster
        self.groups = defaultdict(list)

        # Iterate over each assignment
        for i, label in enumerate(labels):

            # Add the label of the index to the group
            self.groups[label].append(
                self.data_frame.index.values[i]
            )

        # Make a DataFrame which groups those values
        self.group_data_frame = self.data_frame.assign(
            GROUP=labels
        ).groupby(
            "GROUP"
        ).apply(
            self.agg_func
        ).sort_index()

        # Transform from a dist of lists to a list of lists
        self.groups = [
            self.groups[i]
            for i in self.group_data_frame.index.values
        ]

    def order_groups(self):
        """Order the groups based on linkage clustering."""

        # If there are not at least 3 groups
        if len(self.groups) < 3:

            # No need to perform clustering
            return

        # Calculate pairwise distances
        dists = pdist(self.group_data_frame, self.metric)

        # Perform linkage clustering
        Z = hierarchy.linkage(
            dists,
            method=self.method
        )

        # Order the tree so that adjacent leaves are more similar
        Z_ordered = hierarchy.optimal_leaf_ordering(
            Z,
            dists
        )

        # Get the ordered list of leaves
        leaves_list = hierarchy.leaves_list(Z_ordered)

        # Reorder the group_data_frame
        self.group_data_frame = self.group_data_frame.iloc[leaves_list]

        # Reorder the group membership list
        self.groups = [
            self.groups[i]
            for i in leaves_list
        ]

    def concatenate_groups(self):
        """Make an expanded DataFrame with the sorted contents of each group."""

        # If there aren't at least `min_cluster_threshold` rows in `data_frame`
        if self.n <= self.min_cluster_threshold:

            # Make sure that each group is a single row
            for g in self.groups:
                assert len(g) == 1

            # Set the order of the input data using the order of the 'groups'
            self.data_frame = self.group_data_frame

        # Otherwise, make sure to sort the subgroups
        else:
            self.data_frame = pd.concat(
                [
                    self.sort_subgroup(group_members)
                    for i, group_members in enumerate(self.groups)
                ]
            )

    def sort_subgroup(self, group_members):
        """For a subset of the input data, sort and return the sorted DataFrame."""

        # Make a NestedLinkageSorting object for this group
        subgroup = NestedLinkageSorting(
            self.data_frame.reindex(index=group_members).copy(),
            method=self.method,
            metric=self.metric,
            k=self.k,
            min_cluster_threshold=self.min_cluster_threshold,
            agg_func=self.agg_func
        )

        # Sort it
        subgroup.sort()

        # Return the sorted values
        return subgroup.data_frame