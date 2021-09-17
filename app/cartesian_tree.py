import pandas as pd

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

    def extension_x_coords(self):
        """Return a list of x-coordinates to extend the tips of the tree."""

        # Format is a list with each node and its parent, separated by NaN values
        return self._list_link_to_tips(col_name="x")

    def extension_y_coords(self):
        """Return a list of y-coordinates to extend the tips of the tree."""

        # Format is a list with each node and its parent, separated by NaN values
        return self._list_link_to_tips(col_name="y")

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

    def _list_link_to_tips(self, col_name="x"):
        """Internal method: make a list of coordinates to extend the tips of the tree."""

        # col_name may only be x or y
        assert col_name in ['x', 'y'], f"Not recognized: {col_name}"

        # Populate a list which will be output
        output_list = []

        # Iterate over each row
        for _, r in self.df.iterrows():

            # If this is not a leaf
            if not r['is_leaf']:

                # Skip it
                continue

            # If this is a tip
            else:

                # Add the item to the list
                output_list.append(r[col_name])

                # Moving along the x axis
                if col_name == "x":

                    # Extend the tip to the maximum x for the table
                    output_list.append(self.df[col_name].max())

                # Moving along the y axis
                elif col_name == "y":

                    # Extending the tip will keep the same y coordinate
                    output_list.append(r[col_name])

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


def get_node_positions(distance_df, y_offset=-0.5):
    """Construct a CartesianTree object from a set of distances."""

    # If there are more than two genomes/groups
    if distance_df.shape[0] > 2:

        # Format the distances as expected by skbio
        distances_dm = DistanceMatrix(
            distance_df.values, 
            list(map(str, distance_df.index.values))
        )

        # Make a neighbor-joining tree
        tree = nj(distances_dm)

        # Root at midpoint
        tree = tree.root_at_midpoint()

    # If there are only two genomes/groups
    elif distance_df.shape[0] == 2:

        # Get the distance betweeen the genomes/groups
        distance_between = distance_df.values[0, 1]

        # Make a simple tree linking the two
        tree = TreeNode(
            name='root',
            children=[
                TreeNode(
                    name=distance_df.index.values[0],
                    length=distance_between / 2.
                ),
                TreeNode(
                    name=distance_df.index.values[1],
                    length=distance_between / 2.
                )
            ]
        )

    # If there is only one genomes/groups
    elif distance_df.shape[0] == 1:

        # Make a simple tree with a single leaf
        tree = TreeNode(
            name='root',
            children=[
                TreeNode(
                    name=distance_df.index.values[0],
                    length=0
                )
            ]
        )

    # Assign x/y to create a DataFrame
    node_positions = CartesianTree(
        tree,
        y_offset=y_offset,
    )

    # Return that CartesianTree object
    return node_positions


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

    # Get the positions of the nodes in the tree
    node_positions = get_node_positions(
        subset_dists,
        y_offset=-0.5,
    )

    # Lastly, make a mapping from the filenames to the same filename without any underscores
    filename_mapping = {
        fn.replace("_", " "): fn
        for fn in genome_list
    }

    # Now make a new list, ordered by the node_positions, with the complete
    # filenames (including the underscores)
    node_positions.genome_order = list(map(
        lambda s: filename_mapping.get(s, s),
        node_positions.genome_order
    ))

    # Return the layout of the tree
    return node_positions