import plotly.graph_objects as go

class PlotlySubplots:
    """Class used to coordinate the creation of a Plotly figure with optional subplots."""

    def __init__(self):

        # Store a dict of PlotlySubplot objects
        self.subplots = dict()

        # Attach a Plotly figure
        self.fig = go.Figure()

        # The index positions of each axis will be assigned
        # to an arbitrary string ID
        self.x_axes = dict()
        self.y_axes = dict()

    def add(
        self,
        # ID for the subpanel
        id='subpanel',
        # Ordinal position on the horizontal axis
        x_index=0,
        # Ordinal position on the vertical axis
        y_index=0,
        # Relative width of the column
        width=0,
        # Relative height of the row
        height=0,
        # Share the x-coordinates with all other panels in the column
        share_x=True,
        # Share the y-coordinates with all other panels in the row
        share_y=True
    ):
        """Add a PlotlySubplot."""

        # If the same ID has been used
        assert id not in self.subplots, f"Subpanel id '{id}' is not unique"

        # Iterate through all other subplots
        for subplot in self.subplots.values():

            # If the same set of coordinates have been used
            msg = f"Coordinates x={x_index}, y={y_index} are occupied"
            assert x_index != subplot.x_index or y_index != subplot.y_index, msg

            # If the same column is being used
            if x_index == subplot.x_index:

                # The width must be the same
                msg = f"Values for `width` must match for all subplots in column {x_index}"
                assert width == subplot.width, msg

                # The share_x must be the same
                msg = f"Values for `share_x` must match for all subplots in column {x_index}"
                assert share_x == subplot.share_x, msg

            # If the same row is being used
            if y_index == subplot.y_index:

                # The height must be the same
                msg = f"Values for `height` must match for all subplots in row {y_index}"
                assert height == subplot.height, msg

                # The share_y must be the same
                msg = f"Values for `share_y` must match for all subplots in row {y_index}"
                assert share_y == subplot.share_y, msg

        # Add the subplot
        self.subplots[id] = PlotlySubplot(
            id=id,
            x_index=x_index,
            y_index=y_index,
            width=width,
            share_x=share_x,
            share_y=share_y
        )

        # Assign each of the index integers to axis positions
        # (if they have not yet already been assigned)
        self.assign_index(self.x_axes, x_index, "x")
        self.assign_index(self.y_axes, y_index, "y")

    def assign_index(self, n_axes, n_index, base):
        """Assign ordinal indices to axis ID strings."""

        # If the n_index has not yet been assigned to an (arbitrary) axis id
        if n_axes.get(n_index) is None:

            # If this is the first axis
            if len(n_axes) == 0:

                # The key is "x" or "y"
                axis_key = base

            # If this is not the first axis
            else:

                # The key is "{x,y}2", etc.
                axis_key = f"{base}{len(n_axes) + 1}"

            n_axes[n_index] = axis_key

    def get_axis(self, id, ax=None):
        """Return the axis label for the indicated subplot."""

        # Make sure that the id matches an existing subplot
        assert id in self.subplots, f"There is no subplot named {id}"

        assert ax in ['x', 'y'], "ax must be 'x' or 'y'"

        # Return the label assigned to the indicated axis
        if ax == 'x':
            return self.x_axes[self.subplots[id].x_index]
        elif ax == 'y':
            return self.y_axes[self.subplots[id].y_index]

    def plot(self, id=None, trace=None):
        """Add a trace to the indicated subpanel."""

        assert id is not None, "Must provide an `id`"
        assert trace is not None, "Must provide a `trace`"

        # Make sure that the id matches an existing subplot
        assert id in self.subplots, f"There is no subplot named {id}"

        # Assign the trace to the appropriate axis
        trace.xaxis = self.get_axis(id, ax="x")
        trace.yaxis = self.get_axis(id, ax="y")

        # Add the trace to the figure
        self.fig.add_trace(trace)

    def get_xaxis(self, id, ax=None):
        """Return the x/y-axis label used for this subplot."""

        assert ax in ['x', 'y'], "ax must be x or y"

        # Iterate over the subplots
        for subplot in self.subplots:

            # If the ID is a match
            if subplot.id == id:

                # Return the appropriate axis label
                if ax == "x":
                    return self.x_axes[subplot.x_index]

                elif ax == "y":
                    return self.y_axes[subplot.y_index]


class PlotlySubplot:
    """Class to coordinate the location and content of a single Plotly subplot."""

    def __init__(
        self,
        # ID for the subpanel
        id='subpanel',
        # Ordinal position on the horizontal axis
        x_index=0,
        # Ordinal position on the vertical axis
        y_index=0,
        # Relative width of the column
        width=0,
        # Relative height of the row
        height=0,
        # Share the x-coordinates with all other panels in the column
        share_x=True,
        # Share the y-coordinates with all other panels in the row
        share_y=True
    ):

        # Attach the attributes of the subplot
        self.id = id
        self.x_index = x_index
        self.y_index = y_index
        self.width = width
        self.height = height
        self.share_x = share_x
        self.share_y = share_y
    