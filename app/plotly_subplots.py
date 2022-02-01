import pandas as pd
import plotly.graph_objects as go

class PlotlySubplots:
    """Class used to coordinate the creation of a Plotly figure with optional subplots."""

    def __init__(self, logger=None):

        # Store a dict of PlotlySubplot objects
        self.subplots = dict()

        # Attach a Plotly figure
        self.fig = go.Figure()

        # The index positions of each axis will be assigned
        # to an arbitrary string ID
        self.x_axes = dict()
        self.y_axes = dict()

        # Attach the logger, if any
        self.logger = logger

    def add(
        self,
        # ID for the subpanel
        id='subpanel',
        # Ordinal position on the horizontal axis
        x_index=0,
        # Ordinal position on the vertical axis
        y_index=0,
        # Relative width of the column
        width=1,
        # Relative height of the row
        height=1,
        # Optional padding around panel
        padding=0,
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

                # The padding must be the same
                msg = f"Values for `padding` must match for all subplots in column {x_index}"
                assert padding == subplot.padding, msg

                # The share_x must be the same
                msg = f"Values for `share_x` must match for all subplots in column {x_index}"
                assert share_x == subplot.share_x, msg

            # If the same row is being used
            if y_index == subplot.y_index:

                # The height must be the same
                msg = f"Values for `height` must match for all subplots in row {y_index}"
                assert height == subplot.height, msg

                # The padding must be the same
                msg = f"Values for `padding` must match for all subplots in row {y_index}"
                assert padding == subplot.padding, msg

                # The share_y must be the same
                msg = f"Values for `share_y` must match for all subplots in row {y_index}"
                assert share_y == subplot.share_y, msg

        # Add the subplot
        self.subplots[id] = PlotlySubplot(
            id=id,
            x_index=x_index,
            y_index=y_index,
            width=width,
            height=height,
            padding=padding,
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

    def get_axis_label(self, id, ax=None):
        """Return the x/y-axis label used for this subplot."""

        assert ax in ['x', 'y'], "ax must be x or y"

        # Iterate over the subplots
        for subplot in self.subplots.values():

            # If the ID is a match
            if subplot.id == id:

                # Return the appropriate axis label
                if ax == "x":
                    return self.x_axes[subplot.x_index]

                elif ax == "y":
                    return self.y_axes[subplot.y_index]

    def format_axis(self, id:str=None, params:dict=None, ax:str=None, anchor:str=None):
        """Apply formatting to an axis, identified by ordinal position."""

        assert ax in ['x', 'y'], "ax must be x or y"

        assert isinstance(params, dict), "Must provide params as a dict"

        # Get the name of the axis which was used for this subplot
        axis_name = self.get_axis(id, ax=ax)

        # While the axis name may be 'x2', the key used in the layout
        # would then be 'xaxis2'
        axis_name = self.format_axis_name(axis_name)

        # Apply the layout to the appropriate axis
        self.fig.update_layout(
            **{
                axis_name: params
            }
        )

        # If the user has requested that this axis be anchored
        if anchor is not None:

            # Get the name of the axis to set the anchor
            anchor_axis_name = self.get_axis(anchor, ax=ax)
            self.fig.update_layout(
                **{
                    axis_name: dict(
                        anchor=anchor_axis_name
                    )
                }
            )

    def format_axis_name(self, axis_name):
        """Convert y -> yaxis, y2 -> yaxis2, etc."""

        return f"{axis_name[0]}axis{axis_name[1:] if len(axis_name) > 1 else ''}"

    def anchor_xaxis(self, x_index, side=None):
        """Set the anchor for a particular x index on the top or bottom."""

        # `side` must be 'top' or 'bottom'
        assert side in ['top', 'bottom'], f"side must be top or bottom, not '{side}'"

        # Get the list of y axis indices for every subplot with this x axis
        y_indices = [
            subplot.y_index
            for subplot in self.subplots.values()
            if subplot.x_index == x_index
        ]

        # Make sure that there are subplots in this row
        assert len(y_indices) > 0, f"No subplots found with x index {x_index}"

        # If the side is 'top', get the highest index
        if side == "top":
            anchor_index = max(y_indices)
        # If the side is 'bottom', get the lowest index
        elif side == "bottom":
            anchor_index = min(y_indices)

        # Get the name of the two axes
        anchor_axis = self.y_axes[anchor_index]
        plot_axis = self.format_axis_name(self.x_axes[x_index])

        # Update the figure
        self.fig.update_layout(
            **{
                plot_axis: dict(
                    anchor=anchor_axis,
                    side=side
                )
            }
        )
        
    def anchor_yaxis(self, y_index, side=None):
        """Set the anchor for a particular y index on the right or left."""

        # `side` must be 'left' or 'right'
        assert side in ['left', 'right'], f"side must be left or right, not '{side}'"

        # Get the list of x axis indices for every subplot with this y axis
        x_indices = [
            subplot.x_index
            for subplot in self.subplots.values()
            if subplot.y_index == y_index
        ]

        # Make sure that there are subplots in this row
        assert len(x_indices) > 0, f"No subplots found with y index {y_index}"

        # If the side is 'right', get the highest index
        if side == "right":
            anchor_index = max(x_indices)
        # If the side is 'left', get the lowest index
        elif side == "left":
            anchor_index = min(x_indices)

        # Get the name of the two axes
        anchor_axis = self.x_axes[anchor_index]
        plot_axis = self.format_axis_name(self.y_axes[y_index])

        self.log(f"Anchoring {plot_axis} to {anchor_axis} - {side}")

        # Update the figure
        self.fig.update_layout(
            **{
                plot_axis: dict(
                    anchor=anchor_axis,
                    side=side
                )
            }
        )

    def set_domains(self, logger=None):
        """Update the domains for all axes to provide the appropriate height and width."""

        # Set the domain of x and y axes independently
        self.set_domains_axis("x", logger=logger)
        self.set_domains_axis("y", logger=logger)

    def set_domains_axis(self, ax, logger=None):
        """Update the domain for either x or y axes."""

        assert ax in ["x", "y"]
        
        # Get the relative size of each subplot
        span = pd.Series(
            {
                subplot.axis(ax): subplot.span(ax)
                for subplot in self.subplots.values()
            }
        ).sort_index()

        # Get the relative padding around each subplot
        padding = pd.Series(
            {
                subplot.axis(ax): subplot.padding
                for subplot in self.subplots.values()
            }
        ).sort_index()

        # If there is only one axis
        if span.shape[0] == 1:

            # No need to adjust anything
            return

        # Otherwise, there are multiple axes to coordinate

        # Iterate over each axis to apply the padding
        pointer = 0
        pos = []
        for i in span.index.values:

            pos.append(
                dict(
                    ix=i,
                    start=pointer + padding.loc[i],
                    end=pointer + span.loc[i] + padding.loc[i]
                )
            )

            pointer = pointer + span.loc[i] + (padding.loc[i] * 2)

        pos = pd.DataFrame(pos).set_index('ix')

        # Divide by the final value
        pos = pos / pos.end.max()

        # Iterate over each axis
        for axis_index, vals in pos.iterrows():

            # Get the name of the axis
            if ax == "x":
                axis_name = self.x_axes[axis_index]

            elif ax == "y":
                axis_name = self.y_axes[axis_index]

            # Expand to the long form of the axis name
            axis_name = self.format_axis_name(axis_name)

            domain_start = vals.start
            domain_end = vals.end

            self.log(f"Setting domain of {axis_name} to {domain_start} - {domain_end}")

            # Update the domain of the axis
            self.fig.update_layout(
                **{
                    axis_name: dict(
                        domain=(domain_start, domain_end)
                    )
                }
            )

    def log(self, msg):
        if self.logger is not None:
            self.logger.info(msg)
        else:
            print(msg)


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
        width=1,
        # Relative height of the row
        height=1,
        # Optional padding around panel
        padding=0,
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
        self.padding = padding
        self.share_x = share_x
        self.share_y = share_y

    def axis(self, ax):
        """Return either x_index or y_index, based on user input ('x' or 'y')."""

        assert ax in ['x', 'y']

        if ax == 'x':
            return self.x_index
        else:
            return self.y_index

    def span(self, ax):
        """Return either width or height, based on user input ('x' or 'y')."""

        assert ax in ['x', 'y']

        if ax == 'x':
            return self.width
        else:
            return self.height
    