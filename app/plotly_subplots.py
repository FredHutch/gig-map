import json
import pandas as pd
import plotly.graph_objects as go
    

class PlotlySubplots:
    """Class used to coordinate the creation of a Plotly figure with optional subplots."""

    def __init__(self, logger=None):

        # Store a DataFrame with all of the columns used to identify subplots
        self.subplots = pd.DataFrame(
            columns=[
                "id",
                "x_index",
                "y_index",
                "z_index",
                "x_shortname",
                "x_longname",
                "y_shortname",
                "y_longname",
                "x_span",
                "y_span",
                "x_padding",
                "y_padding"
            ]
        )

        # Attach a Plotly figure
        self.fig = go.Figure()

        # Attach the logger, if any
        self.logger = logger

        # Keep a list of any axes which should be made toggleable
        self.toggle_axes = dict()

    def add(
        self,
        # ID for the subpanel
        id='subpanel',
        # Ordinal position on the horizontal axis
        x_index=0,
        # Ordinal position on the vertical axis
        y_index=0,
        # Ordinal position on the z axis
        z_index=0,
        # Relative width of the column
        x_span=1,
        # Relative height of the row
        y_span=1,
        # Optional horizontal padding around panel
        x_padding=0,
        # Optional vertical padding around panel
        y_padding=0,
    ):
        """Add a subplot to the larger display."""

        # If the same ID has been used
        assert id not in self.subplots["id"].values, f"Subpanel id '{id}' is not unique"

        # Iterate through all other subplots
        for _, subplot in self.subplots.iterrows():

            # If the same set of coordinates have been used
            msg = f"Coordinates x={x_index}, y={y_index} are occupied"
            assert x_index != subplot.x_index or y_index != subplot.y_index, msg

            # If the same column is being used
            if x_index == subplot.x_index:

                # The x_span must be the same
                msg = f"Values for `x_span` must match for all subplots in column {x_index}"
                assert x_span == subplot.x_span, msg

                # The x_padding must be the same
                msg = f"Values for `x_padding` must match for all subplots in column {x_index}"
                assert x_padding == subplot.x_padding, msg

            # If the same row is being used
            if y_index == subplot.y_index:

                # The y_span must be the same
                msg = f"Values for `y_span` must match for all subplots in row {y_index}"
                assert y_span == subplot.y_span, msg

               # The y_padding must be the same
                msg = f"Values for `y_padding` must match for all subplots in column {x_index}"
                assert y_padding == subplot.y_padding, msg

        # Based on the x, y, and z indices, generate the short and long
        # axis names which are used by plotly to coordinate the figure layout
        x_shortname, x_longname = self.assign_index(index=x_index, ax="x", z_index=z_index)
        y_shortname, y_longname = self.assign_index(index=y_index, ax="y", z_index=z_index)

        # Format the new subplot
        new_subplot = dict(
            id=id,
            x_index=str(x_index),
            y_index=str(y_index),
            z_index=str(z_index),
            x_shortname=x_shortname,
            x_longname=x_longname,
            y_shortname=y_shortname,
            y_longname=y_longname,
            x_span=x_span,
            y_span=y_span,
            x_padding=x_padding,
            y_padding=y_padding
        )

        # Log it
        self.log(json.dumps(new_subplot))

        # Add the subplot
        self.subplots = pd.concat([
            self.subplots,
            pd.DataFrame([new_subplot]),
        ])

    def query_subplots(self, query_list):
        """Subset the subplots table using any list of queries."""

        # If a single query was provided
        if isinstance(query_list, str):

            # Turn the single query into a list with one member
            query_list = [query_list]

        # The input must otherwise be a list
        assert isinstance(query_list, list)

        # Start with the full table of subplots
        queried_subsets = self.subplots

        # Iterate over each query
        for query in query_list:

            # If there are any subplots remaining
            if queried_subsets.shape[0] > 0:

                # Filter it based on the query string
                try:
                    queried_subsets = queried_subsets.query(query)
                except Exception as e:
                    self.log(f"There was an error while applying the subplot query {query}")
                    self.log(queried_subsets)
                    raise e

        # Return the final queried table
        return queried_subsets

    def assign_index(self, index=None, ax=None, z_index=None):
        """Assign ordinal indices to axis ID strings."""

        # This function will return a tuple of shortname and longname, like
        # x2, xaxis2
        #  or 
        # y, yaxis

        # Either x_index or y_index was provided, not both
        assert index is not None and isinstance(index, int), "Must provide an integer index"
        assert ax is not None and ax in ["x", "y"], "ax must be x or y"

        # The z-axis must be provided
        assert z_index is not None, "Must provide z_index"

        # First check to see if this axis already exists
        overlap_df = self.query_subplots([
            f"z_index == '{z_index}'",
            f"{ax}_index == '{index}'"
        ])
        
        # If there are any axes already defined
        if overlap_df.shape[0] > 0:

            # They must all have the same name
            assert overlap_df[f"{ax}_shortname"].unique().shape[0] == 1, overlap_df
            assert overlap_df[f"{ax}_longname"].unique().shape[0] == 1, overlap_df

            # Return the predefined shortname and longname
            return overlap_df[f"{ax}_shortname"].values[0], overlap_df[f"{ax}_longname"].values[0]

        # If not, a new axis name needs to be assigned
        else:

            # Get the number of previous axes
            n_prev = self.subplots[f"{ax}_shortname"].unique().shape[0]

            # Format the new shortname and longname
            shortname = f"{ax}{'' if n_prev == 0 else n_prev + 1}"
            longname = f"{ax}axis{'' if n_prev == 0 else n_prev + 1}"

            return shortname, longname

    def get_attr(self, id, attr):
        """Return the axis label for the indicated subplot."""

        # Make sure that the id matches an existing subplot
        assert id in self.subplots["id"].values, f"There is no subplot named {id}"

        # Make sure that the requested attribute is valid
        assert attr in self.subplots.columns.values, f"No subplot attribute named '{attr}'"

        return self.subplots.set_index("id").loc[id, attr]

    def plot(self, id=None, trace=None):
        """Add a trace to the indicated subpanel."""

        assert id is not None, "Must provide an `id`"
        assert trace is not None, "Must provide a `trace`"

        # Assign the trace to the appropriate axis
        trace.xaxis = self.get_attr(id, "x_shortname")
        trace.yaxis = self.get_attr(id, "y_shortname")

        # Add the trace to the figure
        self.fig.add_trace(trace)

    def format_axis(self, id:str=None, params:dict=None, ax:str=None, log:bool=False):
        """Apply formatting to an axis, identified by ordinal position."""

        assert ax in ['x', 'y'], "ax must be x or y"

        assert isinstance(params, dict), "Must provide params as a dict"

        # Get the name of the axis which was used for this subplot
        axis_longname = self.get_attr(id, f"{ax}_longname")

        # If the log flag has been set
        if log:

            # Log the layout options which are being applied
            for k, v in params.items():
                self.log(f"Setting {axis_longname}: {k} = {v}")

        # Apply the layout to the appropriate axis
        self.fig.update_layout(
            **{
                axis_longname: params
            }
        )

    def get_index_sizes(self, ax, attr):
        """Get all of the x_span values sorted by x_index, for example."""
        index_cname = f"{ax}_index"
        assert index_cname in self.subplots.columns.values, f"Not a valid axis: {ax}"
        attr_cname = f"{ax}_{attr}"
        assert attr_cname in self.subplots.columns.values, f"Not a valid attribute: {attr}"

        # Return a vector, sorted by index
        return self.subplots.reindex(
            columns=[index_cname, attr_cname]
        ).drop_duplicates(
        ).set_index(
            index_cname
        )[
            attr_cname
        ].sort_index()

    def set_domains(self):
        """Update the domains for all axes to provide the appropriate height and width."""

        # Set the domains for each axis
        for ax in ["x", "y"]:
        
            # Get the relative size of each subplot, based on the index position
            span = self.get_index_sizes(ax, "span")

            # Get the relative padding around each subplot
            padding = self.get_index_sizes(ax, "padding")

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

            # Iterate over each index position
            for axis_index, vals in pos.iterrows():

                # Get the name of all axes which share this index
                for axis_name in self.query_subplots(
                    f"{ax}_index == '{axis_index}'"
                )[
                    f"{ax}_longname"
                ].unique():

                    self.log(f"Setting domain of {axis_name} to {vals.start} - {vals.end}")

                    # Update the domain of the axis
                    self.fig.update_layout(
                        **{
                            axis_name: dict(
                                domain=(vals.start, vals.end)
                            )
                        }
                    )

    def open_column(self, y_index, side="right"):
        """Return the first open position on the right or left of a row."""

        # Get the list of x axis indices for every subplot with this y axis
        x_indices = self.query_subplots(
            f"y_index == '{y_index}'"
        )["x_index"]

        msg = f"No subplots found in the row {y_index}"
        assert x_indices.shape[0] > 0, msg

        if side == "right":
            return x_indices.max() + 1
        else:
            assert side == "left", f"side can be 'right' or 'left', not '{side}'"
            return x_indices.min() - 1

    def open_row(self, x_index, side="bottom"):
        """Return the first open position on the top or bottom of a column."""

        # Get the list of y axis indices for every subplot with this x axis
        y_indices = self.query_subplots(
            f"x_index == '{x_index}'"
        )["y_index"]

        msg = f"No subplots found in the column {x_index}"
        assert y_indices.shape[0] > 0, msg

        if side == "top":
            return y_indices.max() + 1
        else:
            assert side == "bottom", f"side can be 'top' or 'bottom', not '{side}'"
            return y_indices.min() - 1

    def log(self, msg):
        if self.logger is not None:
            self.logger.info(msg)
        else:
            print(msg)

    def toggle_axis(self, axis, label):
        """Enable toggle functionality for the tick labels on a particular axis."""

        self.log(f"Adding toggle buttons for the {axis} axis showing {label} values")

        # Add the axis and label to the dict of axes to toggle
        self.toggle_axes[axis] = label

    def untoggle_axis(self, axis):
        """Disable toggle functionality for the tick labels on a particular axis."""

        # Remove the axis from the dict of axes to toggle
        del self.toggle_axes[axis]

    def add_interactivity(self):
        """Add any interactive buttons or sliders which have been defined."""

        self.fig.update_layout(
            updatemenus=[
                dict(
                    type="buttons",
                    active=-1,
                    showactive=True,
                    pad=dict(t=20),
                    buttons=[
                        i
                        for axis, label in self.toggle_axes.items()
                        for i in [
                            dict(
                                label=f"Hide {label} labels",
                                method="relayout",
                                args=[{f"{axis}.showticklabels": False}]
                            ),
                            dict(
                                label=f"Show {label} labels",
                                method="relayout",
                                args=[{f"{axis}.showticklabels": True}]
                            )
                        ]
                    ]
                )
            ]
        )

