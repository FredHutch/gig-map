#!/usr/bin/env python3

import argparse
from collections import defaultdict
import logging
import pandas as pd
from plotly_subplots import PlotlySubplots
import sys
from types import FunctionType
from uuid import uuid4


class FigureArgument:
    """Define an argument used to drive the creation of a figure."""

    def __init__(
        self,
        key:str=None,
        description:str=None,
        type=str,
        default=None,
    ):

        assert isinstance(key, str), "key must be a string"
        assert isinstance(description, str), "description must be a string"

        assert not " " in key, "key cannot contain spaces"

        self.key = key
        self.description = description
        self.type = type
        self.default = default


class FigureElement:
    """Define a component of the overall FigureBuilder."""

    def __init__(
        self,
        id:str=None,
        args:list=None,
        read_f:FunctionType=None,
        plot_f:FunctionType=None,
    ):
        """
        Instantiate a figure element.
        id:     An identifier for the element;
        args:   A list of FigureArgument objects;
        read_f: A function which returns a data object,
                based on the values provided for those arguments;
        plot_f: A function which adds to a PlotlySubplot object,
                based on the data object and arguments;
        """
        
        assert isinstance(id, str), "id must be provided as a str"
        assert isinstance(args, list), "args must be provided as a list"
        for arg in args:
            assert isinstance(arg, FigureArgument), "args must be a list of FigureArgument objects"

        self.id = id
        self.args = args

        # If read_f was provided
        if read_f is not None:
            assert isinstance(read_f, FunctionType), "read_f must be a function"
            self.read_f = read_f

        # If plot_f was provided
        if plot_f is not None:
            assert isinstance(plot_f, FunctionType), "plot_f must be a function"
            self.plot_f = plot_f

        # By default, the element is enabled
        self.enabled = True

    def disable(self):
        """Prevent this element from being plotted."""

        self.enabled = False


class FigureAxis:
    """Configure an axis."""

    def __init__(self, id):

        # By default, no values are provided
        self.id = id
        self.axis = None
        self.exists = False
        self.is_fixed = False

    def set(self, axis):
        """Set the order of an axis with a pd.Series object."""

        msg = "Axis must be set by a pandas Series object"
        assert isinstance(axis, pd.Series), msg

        # Attach the Series
        self.axis = axis

        # Mark that the axis has been instantiated
        self.exists = True

    def member_set(self):
        """Return a set of the elements on the axis."""

        assert self.exists, f"Axis ({self.id}) has not yet been populated"

        return set(self.axis.index.values)

    def set_order(self, ordered_list):
        """Update the order of the axis."""
    
        assert self.exists, f"Axis ({self.id}) has not yet been populated"

        # Make sure that every element of the list is in the axis
        member_set = self.member_set()

        for item in ordered_list:

            assert item in member_set, f"Item ('{item}') not found in {self.id} axis"

        # Reorder the Series
        self.axis = self.axis.reindex(index=ordered_list)

        self.is_fixed = True

    def length(self):

        assert self.exists, f"Axis ({self.id}) has not yet been populated"

        return self.axis.shape[0]

    def labels(self):

        assert self.exists, f"Axis ({self.id}) has not yet been populated"

        return self.axis.values

    def order(self):

        assert self.exists, f"Axis ({self.id}) has not yet been populated"

        return self.axis.index.values

    def label_dict(self):

        assert self.exists, f"Axis ({self.id}) has not yet been populated"

        return self.axis.to_dict()


class FigureBuilder:
    """
    Class to coordinate the construction of a figure from flexible data inputs.

    The order of operations is:
        - read_data: 'global' first, and then each element in turn
        - make_plots: Elements first, and then 'global' last
    """

    def __init__(
        self,
        logging_id:str="figure-builder",
        description:str="Figure description",
        args:list=[],
        read_f:FunctionType=None,
        plot_f:FunctionType=None,
        elements:list=None
    ):

        # Store the logging_id associated with this figure builder
        self.logging_id = logging_id

        # Store the description used for the argument parser
        self.description = description

        msg = "Global arguments (args=) must be a list"
        assert isinstance(args, list), msg
        self.args = args

        msg = "Each argument must be a FigureArgument"
        for arg in self.args:
            assert isinstance(arg, FigureArgument), msg

        # If read_f was provided
        if read_f is not None:
            msg = "read_f must be a function"
            assert isinstance(read_f, FunctionType)
            self.read_f = read_f

        # If plot_f was provided
        if plot_f is not None:
            msg = "plot_f must be a function"
            assert isinstance(plot_f, FunctionType)
            self.plot_f = plot_f

        # Validate and store all of the `elements`
        self.elements = self.validate_elements(elements)

        # Create a self.logger object
        self.create_logger()

        # Store all axes as a dict
        self.axes = dict()

    def validate_elements(self, elements):
        """Validate that the input `elements` are unique FigureElements."""

        assert isinstance(elements, list), "`elements` must be a list"
        
        # All of the element `id` string must be unique
        all_ids = set()

        # All of the items in `elements` must be FigureElement
        msg = "Figure elements must be FigureElement objects"
        for element in elements:
            assert isinstance(element, FigureElement), msg

            assert element.id not in all_ids, f"Element id '{element.id}' is not unique"
            all_ids.add(element.id)

        return elements

    def create_logger(self):
        """Create a logging instance."""

        # Set the level of the logger to INFO
        logFormatter = logging.Formatter(
            f'%(asctime)s %(levelname)-8s [{self.logging_id}] %(message)s'
        )
        self.logger = logging.getLogger(str(uuid4()))
        self.logger.setLevel(logging.INFO)

        # Write to STDOUT
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        self.logger.addHandler(consoleHandler)

        # Print the description of the figure builder
        self.log(self.description)

        # Print the version of Python being used
        self.log(f"Python version: {sys.version}")
        self.log(sys.version_info)

        # Print the version of Pandas being used
        self.log(f"Pandas version: {pd.__version__}")

    def log(self, msg):
        """Emit a message in the logging stream."""
        
        # By default, all logging message will be INFO
        self.logger.info(msg)

    def parse_args(self):
        """Parse arguments from the command line."""

        # Instantiate the argument parser
        self.parser = argparse.ArgumentParser(
            description=self.description
        )

        # Iterate over each of the global-level arguments
        for arg in self.args:

            # Add this argument to the parser
            self.parser.add_argument(
                f"--{arg.key}",
                type=arg.type,
                default=arg.default,
                help=arg.description
            )

        # Iterate over each of the elements of the figure
        for element in self.elements:

            # Create an argument group
            arg_group = self.parser.add_argument_group(
                element.id,
            )

            # Iterate over each of the arguments for that element
            for arg in element.args:

                # Add this argument to the parser
                arg_group.add_argument(
                    f"--{element.id}-{arg.key}",
                    type=arg.type,
                    default=arg.default,
                    help=arg.description
                )

        # All arguments from the command line will be stored
        # in a nested dict, where the first key is either
        # 'global' or an `element_id`
        self.params = defaultdict(dict)

        # Parse the arguments, iterating over each one
        for kw, val in self.parser.parse_args().__dict__.items():

            self.add_param(kw, val)

        # Log the params
        self.log("Parameters:")
        for param_group, params in self.params.items():
            self.log(f"Group: {param_group}")
            for kw, val in params.items():
                self.log(f"       {kw}: {val}")
        self.log("End of Parameters")

    def add_param(self, kw, val):
        """Add a single argument to either the 'global' or `element_id` namespace."""

        # First iterate over the global arguments
        for arg in self.args:

            # If the keyword matches
            if arg.key.replace("-", "_") == kw:

                # Then add this value to the 'global' namespace
                self.params['global'][kw] = val
                return

        # Next, iterate over each of the elements of the figure
        for element in self.elements:

            # Iterate over each of the arguments for this element
            for arg in element.args:

                # If the keyword matches
                if f"{element.id}_{arg.key.replace('-', '_')}" == kw:

                    # Then add this value to the `element.id` namespace
                    self.params[
                        element.id
                    ][
                        arg.key.replace('-', '_')
                    ] = val
                    return
        
        # At this point, no match was found
        msg = f"Could not find a match for argument: {kw}"
        raise Exception(msg)

    def read_data(self):
        """Read in data using the `read_f` functions provided."""

        # First, read in the data for the global namespace
        self.log("Reading in data from the global namespace")
        self.read_f(
            self,
            **self.params['global']
        )

        # Next, iterate over each figure element
        for element in self.elements:

            # Read in the data for each of those elements
            self.log(f"Reading in data for element: {element.id}")
            element.read_f(
                self,
                **self.params[element.id]
            )

    def make_plots(self):
        """Make a plot using the arguments, data, and functions defined."""

        # Instantiate a subplot object
        self.subplots = PlotlySubplots(logger=self.logger)

        # For each element
        for element in self.elements:

            # Only make the plot if the element is enabled
            if element.enabled:

                # Call the plotting function
                element.plot_f(self)

        # Apply the global plotting function
        self.plot_f(self)

        # Arrange the axis domains
        self.subplots.set_domains()

    def axis(self, axis_name):
        """Return a FigureAxis object corresponding to axis_name."""

        if axis_name not in self.axes:

            self.axes[axis_name] = FigureAxis(axis_name)

        return self.axes[axis_name]

    def write_html(self, fp):
        """Write the figure to HTML"""

        self.log(f"Writing HTML to {fp}")

        self.subplots.fig.write_html(fp)