import types

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
        args:list=None,
        read_f:function=None,
        plot_f:function=None,
    ):
        """
        Instantiate a figure element.
        args: A list of FigureArgument objects;
        read_f: A function which returns a data object,
                based on the values provided for those arguments;
        plot_f: A function which adds to a PlotlySubplot object,
                based on the data object and arguments;
        """
        
        assert isinstance(args, list), "args must be provided as a list"
        for arg in args:
            assert isinstance(arg, FigureArgument), "args must be a list of FigureArgument objects"
        assert isinstance(read_f, types.FunctionType), "read_f must be a function"
        assert isinstance(plot_f, types.FunctionType), "plot_f must be a function"

        self.args = args
        self.read_f = read_f
        self.plot_f = plot_f
    

class FigureBuilder:
    """Class to coordinate the construction of a figure from flexible data inputs."""

    def __init__(self, **kwargs):

        # Store all of the elements of the figure associated with a keyword
        self.elements = dict()

        # Add all of the elements provided for this figure
        for key, element in kwargs.items():
            self.add_element(key, element)

    def add_element(self, key:str, element:FigureElement):
        """Add a FigureElement."""

        assert isinstance(key, str), "Figure elements must be identified with a string"
        assert isinstance(element, FigureElement), "Figure elements must be FigureElement objects"

        self.elements[key] = element
