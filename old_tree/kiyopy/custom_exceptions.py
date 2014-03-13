"""Custom Exceptions.

These mostly don't do anything special, but are defined such that the
exceptions I raise don't conflict with generic exceptions raised by built-ins.
"""

class DataError(Exception) :
    """Exception to raise if data of some sort is invalid or does not have
    expected properties."""
    pass

class NextIteration(Exception) :
    """Exception raised to skip iterations in a nested loop in a controled way.
    """
    pass

class FileParameterTypeError(TypeError) :
    """Exception to raise if a parameter read from file should be a certain
    type and is not."""
    pass

class ParameterFileError(Exception) :
    """Exception to raise if reading a parameter file fails."""
    pass
