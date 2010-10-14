"""Module contains the base classes for fits file readers of variouse formats
of data (time stream, maps, etc.).
"""

import scipy as sp
import numpy.ma as ma
import pyfits

import kiyopy.custom_exceptions as ce
import base_data as bd

# These globals are the cards for (our custom) history entries in a fits header
card_hist = 'DB-HIST'
card_detail = 'DB-DET'


class Reader(object) :
    """Base class for fits readers."""

    # Note to programmer: assignments of a subset of one array to another
    # are often by reference, not by value.  To be safe, any array you are
    # going to modify should be forced to assign by value with the
    # sp.array(an_array) function.

    # Overwrite the following attribute for classes inheriting from this one.
    field_and_axis = {}

    # Cards to 



