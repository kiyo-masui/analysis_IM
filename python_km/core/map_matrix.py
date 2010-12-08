"""Data container for matrices of size npix x npix x nfreq."""

import scipy as sp

import base_data
import kiyopy.custom_exceptions as ce

class MapMatrix(base_data.BaseData) :
    """This class holds data ."""

    axes = ('long1', 'lat1', 'long2', 'lat2', 'freq')



