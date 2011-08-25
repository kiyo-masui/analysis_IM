"""Dirty map making module."""

import scipy as sp
import numpy.ma as ma

import tools

class Pointing(object):
    """Class represents the pointing operator."""

    def __init__(self, Blocks, map, scheme='nearest'):
        pass

    def apply_to_time_stream(self, time_stream, map_out=None):
        pass

    def get_matrix(self):
        pass


class Noise_independant_channels(object):

    def __init__(self):
        pass

    def apply_inverse(self):
        pass

    def transform_map_space(self, Pointing, matrix_out=None):
        pass


class Noise(Noise_independant_channels):

    def __init__(self):
        pass


