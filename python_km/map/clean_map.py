"""Converter from a dirty map to a clean map."""

import scipy as sp
import scipy.linalg as linalg

from core import algebra
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce


params_init = {'input_root' : './',
               'dirty_map_files' : ('dirty_map_I.npy',),
               'noise_inverse_files' : ('noise_inv_I.npy',),
               'output_root' : './'
               }
prefix = 'cm_'

class CleanMapMaker(object) :
    """Converts a Dirty map to a clean map."""

    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix=prefix)
        self.feedback = feedback

    def execute(self, nprocesses=1) :
        """Worker funciton."""
        params = self.params
        # Make parent directory and write parameter file.
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix='mm_')
        # Loop over files to process.
        for ii in range(len(params['dirty_map_files'])) :
            dmap_fname = params['dirty_map_files'][ii]
            noise_fname = params['noise_inverse_files'][ii]
            # Load the dirty map and the noise matrix.
            dirty_map = algebra.load(params['input_root'] + dmap_fname)
            dirty_map = algebra.make_vect(dirty_map)
            if dirty_map.axes != ('freq', 'ra', 'dec') :
                raise ce.DataError("Expeced dirty map to have axes "
                                   "('freq', 'ra', 'dec'), but it has axes: "
                                   + str(dirty_map.axes))
            noise_inv = algebra.open_memmap(params['input_root'] + noise_fname,
                                           'r')
            noise_inv = algebra.make_mat(noise_inv)
            # Two cases for the noise.  If its the same shape as the map then
            # the noise is diagonal.  Otherwise, it should be block diagonal in
            # frequency.
            if noise_inv.ndim == 3 :
                if noise_inv.axes != ('freq', 'ra', 'dec') :
                    raise ce.DataError("Expeced noise matrix to have axes "
                                       "('freq', 'ra', 'dec'), but it has: "
                                       + str(dirty_map.axes))
                # Noise inverse can fit in memory, so copy it.
                noise_inv_memory = sp.array(noise_inv, copy=True)
                # Find the non-singular (covered) pixels.
                max_information = noise_inv_memory.max()
                good_data = noise_inv_memory < 1.0e-5*max_information
                # Initialize the clean map.
                clean_map = algebra.info_array(sp.zeros(dirty_map.shape))
                clean_map.info = dict(dirty_map.info)
                clean_map = algebra.make_vect(clean_map)
                # Make the clean map.
                clean_map[good_data] = (dirty_map[good_data] 
                                        / noise_inv_memory[good_data])
            elif noise_inv.ndim == 5 :
                pass
            elif noise_inv.ndim == 6 :
                raise NotImplementedError("Full noise matrix not yet "
                                          "implemented.  Best we can do is "
                                          "block diagonal in frequency.")
            else :
                raise ce.DataError("Noise matrix has bad shape.")




