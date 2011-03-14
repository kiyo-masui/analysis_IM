"""Converter from a dirty map to a clean map."""

import scipy as sp
import scipy.linalg as linalg

from core import algebra
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce


params_init = {'input_root' : './',
               'polarizations' : ('I',),
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
        in_root = params['input_root']
        # Loop over files to process.
        for pol_str in params['polarizations']:
            dmap_fname = in_root + 'dirty_map_' + pol_str + '.npy'
            noise_fname = in_root + 'noise_inv_' + pol_str + '.npy'
            # Load the dirty map and the noise matrix.
            dirty_map = algebra.load(dmap_fname)
            dirty_map = algebra.make_vect(dirty_map)
            if dirty_map.axes != ('freq', 'ra', 'dec') :
                raise ce.DataError("Expeced dirty map to have axes "
                                   "('freq', 'ra', 'dec'), but it has axes: "
                                   + str(dirty_map.axes))
            shape = dirty_map.shape
            noise_inv = algebra.open_memmap(noise_fname, 'r')
            noise_inv = algebra.make_mat(noise_inv)
            # Initialize the clean map.
            clean_map = algebra.info_array(sp.zeros(dirty_map.shape))
            clean_map.info = dict(dirty_map.info)
            clean_map = algebra.make_vect(clean_map)
            # Two cases for the noise.  If its the same shape as the map then
            # the noise is diagonal.  Otherwise, it should be block diagonal in
            # frequency.
            if noise_inv.ndim == 3 :
                if noise_inv.axes != ('freq', 'ra', 'dec') :
                    raise ce.DataError("Expeced noise matrix to have axes "
                                       "('freq', 'ra', 'dec'), but it has: "
                                       + str(noise_inv.axes))
                # Noise inverse can fit in memory, so copy it.
                noise_inv_memory = sp.array(noise_inv, copy=True)
                # Find the non-singular (covered) pixels.
                max_information = noise_inv_memory.max()
                good_data = noise_inv_memory < 1.0e-5*max_information
                # Make the clean map.
                clean_map[good_data] = (dirty_map[good_data] 
                                        / noise_inv_memory[good_data])
            elif noise_inv.ndim == 5 :
                if noise_inv.axes != ('freq', 'ra', 'dec', 'ra', 'dec') :
                    raise ce.DataError("Expeced noise matrix to have axes "
                                       "('freq', 'ra', 'dec', 'ra', 'dec'), "
                                       "but it has: "
                                       + str(noise_inv.axes))
                # Arrange the dirty map as a vector.
                dirty_map_vect = sp.array(dirty_map) # A view.
                dirty_map_vect.shape = (shape[0], shape[1]*shape[2])
                # Block diagonal in frequency so loop over frequencies.
                for ii in xrange(dirty_map.shape[0]) :
                    noise_inv_freq = sp.array(noise_inv[ii, ...], copy=True)
                    noise_inv_freq.shape = (shape[1]*shape[2],
                                            shape[1]*shape[2])
                    # Solve the map making equation by diagonalization.
                    noise_inv_diag, Rot = sp.linalg.eigh(noise_inv_freq, 
                                                         overwrite_a=True)
                    map_rotated = sp.dot(Rot.T, dirty_map_vect[ii])
                    # Zero out infinite noise modes.
                    bad_modes = noise_inv_diag < 1.0e5*noise_inv_diag.max()
                    map_rotated[bad_modes] = 0.
                    noise_inv_diag[bad_modes] = 1.0
                    # Solve for the clean map and rotate back.
                    map_rotated /= noise_inv_diag
                    map = sp.dot(Rot, map_rotated)
                    # Fill the clean array.
                    map.shape = (shape[1], shape[2])
                    clean_map[ii, ...] = map
            elif noise_inv.ndim == 6 :
                raise NotImplementedError("Full noise matrix not yet "
                                          "implemented.  Best we can do is "
                                          "block diagonal in frequency.")
            else :
                raise ce.DataError("Noise matrix has bad shape.")
            # Write the clean map to file.
            out_fname = params['output_root'] + 'clean_map_' + pol_str + '.npy'
            algebra.save(out_fname, clean_map)

