"""Converter from a dirty map to a clean map."""

import sys
import glob

import scipy as sp
import scipy.linalg as linalg

from core import algebra, hist
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
import _cholesky as _c


params_init = {'input_root' : './',
               'polarizations' : ('I',),
               'output_root' : './',
               'save_noise_diag' : False,
               'bands' : ()
               }
prefix = 'cm_'

class CleanMapMaker(object) :
    """Converts a Dirty map to a clean map."""

    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix=prefix, feedback=feedback)
        self.feedback = feedback

    def execute(self, nprocesses=1) :
        """Worker funciton."""
        params = self.params
        # Make parent directory and write parameter file.
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        save_noise_diag = params['save_noise_diag']
        in_root = params['input_root']
        all_out_fname_list = []
        all_in_fname_list = []
        # Figure out what the band names are.
        bands = params['bands']
        if not bands:
            map_files = glob.glob(in_root + 'dirty_map_' + pol_str + "_*.npy")
            bands = []
            root_len = len(in_root + 'dirty_map_')
            for file_name in map_files:
                bands.append(file_name[root_len:-4])
        # Loop over files to process.
        for pol_str in params['polarizations']:
            for band in bands:
                dmap_fname = (in_root + 'dirty_map_' + pol_str + "_" +
                              repr(band) + '.npy')
                noise_fname = (in_root + 'noise_inv_' + pol_str + "_" +
                               repr(band) + '.npy')
                if self.feedback > 1:
                    print "Using dirty map: " + dmap_fname
                    print "Using noise inverse: " + noise_fname
                all_in_fname_list.append(
                    kiyopy.utils.abbreviate_file_path(dmap_fname))
                all_in_fname_list.append(
                    kiyopy.utils.abbreviate_file_path(noise_fname))
                # Load the dirty map and the noise matrix.
                dirty_map = algebra.load(dmap_fname)
                dirty_map = algebra.make_vect(dirty_map)
                if dirty_map.axes != ('freq', 'ra', 'dec') :
                    msg = ("Expeced dirty map to have axes ('freq',"
                           "'ra', 'dec'), but it has axes: "
                           + str(dirty_map.axes))
                    raise ce.DataError(msg)
                shape = dirty_map.shape
                noise_inv = algebra.open_memmap(noise_fname, 'r')
                noise_inv = algebra.make_mat(noise_inv)
                # Initialize the clean map.
                clean_map = algebra.info_array(sp.zeros(dirty_map.shape))
                clean_map.info = dict(dirty_map.info)
                clean_map = algebra.make_vect(clean_map)
                # If needed, initialize a map for the noise diagonal.
                if save_noise_diag :
                    noise_diag = algebra.zeros_like(clean_map)
                # Two cases for the noise.  If its the same shape as the map
                # then the noise is diagonal.  Otherwise, it should be
                # block diagonal in frequency.
                if noise_inv.ndim == 3 :
                    if noise_inv.axes != ('freq', 'ra', 'dec') :
                        msg = ("Expeced noise matrix to have axes "
                                "('freq', 'ra', 'dec'), but it has: "
                                + str(noise_inv.axes))
                        raise ce.DataError(msg)
                    # Noise inverse can fit in memory, so copy it.
                    noise_inv_memory = sp.array(noise_inv, copy=True)
                    # Find the non-singular (covered) pixels.
                    max_information = noise_inv_memory.max()
                    good_data = noise_inv_memory < 1.0e-10*max_information
                    # Make the clean map.
                    clean_map[good_data] = (dirty_map[good_data] 
                                            / noise_inv_memory[good_data])
                    if save_noise_diag :
                        noise_diag[good_data] = 1/noise_inv_memory[good_data]
                elif noise_inv.ndim == 5 :
                    if noise_inv.axes != ('freq', 'ra', 'dec', 'ra', 'dec') :
                        msg = ("Expeced noise matrix to have axes "
                               "('freq', 'ra', 'dec', 'ra', 'dec'), "
                               "but it has: " + str(noise_inv.axes))
                        raise ce.DataError(msg)
                    # Arrange the dirty map as a vector.
                    dirty_map_vect = sp.array(dirty_map) # A view.
                    dirty_map_vect.shape = (shape[0], shape[1]*shape[2])
                    frequencies = dirty_map.get_axis('freq')/1.0e6
                    # Allowcate memory only once.
                    noise_inv_freq = sp.empty((shape[1], shape[2], shape[1],
                                               shape[2]), dtype=float)
                    if self.feedback > 1 :
                        print "Inverting noise matrix."
                    # Block diagonal in frequency so loop over frequencies.
                    for ii in xrange(dirty_map.shape[0]) :
                        if self.feedback > 1:
                            print "Frequency: ", "%5.1f"%(frequencies[ii]),
                        if self.feedback > 2:
                            print ", start mmap read:",
                            sys.stdout.flush()
                        noise_inv_freq[...] = noise_inv[ii, ...]
                        if self.feedback > 2:
                            print "done, start eig:",
                            sys.stdout.flush()
                        noise_inv_freq.shape = (shape[1]*shape[2],
                                                shape[1]*shape[2])
                        # Solve the map making equation by diagonalization.
                        noise_inv_diag, Rot = sp.linalg.eigh(noise_inv_freq, 
                                                             overwrite_a=True)
                        if self.feedback > 2:
                            print "done",
                        map_rotated = sp.dot(Rot.T, dirty_map_vect[ii])
                        # Zero out infinite noise modes.
                        bad_modes = (noise_inv_diag
                                     < 1.0e-5 * noise_inv_diag.max())
                        if self.feedback > 1:
                            print ", discarded: ",
                            print "%4.1f" % (100.0 * sp.sum(bad_modes) 
                                             / bad_modes.size),
                            print "% of modes",
                        if self.feedback > 2:
                            print ", start rotations:",
                            sys.stdout.flush()
                        map_rotated[bad_modes] = 0.
                        noise_inv_diag[bad_modes] = 1.0
                        # Solve for the clean map and rotate back.
                        map_rotated /= noise_inv_diag
                        map = sp.dot(Rot, map_rotated)
                        if self.feedback > 2:
                            print "done",
                            sys.stdout.flush()
                        # Fill the clean array.
                        map.shape = (shape[1], shape[2])
                        clean_map[ii, ...] = map
                        if save_noise_diag :
                            # Using C = R Lambda R^T 
                            # where Lambda = diag(1/noise_inv_diag).
                            temp_noise_diag = 1/noise_inv_diag
                            temp_noise_diag[bad_modes] = 0
                            # Multiply R by the diagonal eigenvalue matrix.
                            # Broadcasting does equivalent of mult by diag
                            # matrix.
                            temp_mat = Rot*temp_noise_diag
                            # Multiply by R^T, but only calculate the diagonal
                            # elements.
                            for jj in range(shape[1]*shape[2]) :
                                temp_noise_diag[jj] = sp.dot(temp_mat[jj,:], 
                                                             Rot[jj,:])
                            temp_noise_diag.shape = (shape[1], shape[2])
                            noise_diag[ii, ...] = temp_noise_diag
                        # Return workspace memory to origional shape.
                        noise_inv_freq.shape = (shape[1], shape[2],
                                                shape[1], shape[2])
                        if self.feedback > 1:
                            print ""
                            sys.stdout.flush()
                elif noise_inv.ndim == 6 :
                    if save_noise_diag:
                        clean_map, noise_diag = solve(noise_inv, dirty_map,
                                                True, feedback=self.feedback)
                    else:
                        clean_map = solve(noise_inv, dirty_map, False,
                                          feedback=self.feedback)
                else :
                    raise ce.DataError("Noise matrix has bad shape.")
                # Write the clean map to file.
                out_fname = (params['output_root'] + 'clean_map_'
                             + pol_str + "_" + repr(band) + '.npy')
                if self.feedback > 1:
                    print "Writing clean map to: " + out_fname
                algebra.save(out_fname, clean_map)
                all_out_fname_list.append(
                    kiyopy.utils.abbreviate_file_path(out_fname))
                if save_noise_diag :
                    noise_diag_fname = (params['output_root'] + 'noise_diag_'
                                        + pol_str + "_" + repr(band) + '.npy')
                    algebra.save(noise_diag_fname, noise_diag)
                    all_out_fname_list.append(
                        kiyopy.utils.abbreviate_file_path(noise_diag_fname))
            # This needs to be added to the new dirty map maker before I can
            # add it here.
            # Finally update the history object.
            #history = hist.read(in_root + 'history.hist')
            #history.add('Read map and noise files:', all_in_fname_list)
            #history.add('Converted dirty map to clean map.',
            #            all_out_fname_list)
            #h_fname = params['output_root'] + "history.hist"
            #history.write(h_fname)


def solve(noise_inv, dirty_map, return_noise_diag=False, feedback=0):
    """Solve for the clean map.

    Matrix and vector passed in as algebra objects with no block diagonality.
    The system is solved using a GPU accelerated cholesky decomposition.
    Optionally, the diagonal on the noise matrix is returned.
    """
    
    # Put into the 2D matrix shape.
    expanded = noise_inv.view()
    side_size = noise_inv.shape[0] * noise_inv.shape[1] * noise_inv.shape[2]
    expanded.shape = (side_size,) * 2
    # Allowcate memory for the cholesky and copy the upper triangular data.
    if feedback > 1:
        print "Copying matrix."
    tri_copy = _c.up_tri_copy(expanded)
    # Cholesky decompose it.
    if feedback > 1:
        print "Cholesky decomposition."
    _c.call_cholesky(tri_copy)
    # Solve for the clean map.
    flat_map = dirty_map.view()
    flat_map.shape = (flat_map.size,)
    if feedback > 1:
        print "Solving for clean map."
    clean_map = linalg.cho_solve((tri_copy, False), flat_map)
    # Reshape and cast as a map.
    clean_map.shape = dirty_map.shape
    clean_map = algebra.as_alg_like(clean_map, dirty_map)
    if not return_noise_diag:
        return clean_map
    else:
        if feedback > 1:
            print "Getting noise diagonal."
        noise_diag = sp.empty(side_size, dtype=float)
        _c.inv_diag_from_chol(tri_copy, noise_diag)
        noise_diag.shape = dirty_map.shape
        noise_diag = algebra.as_alg_like(noise_diag, dirty_map)
        return clean_map, noise_diag



