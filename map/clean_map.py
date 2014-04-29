"""Converter from a dirty map to a clean map."""

import sys
import glob
import time

import numpy as np
import scipy as sp
import scipy.linalg as linalg
from ast import literal_eval as ev

from core import algebra, hist
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
import constants
import h5py


params_init = {'input_root' : './',
               'polarizations' : ('I',),
               'output_root' : './',
               'save_noise_diag' : False,
               'save_noise_inv_diag' : False,
               'save_cholesky' : False,
               'from_eig' : False,
               'bands' : ()
               'mem_lim' : ()
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
                if band == -1:
                    band_str = ''
                else:
                    band_str =  "_" + repr(band)
                dmap_fname = (in_root + 'dirty_map_' + pol_str + 
                              band_str + '.npy')
                all_in_fname_list.append(
                    kiyopy.utils.abbreviate_file_path(dmap_fname))
                # Load the dirty map and the noise matrix.
                dirty_map = algebra.load(dmap_fname)
                dirty_map = algebra.make_vect(dirty_map)
                if dirty_map.axes != ('freq', 'ra', 'dec') :
                    msg = ("Expeced dirty map to have axes ('freq',"
                           "'ra', 'dec'), but it has axes: "
                           + str(dirty_map.axes))
                    raise ce.DataError(msg)
                shape = dirty_map.shape
                # Initialize the clean map.
                clean_map = algebra.info_array(sp.zeros(dirty_map.shape))
                clean_map.info = dict(dirty_map.info)
                clean_map = algebra.make_vect(clean_map)
                # If needed, initialize a map for the noise diagonal.
                if save_noise_diag :
                    noise_diag = algebra.zeros_like(clean_map)
                if params["from_eig"]:
                    # Solving from eigen decomposition of the noise instead of
                    # the noise itself.
                    # Load in the decomposition.
                    evects_fname = (in_root + 'noise_evects_' + pol_str +
                                    + band_str + '.npy')
                    if self.feedback > 1:
                        print "Using dirty map: " + dmap_fname
                        print "Using eigenvectors: " + evects_fname
                    evects = algebra.open_memmap(evects_fname, 'r')
                    evects = algebra.make_mat(evects)
                    evals_inv_fname = (in_root + 'noise_evalsinv_' + pol_str
                                       + "_" + repr(band) + '.npy')
                    evals_inv = algebra.load(evals_inv_fname)
                    evals_inv = algebra.make_mat(evals_inv)
                    # Solve for the map.
                    if params["save_noise_diag"]:
                        clean_map, noise_diag = solve_from_eig(evals_inv,
                                    evects, dirty_map, True, self.feedback)
                    else:
                        clean_map = solve_from_eig(evals_inv,
                                    evects, dirty_map, False, self.feedback)
                    # Delete the eigen vectors to recover memory.
                    del evects
                else:
                    # Solving from the noise.
                    noise_fname1 = (in_root + 'noise_inv_' + pol_str +
                                   band_str + '.npy')
                    noise_fname2 = (in_root + 'noise_inv_' + pol_str +
                                   band_str + '.hdf5')
                    noise_fnames = glob.glob(noise_fname1)+glob.glob(noise_fname2)
                    if len(noise_fnames) == 0:
                        raise ValueError("Couldn't find noise_inv file.")
                    if len(noise_fnames) > 1:
                        raise ValueError("Multiple files have naming pattern of noise_inv")
                    noise_fname = noise_fnames[0]
                    if self.feedback > 1:
                        print "Using dirty map: " + dmap_fname
                        print "Using noise inverse: " + noise_fname
                    all_in_fname_list.append(
                        kiyopy.utils.abbreviate_file_path(noise_fname))

                    if noise_fname.split('.')[-1] == 'npy':
                        noise_inv = algebra.open_memmap(noise_fname, 'r')
                        noise_inv = algebra.make_mat(noise_inv)
                    elif noise_fname.split('.')[-1] == 'hdf5':
                        noise_h5 = h5py.File(noise_fname, 'r')
                        noise_inv = algebra.load_h5_memmap(noise_h5, 'inv_cov', noisefname + ".npy" , params["mem_lim"])
                        noise_inv = algebra.make_mat(noise_inv)
                    else:
                        raise ValueError("Noise file is of unsupported type, neither .npy or .hdf5.")
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
                            noise_diag[good_data] = \
                                    1/noise_inv_memory[good_data]
                    elif noise_inv.ndim == 5 :
                        if noise_inv.axes != ('freq', 'ra', 'dec', 'ra',
                                              'dec'):
                            msg = ("Expeced noise matrix to have axes "
                                   "('freq', 'ra', 'dec', 'ra', 'dec'), "
                                   "but it has: " + str(noise_inv.axes))
                            raise ce.DataError(msg)
                        # Arrange the dirty map as a vector.
                        dirty_map_vect = sp.array(dirty_map) # A view.
                        dirty_map_vect.shape = (shape[0], shape[1]*shape[2])
                        frequencies = dirty_map.get_axis('freq')/1.0e6
                        # Allowcate memory only once.
                        noise_inv_freq = sp.empty((shape[1], shape[2], 
                                        shape[1], shape[2]), dtype=float)
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
                            noise_inv_diag, Rot = sp.linalg.eigh(
                                noise_inv_freq, overwrite_a=True)
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
                                # Multiply by R^T, but only calculate the
                                # diagonal elements.
                                for jj in range(shape[1]*shape[2]) :
                                    temp_noise_diag[jj] = sp.dot(
                                        temp_mat[jj,:], Rot[jj,:])
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
                            # OLD WAY.
                            #clean_map, noise_diag, chol = solve(noise_inv,
                            #        dirty_map, True, feedback=self.feedback)
                            # NEW WAY.
                            clean_map, noise_diag, noise_inv_diag, chol = \
                                      solve(noise_fname, noise_inv, dirty_map,
                                      True, feedback=self.feedback)
                        else:
                            # OLD WAY.
                            #clean_map, chol = solve(noise_inv, dirty_map, 
                            #            False, feedback=self.feedback)
                            # NEW WAY.
                            clean_map, noise_inv_diag, chol = \
                                      solve(noise_fname, noise_inv, dirty_map,
                                      False, feedback=self.feedback)
                        if params['save_cholesky']:
                            chol_fname = (params['output_root'] + 'chol_'
                                        + pol_str + band_str + '.npy')
                            sp.save(chol_fname, chol)
                        if params['save_noise_inv_diag']:
                            noise_inv_diag_fname = (params['output_root'] +
                                       'noise_inv_diag_' + pol_str + band_str 
                                       + '.npy')
                            algebra.save(noise_inv_diag_fname, noise_inv_diag)
                        # Delete the cholesky to recover memory.
                        del chol
                    else :
                        raise ce.DataError("Noise matrix has bad shape.")
                    # In all cases delete the noise object to recover memeory.
                    del noise_inv
                # Write the clean map to file.
                out_fname = (params['output_root'] + 'clean_map_'
                             + pol_str + band_str + '.npy')
                if self.feedback > 1:
                    print "Writing clean map to: " + out_fname
                algebra.save(out_fname, clean_map)
                all_out_fname_list.append(
                    kiyopy.utils.abbreviate_file_path(out_fname))
                if save_noise_diag :
                    noise_diag_fname = (params['output_root'] + 'noise_diag_'
                                        + pol_str + band_str + '.npy')
                    algebra.save(noise_diag_fname, noise_diag)
                    all_out_fname_list.append(
                        kiyopy.utils.abbreviate_file_path(noise_diag_fname))
                # Check the clean map for faileur.
                if not sp.alltrue(sp.isfinite(clean_map)):
                    n_bad = sp.sum(sp.logical_not(sp.isfinite(clean_map)))
                    msg = ("Non finite entries found in clean map. Solve"
                           " failed. %d out of %d entries bad" 
                           % (n_bad, clean_map.size)) 
                    raise RuntimeError(msg)
            # This needs to be added to the new dirty map maker before I can
            # add it here.
            # Finally update the history object.
            #history = hist.read(in_root + 'history.hist')
            #history.add('Read map and noise files:', all_in_fname_list)
            #history.add('Converted dirty map to clean map.',
            #            all_out_fname_list)
            #h_fname = params['output_root'] + "history.hist"
            #history.write(h_fname)

def solve(noise_inv_filename, noise_inv, dirty_map, 
                                         return_noise_diag=False, feedback=0):
    """Solve for the clean map.

    Matrix and vector passed in as algebra objects with no block diagonality.
    The system is solved using a GPU accelerated cholesky decomposition.
    Optionally, the diagonal on the noise matrix is returned.
    """
    
    # Import the cython stuff locally so some clean maps can be made on any
    # machine.
    import _cholesky as _c
    # Put into the 2D matrix shape.
    expanded = noise_inv.view()
    side_size = noise_inv.shape[0] * noise_inv.shape[1] * noise_inv.shape[2]
    expanded.shape = (side_size,) * 2
    # Allowcate memory for the cholesky and copy the upper triangular data.
    # Instead of copying, open the matrix in copy on write mode.
    if feedback > 1:
        print "Copying matrix."
    time_before = time.time() / 60.
    # OLD WAY.
    # tri_copy = _c.up_tri_copy(expanded)
    # NEW WAY.
    tri_copy = up_tri_copy_from_file(noise_inv_filename)
    time_after = time.time() / 60.
    if feedback > 1:
        print "\nMatrix copying time: %.2f minutes." % (time_after-time_before) 
    #tri_copy = expanded
    # Get the diagonal of the giant noise inverse.
    if feedback > 1:
        print "Getting noise inverse diagonal."
    time_before = time.time() / 60.
    tri_copy.shape = side_size*side_size
    noise_inv_diag = tri_copy[::side_size+1]
    tri_copy.shape = (side_size,side_size)
    noise_inv_diag.shape = dirty_map.shape
    noise_inv_diag = algebra.as_alg_like(noise_inv_diag, dirty_map)
    time_after = time.time() / 60.
    if feedback > 1:
        print "\nNoise inverse diagonal gotten in: %.2f minutes." \
                                                 % (time_after-time_before)
    # Cholesky decompose it.
    if feedback > 1:
        print "Cholesky decomposition."
    time_before = time.time() / 60.
    _c.call_cholesky(tri_copy)
    time_after = time.time() / 60.
    if feedback > 1:
        print "\nCholesky decomposition time: %.2f minutes." \
                                                 % (time_after-time_before) 
    # Solve for the clean map.
    flat_map = dirty_map.view()
    flat_map.shape = (flat_map.size,)
    if feedback > 1:
        print "Solving for clean map."
    time_before = time.time() / 60.
    clean_map = _c.cho_solve(tri_copy, flat_map)
    time_after = time.time() / 60.
    if feedback > 1:
        print "\nClean map solving time: %.2f minutes." % (time_after-time_before) 
    # Reshape and cast as a map.
    clean_map.shape = dirty_map.shape
    clean_map = algebra.as_alg_like(clean_map, dirty_map)
    if not return_noise_diag:
        return clean_map, noise_inv_diag, tri_copy
    else:
        if feedback > 1:
            print "Getting noise diagonal."
        noise_diag = sp.empty(side_size, dtype=float)
        time_before = time.time() / 60.
        _c.inv_diag_from_chol(tri_copy, noise_diag)
        time_after = time.time() / 60.
        if feedback > 1:
            print "\nNoise diagonal gotten in: %.2f minutes." \
                                                 % (time_after-time_before) 
        noise_diag.shape = dirty_map.shape
        noise_diag = algebra.as_alg_like(noise_diag, dirty_map)
        return clean_map, noise_diag, noise_inv_diag, tri_copy

def solve_from_eig(noise_evalsinv, noise_evects, dirty_map,
                   return_noise_diag=False, feedback=0):
    """Converts a dirty map to a clean map using the eigen decomposition of the
    noise inverse.
    """
    
    # Check the shapes.
    if noise_evects.ndim != 4:
        raise ValueError("Expected 4D array for 'noise_evects`.")
    if noise_evalsinv.shape != (noise_evects.shape[-1],):
        raise ValueError("Wrong number of eigenvalues.")
    if dirty_map.shape != noise_evects.shape[:-1]:
        raise ValueError("Dirty map and noise don't have matching dimensions.")
    if dirty_map.size != noise_evects.shape[-1]:
        raise ValueError("Eigen space not the same total size as map space.")
    n = noise_evalsinv.shape[0]
    nf = dirty_map.shape[0]
    nr = dirty_map.shape[1]
    nd = dirty_map.shape[2]
    # Copy the eigenvalues.
    noise_evalsinv = noise_evalsinv.copy()
    # Find poorly constrained modes and zero them out.
    bad_inds = noise_evalsinv < 1./constants.T_huge**2
    n_bad = sp.sum(bad_inds)
    if feedback > 1:
        print ("Discarding %d modes of %d. %f percent." 
               % (n_bad, n, 100. * n_bad / n))
    noise_evalsinv[bad_inds] = 1.
    # Rotate the dirty map into the diagonal noise space.
    if feedback > 1:
        print "Rotating map to eigenspace."
    map_rot = sp.zeros(n, dtype=sp.float64)
    for ii in xrange(nf):
        for jj in xrange(nr):
            for kk in xrange(nd):
                tmp = noise_evects[ii,jj,kk,:].copy()
                map_rot += dirty_map[ii,jj,kk] * tmp
    # Multiply by the (diagonal) noise (inverse, inverse).  Zero out any poorly
    # constrained modes.
    map_rot[bad_inds] = 0
    # Take inverse and multiply.
    map_rot = map_rot / noise_evalsinv
    # Now rotate back to the origional space.
    if feedback > 1:
        print "Rotating back to map space."
    clean_map = algebra.zeros_like(dirty_map)
    for ii in xrange(nf):
        for jj in xrange(nr):
            for kk in xrange(nd):
                tmp = noise_evects[ii,jj,kk,:].copy()
                clean_map[ii,jj,kk] = sp.sum(map_rot * tmp)
    if return_noise_diag:
        if feedback > 1:
            print "Getting noise diagonal."
        noise_diag = algebra.zeros_like(dirty_map)
        noise_evals = 1. / noise_evalsinv
        noise_evals[bad_inds] = constants.T_huge**2
        for ii in xrange(nf):
            for jj in xrange(nr):
                for kk in xrange(nd):
                    tmp = noise_evects[ii,jj,kk,:].copy()
                    noise_diag[ii,jj,kk] = sp.sum(tmp**2 * noise_evals)
        return clean_map, noise_diag
    else:
        return clean_map


def read_n_numbers(f,arr,dt,n,row_len,row_pos):
    '''Read n numbers from an open file stream, f, and add them to the next
    available row, row_pos, in the 2D array, arr. The data type, dt, is needed
    to convert from bytes to usable numbers. row_len and row_pos are needed
    so that only the upper triangular part of a row is added to arr.'''
    bytes_per_num = 0
    # Assign how long a number is in memory.
    if ((dt.name == 'int64') or (dt.name == 'float64')):
        bytes_per_num = 8
    if ((dt.name == 'int32') or (dt.name == 'float32')):
        bytes_per_num = 4
    # Check if possible to run.
    if (bytes_per_num == 0):
        msg =  "%s is not a supported number type right now." % (dt.name)
        raise ce.DataError(msg)
    # Read in n numbers.
    raw_nums = f.read(bytes_per_num*n)
    # Convert numbers in bytes to useable numbers all at once.
    new_slice = np.fromstring(raw_nums,dtype=dt,count=n)    
    # n is a multiple of rows.
    rows_to_read = int(n/row_len)
    # Copy each row from the buffer (in memory) to the array (in memory).
    # We are copying in only the upper triangular part of the array so
    # we must only write in the appropriate values.
    for i in range(rows_to_read):
        # Every row down the matrix, we want to ignore one more number
        # at the beginning of a row.
        arr[row_pos+i,row_pos+i:row_len] = \
             new_slice[i*row_len+row_pos+i:i*row_len+row_len]


def up_tri_copy_from_file(filename):
    '''Return the upper triagular part of the array in filename. The file
    must be a binary .npy file.'''
    f = open(filename, 'rb')
    # The first 6 bytes are a magic string: exactly "\x93NUMPY".
    numpy_string = f.read(6)
    # The next 1 byte is an unsigned byte: the major version number
    # of the file format, e.g. \x01.
    major_ver = ord(f.read(1))
    # The next 1 byte is an unsigned byte: the minor version number
    # of the file format, e.g. \x00.
    minor_ver = ord(f.read(1))
    # Check that the version of the file used is what this code can handle.
    if ((numpy_string != "\x93NUMPY") or (major_ver != 1) or (minor_ver != 0)):
        msg = "Array can only be read from a '\93NUMPY' version 1.0 .npy file "
        msg += "not %s %d.%d" % (numpy_string, major_ver, minor_ver)
        raise ce.DataError(msg)
    # The next 2 bytes form a little-endian unsigned short int: the
    # length of the header data HEADER_LEN.
    byte2,byte1 = f.read(2)
    # Get the value from a short int (2 bytes) to decimal.
    header_len = (16**2)*ord(byte1) + ord(byte2)
    # The next HEADER_LEN bytes form the header data describing the
    # array's format. It is an ASCII string which contains a Python
    # literal expression of a dictionary. It is terminated by a newline
    # ('\n') and padded with spaces ('\x20') to make the total length of
    # the magic string + 4 + HEADER_LEN be evenly divisible by 16 for
    # alignment purposes.
    raw_dict = f.read(header_len)
    # Take off trailing uselessness.
    raw_dict = raw_dict.rstrip('\n')
    raw_dict = raw_dict.rstrip()
    header_dict = ev(raw_dict)
    # Check that it is possible to read array.
    # "fortran_order" : bool
    #     Whether the array data is Fortran-contiguous or not.
    if (header_dict['fortran_order']):
        msg = "Array must be in C order, not fortran order."
        raise ce.DataError(msg)
    # "shape" : tuple of int
    #     The shape of the array.
    arr_shape = header_dict['shape']
    # "descr" : dtype.descr
    #     An object that can be passed as an argument to the
    #     numpy.dtype() constructor to create the array's dtype.
    dt = np.dtype(header_dict['descr'])
    # Where to save all the data. Note this does not make a new array in memory.
    num_nums = np.product(arr_shape)
    # Assume array is square
    row_len = int(np.sqrt(num_nums))
    num_rows = row_len
    loaded_arr = np.empty((row_len,row_len))
    print "Rows copied:"
    # Load the array reading n numbers at a time. Reading one number at
    # a time is best for memory but requires too many disk seeks for a
    # memory map. Loading too many numbers at once requires too much memory,
    # so choose a happy medium.
    # Assume array is square. Reads n rows of numbers at a time.
    rows_to_read = 1000
    n = rows_to_read*row_len
    row_pos = 0
    while (row_pos < num_rows):
        # n may not divide evenly into the total number of numbers.
        # This statement can only be True at the end of the array.
        if ((num_rows-row_pos) < rows_to_read):
            rows_to_read = num_rows-row_pos
            n = rows_to_read*row_len
        read_n_numbers(f,loaded_arr,dt,n,row_len,row_pos)
        # Changing row_pos in the function called does not change it here.
        row_pos += rows_to_read
        sys.stderr.write(str(row_pos)+' ')
    print
    f.close()
    return loaded_arr

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        par_file = sys.argv[1]
        nproc = 1
        CleanMapMaker(par_file).execute(nproc)
