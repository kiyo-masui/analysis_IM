import numpy as np
from mpi4py import MPI

from pyscalapack import core as pscore
from pyscalapack import npyutils as npyutils
from pyscalapack import routines as psroutine

from core import algebra
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce

params_init = {'input_root' : './',
               'polarizations' : ('I',),
			   'bands' : (),
			   'output_root' : './',
			   'blocksize' : 512
              }

prefix = 'en_'

# For now.
comm = MPI.COMM_WORLD

class EigNoise(object):

    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Only have the first node report.
        rank = comm.Get_rank()
        if rank == 0:
            pass
        else:
            feedback = 0
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix=prefix, feedback=feedback)
        self.feedback = feedback

    def execute(self, nprocesses=1):
        #nprocesses doesn't do anything.
	    
        params = self.params
        # Matrix blocking.
        blocksize = (params['blocksize'], params['blocksize'])
        # Set up mpi and the process grid.
        nproc = comm.Get_size()
        rank = comm.Get_rank()
        # Number of processors per side on grid.
        # usually this should be square.
        npx = int(nproc**0.5)  # Round down intensionally.
        npy = npx
        pscore.initmpi(gridsize=[npx, npy], blocksize=blocksize)
        # Make parent directory and write parameter file.
        if rank == 0:
            kiyopy.utils.mkparents(params['output_root'])
            parse_ini.write_params(params, params['output_root'] + 'params.ini',
                                   prefix=prefix)
        in_root = params['input_root']
        out_root = params['output_root']
        # Figure out what the band names are.
        bands = params['bands']
        for pol_str in params['polarizations']:
            for band in bands:
                noise_fname = (in_root + 'noise_inv_' + pol_str + "_" +
                               repr(band) + '.npy')
                if self.feedback > 1 and rank == 0:
                    print "Using noise inverse: " + noise_fname
                # Read the array header as there is some information it it we
                # need.
                shape, fortran_order, dtype, header_length = \
                        npyutils.read_header_data(noise_fname)
                print "shape: %s || fortran_order: %s || dtype: %s" % (repr(shape),repr(fortran_order),repr(dtype))
                if len(shape) != 6 or fortran_order or dtype != '<f8':
                    msg = "Noise matrix header not as expected."
                    raise ce.DataError(msg)
                if (shape[0] != shape[3] or shape[1] != shape[4]
                    or shape[2] != shape[5]):
                    msg = "Noise matrix not square."
                    raise ce.DataError(msg)
                n = shape[0] * shape[1] * shape[2]
                mat_shape = (n, n)
                print "Making Distributed Matrix."
                Mat = pscore.DistributedMatrix.from_npy(noise_fname,
                        blocksize=blocksize, shape_override=mat_shape)
                # Replace with a Scalapack call.
                print "Eval decomposition."
                evals, evects = psroutine.pdsyevd(Mat, destroy=True)
                print "Done decomposition."
                #evals = np.zeros(n)
                #evects = Mat
                evals, evects = psroutine.pdsyevd(Mat, destroy=True)
                # Now construct meta data for these objects.
                evals = algebra.make_mat(evals, axis_names=('mode',),
                                         row_axes=(0,), col_axes=(0,))
                evects_mat_shape = (shape[0], shape[1], shape[2], n)
                evects_mat_axis_names = ('freq', 'ra', 'dec', 'mode')
                evects_mat_rows = (0, 1, 2)
                evects_mat_cols = (3,)
                evects_mat_info = {'type': 'mat',
                                   'axes': evects_mat_axis_names,
                                   'rows': evects_mat_rows,
                                   'cols': evects_mat_cols}
                # Figure out file names:
                evects_fname = (out_root + 'noise_evects_' + pol_str + "_" +
                               repr(band) + '.npy')
                evals_fname = (out_root + 'noise_evalsinv_' + pol_str + "_" +
                               repr(band) + '.npy')
                # Write out.
                if rank == 0:
                    if self.feedback > 1:
                        print "Writing eigenvectors to " + evects_fname
                    # Eigenvalues.
                    algebra.save(evals_fname, evals)
                    # Eigenvector meta data.
                    meta_fid = open(evects_fname + ".meta", 'w')
                    try :
                        meta_fid.write(repr(evects_mat_info))
                    finally :
                        meta_fid.close()
                # Eigenvectors.
                evects.to_npy(evects_fname, fortran_order=False,
                              shape_override=evects_mat_shape)


		    
# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    EigNoise(str(sys.argv[1])).execute()


