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
        if rank == 0:
            print "NUMBER OF PROCESSES: %d" % nproc
            print "BLOCKSIZE: %.2f" % blocksize[0]
        # Number of processors per side on grid.
        # usually this should be square.
        npx = int(nproc**0.5)  # Round down intensionally.
        npy = npx
        pscore.initmpi(gridsize=[npx, npy], blocksize=blocksize)
        # Make parent directory and write parameter file.
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
                dirty_fname = (in_root + 'dirty_map_' + pol_str + "_" +
                               repr(band) + '.npy')
                if self.feedback > 1 and rank == 0:
                    print "Using noise inverse: " + noise_fname
                    print "Using dirty map: " + dirty_fname

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
                print "Making Distributed Noise Matrix."
                Mat = pscore.DistributedMatrix.from_npy(noise_fname,
                        blocksize=blocksize, shape_override=mat_shape)

                # Get the dirty map, too.
                d_shape, d_fortran_order, d_dtype, d_header_length = \
                        npyutils.read_header_data(dirty_fname)
                if len(d_shape) != 3 or d_fortran_order or d_dtype != '<f8':
                    msg = "Dirty map header not as expected."
                    raise ce.DataError(msg)
                if (d_shape[0] != shape[0] or d_shape[1] != shape[1]
                    or d_shape[2] != shape[2]):
                    msg = "First three dimensions of Noise matrix and "
                    msg += "Dirty map must be the same."
                    raise ce.DataError(msg)
                d_n = d_shape[0] * d_shape[1] * d_shape[2]
                d_mat_shape = (d_n, 1)
                print "Making Distributed Dirty Map Vector/Matrix."
                Rhs = pscore.DistributedMatrix.from_npy(dirty_fname,
                        blocksize=blocksize, shape_override=d_mat_shape)

                # Factorize. Upper triangular by default.
                print "Chol decomposition."
                chol = psroutine.pdpotrf(Mat, destroy=True)
                print "Done decomposition."
                # Solve. Assumes upper triangular factorization by default.
                print "Solving for clean map."
                clean_map = psroutine.pdpotrs(chol,Rhs,destroy=True)
                print "Done solving."

                #evals = np.zeros(n)
                #evects = Mat
                # Now construct meta data for these objects.
                #evals = algebra.make_mat(evals, axis_names=('mode',),
                #                         row_axes=(0,), col_axes=(0,))
                #evects_mat_shape = (shape[0], shape[1], shape[2], n)
                #evects_mat_axis_names = ('freq', 'ra', 'dec', 'mode')
                #evects_mat_rows = (0, 1, 2)
                #evects_mat_cols = (3,)
                #evects_mat_info = {'type': 'mat',
                #                   'axes': evects_mat_axis_names,
                #                   'rows': evects_mat_rows,
                #                   'cols': evects_mat_cols}
                # Figure out file names:
                chol_fname = (out_root + 'noise_inv_chol_' + pol_str + "_" +
                               repr(band) + "_" + str(blocksize[0]) + "_" + str(nproc) + '.npy')
                evals_fname = (out_root + 'noise_evalsinv_' + pol_str + "_" +
                               repr(band) + '.npy')
                clean_fname = (out_root + 'clean_map_' + pol_str + "_" +
                               repr(band) + "_" + str(blocksize[0]) + "_" + str(nproc) + '.npy')

                # Write out.
      #          if rank == 0 and self.feedback > 1:
      #              print "Writing Cholesky to " + chol_fname
      #          chol.to_npy(chol_fname, fortran_order=False)#,
#                   shape_override=evects_mat_shape)
                if rank == 0 and self.feedback > 1:
                    print "Writing clean map to " + clean_fname
                clean_map.to_npy(clean_fname, fortran_order=False,
                    shape_override=d_shape)

                # Write meta data.
                if rank == 0:
#                    if self.feedback > 1:
#                        print "Writing cholesky to " + chol_fname
#                        print "Writing clean map to " + clean_fname
#                    # Eigenvalues.
#                    #algebra.save(evals_fname, evals)
#                    # Eigenvector meta data.
#                    meta_fid = open(clean_fname + ".meta", 'w')
#                    try :
#                        clean_map_info = clean_map.info
#                        meta_fid.write(repr(clean_map_info))
#                        #placeholder_uselessness='just_to_fill_in_the_try_block'
#                    finally :
#                        meta_fid.close()
                    if self.feedback > 1:
                        print "Writing clean map meta"
                    try:
                        d_meta_fid = open(dirty_fname + ".meta", 'r')
                        meta_fid = open(clean_fname + ".meta", 'w')
                        meta_fid.write(d_meta_fid.read())
                    finally:
                        meta_fid.close()
                        d_meta_fid.close()
                    if self.feedback > 1:
                        print "Done."
                    

#                # Cholesky.
#                chol.to_npy(chol_fname, fortran_order=False)#,
#                              shape_override=evects_mat_shape)
                


		    
# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    EigNoise(str(sys.argv[1])).execute()


