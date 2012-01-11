import numpy as np
from mpi4py import MPI

from pyscalapack import core as pscore
from pyscalapack import routines as psroutine

from core import algebra
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce

params_init = {'input_root' : './',
               'polarizations' : ('I',),
			   'bands' : (),
			   'output_root' : './',
			   'blocking_factor' : 512
              }

prefix = 'en_'

class EigNoise(object):

    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix=prefix, feedback=feedback)
        self.feedback = feedback

    def execute(self, nprocesses=1):
        #nprocesses doesn't do anything.
	    
        # Set up mpi and the process grid.
        comm = MPI.COMM_WORLD
        nproc = comm.Get_size()
        rank = comm.Get_rank()
        # Number of processors per side on grid.
        # usually this should be square.
        npx = int(nproc**0.5)  # Round down intensionally.
        npy = npx
        params = self.params
        # Make parent directory and write parameter file.
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        in_root = params['input_root']        
        # Figure out what the band names are.
        bands = params['bands']
        for pol_str in params['polarizations']:
            for band in bands:
                noise_fname = (in_root + 'noise_inv_' + pol_str + "_" +
                               repr(band) + '.npy')
                if self.feedback > 1:
                    print "Using noise inverse: " + noise_fname
                noise_inv = algebra.open_memmap(noise_fname, 'r')
                noise_inv = algebra.make_mat(noise_inv)
                if noise_inv.axes != ('freq', 'ra', 'dec', 'freq', 'ra', 
                                      'dec'):
                    msg = ("Expeced noise matrix to have axes "
                            "('freq', 'ra', 'dec', 'freq', 'ra', 'dec),"
                            "but it has: " + str(noise_inv.axes))
                    raise ce.DataError(msg)
                print str(nproc) + "   " + str(rank)

		    
# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    EigNoise(str(sys.argv[1])).execute()


