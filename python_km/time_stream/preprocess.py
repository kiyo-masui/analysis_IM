"""This module performs time stream proprocessing.

In this module, only perform operations that you would probably want to do
in every concievable pipline.
"""

import numpy.ma as ma

from kiyopy import parse_ini
#from kiyopy import custom_exceptions as ce
import core.data_block

params_init = {
              # All parameters prefixed with pp_
              # IO:
              'pp_input_root' : './',
              # List of strings that get concatenated with input_root and ends.
              'pp_input_middles' : ["GBTdata"],
              'pp_input_ends' : ".raw.acs.fits",
              'pp_output_dir' : "./",
              'pp_output_root' : "pp_",
              # What data to process within each file.
              'pp_scans' : [],
              'pp_IFs' : [],
              # What to do:
              'pp_hanning' : True
              }

def preprocess(parameter_ini) :
    """Preprocesses GBT data.

    Driver that performs the variouse operations included in this module on
    data.
    """
    pass




def hanning(Data) :
    """Perform Hanning smoothing.

    This function accepts a DataBlock class and returns nothing.  It changes
    the data attribute of the passed Data.
    """
    
    # Perform Hanning smoothing frequency space convolution.
    # Masked array automatically flags data adjacent to masked
    # data.
    if len(Data.axes) != 4 or Data.axes[-1] != 'freq' :
        raise DataError('Hanning smoothing expects data be 4D with frequency '
                        'being the last axis.') 
    Data.data[:,:,:,1:-1] = ( 0.25*Data.data[:,:,:,:-2] + 
                         0.50*Data.data[:,:,:,1:-1] + 
                         0.25*Data.data[:,:,:,2:] )
    # End points where not smoothed.
    Data.data[:,:,:,0] = ma.masked
    Data.data[:,:,:,-1] = ma.masked

def pipe_driver(parameter_ini) :
    """Conventionally named function that just calls the preprocess function.
    
    A function named pipe_driver is expected by the pipline module.  It just
    pointes to the main function of this module.
    """
    preprocess(parameter_ini)
    


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    preprocess(string(sys.argv[1]))

