#!/usr/bin/python
"""This module performs time stream preprocessing.

In this module, only perform operations that you would probably want to do
in every conceivable pipeline such as Hanning smoothing. Steps included here
should not be 'controversial'.
"""

import numpy.ma as ma

from kiyopy import parse_ini
from kiyopy.utils import mkdir_p
#from kiyopy import custom_exceptions as ce

from core import data_block
from core import fitsGBT

params_init = {
              # All parameters prefixed with pp_ for PreProcess.
              # IO:
              'pp_input_root' : './',
              # List of strings that get concatenated with input_root and ends.
              'pp_input_middles' : ["GBTdata"],
              'pp_input_end' : ".raw.acs.fits",
              'pp_output_dir' : "./",
              'pp_output_root' : "pp_",
              'pp_output_end' : ".fits",
              # What data to process within each file.
              'pp_scans' : [],
              'pp_IFs' : [],
              # What to do:
              'pp_hanning' : True
              }

def preprocess(parameter_file_or_dict) :
    """Preprocesses GBT data.

    Driver that performs the various operations included in this module on
    the data.
    """
    
    # Read in the parameters.
    params = parse_ini.parse(parameter_file_or_dict, params_init, 
                             checking=10*self.feedback + 2)
    mkdir_p(params['pp_output_dir'])
    parse_ini.write_params(params, params['pp_output_dir'] + '/params.ini')

    # Loop over the files to process.
    for file_middle in params['pp_input_middles'] :
        input_fname = (params['pp_input_root'] + file_middle +
                       params['pp_input_end'])
        output_fname = (params['pp_output_dir'] + '/' + file_middle +
                        params['pp_output_end'])
        Writer = fitsGBT.Writer()
        
        # Read in the data, and loop over data blocks.
        Reader = fitsGBT.Reader(input_fname)
        Blocks = Reader.read(params['pp_scans'], params['pp_IFs'])
        for Data in Blocks :
            # Now process the data.
            if params['pp_hanning'] :
                hanning(Data)

            # Check that the data is valid and add it to the writer.
            Data.verify()
            Writer.add_data(Data)

        # Finally write the data back to file.
        Writer.write(output_fname)

def pipe_driver(parameter_file_or_dict) :
    """Conventionally named function that calls the preprocess function.
    
    A function named pipe_driver is expected by the pipeline module.  It just
    points to the main function of this module.
    """
    preprocess(parameter_file_or_dict)

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
    # Channels adjacent to masked data automatically masked.
    Data.data[:,:,:,1:-1] = ( 0.25*Data.data[:,:,:,:-2] + 
                         0.50*Data.data[:,:,:,1:-1] + 
                         0.25*Data.data[:,:,:,2:] )
    # End points where not smoothed, mask them as bad data.
    Data.data[:,:,:,0] = ma.masked
    Data.data[:,:,:,-1] = ma.masked


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    preprocess(string(sys.argv[1]))

