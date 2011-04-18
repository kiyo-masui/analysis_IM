"""Measure the noise parameters of the data."""

import scipy as sp
import numpy.ma as ma

from kiyopy import parse_ini
import kiyopy.utils
import core.fitsGBT


params_init = {
               # IO.
               "input_root" : "./",
               "file_middles" : ("testfile_GBTfits",),
               "input_end" : ".fits",
               "output_root" : "./",
               # Algorithm.
               "model" : "scan_var"
               }

prefix = 'mn_'

class Measure(object) :
    """Measures the noise of data files.
    """

    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read the parameter file, store in dictionary named parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                      prefix=prefix, feedback=feedback)
        self.feedback = feedback

    def execute(self, nprocesses=1) :
        
        params = self.params
        model = params["model"]
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        # Loop over files to process.
        for file_middle in params['file_middles'] :
            input_fname = (params['input_root'] + file_middle +
                           params['input_end'])
            Reader = core.fitsGBT.Reader(input_fname, feedback=self.feedback)
            output_fname = params["output_root"] + file_middle + ".npy"
            if model == "scan_var" :
                n_scans = len(Reader.scan_set)
                n_IFs = len(Reader.IF_set)
                first_block = True
                for ii in range(n_scans) :
                    for jj in range(n_IFs) :
                        Data = Reader.read(ii, jj)
                        if first_block :
                            out_shape = (n_scans, n_IFs) + Data.dims[1:]
                            out_arr = sp.empty(out_shape, dtype=float)
                        out_arr[ii, jj, ...] = ma.var(Data.data, 0)
                sp.save(output_fname, out_arr)
            else :
                raise ValueError("Invalid noise model: " + model)


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Measure(str(sys.argv[1])).execute()
               
