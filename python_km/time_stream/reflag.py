"""Time stream module reprocesses data after comparing with map.

Right now just reflags data but will eventually do a fine calibration maybe?"""

import scipy as sp
import numpy.ma as ma

from core import fitsGBT
import base_single
import kiyopy.custom_exceptions as ce

params_init = {
               'thres' : 3.0,
               'subtracted_input_root' : './',
               'subtracted_output_root' : './subtracted_'
              }

prefix = 'sf_'

class ReFlag(base_single.BaseSingle) :

    prefix = prefix
    params_init = params_init
    
    def process_file(self, file_ind) :
        params = self.params
        file_middle = params['file_middles'][file_ind]
        input_fname = (params['input_root'] + file_middle +
                       params['input_end'])
        sub_input_fname = (params['subtracted_input_root'] + file_middle
                           + params['input_end'])
        output_fname = (params['output_root']
                        + file_middle + params['output_end'])
        sub_output_fname = (params['subtracted_output_root']
                            + file_middle + params['output_end'])
        Writer = fitsGBT.Writer(feedback=self.feedback)
        SubWriter = fitsGBT.Writer(feedback=self.feedback)
        
        # Read in the data, and loop over data blocks.
        Reader = fitsGBT.Reader(input_fname, feedback=self.feedback)
        SubReader = fitsGBT.Reader(sub_input_fname, feedback=self.feedback)
        if (sp.any(Reader.scan_set != SubReader.scan_set)
            or sp.any(Reader.IF_set != SubReader.IF_set)) :
            raise ce.DataError("IFs and scans don't match signal subtracted"
                               " data.")
        # Get the number of scans if asked for all of them.
        scan_inds = params['scans']
        if len(scan_inds) == 0 or scan_inds is None :
            scan_inds = range(len(Reader.scan_set))
        if_inds = params['IFs']
        if len(if_inds) == 0 or scan_inds is None :
            if_inds = range(len(Reader.IF_set))
        # Loop over scans and IFs
        for thisscan in scan_inds :
            for thisIF in if_inds :
                Data = Reader.read(thisscan, thisIF)
                SubData = SubReader.read(thisscan, thisIF)
                # Now do the flagging.
                flag(Data, SubData, params['thres'])
                Writer.add_data(Data)
                SubWriter.add_data(SubData)
        # Finally write the data back to file.
        Writer.write(output_fname)
        SubWriter.write(sub_output_fname)

    
def flag(Data, NoiseData, thres=3.0) :
    """Flags data for outliers using a signal subtracted data set.
    
    Flags outliers of in a time stream data by looking at a version of the data
    that has had the signal subtracted out of it.  Each frequency channel,
    polarization and cal state are treated separately.

    Parameters
    ----------
    Data : DataBlock Object
        Data to be flaged.  Upon exit, this object will have new flags.
    NoiseData : DataBlock Object
        Version of `Data` with the signal subtracted.
    thres : float
        Threshold for flagging in units of sigma (default is 3.0).
    """
    
    # Get the deviation from the mean.
    residuals = ma.anom(NoiseData.data, 0)
    # Get indices above the threshold.
    mask = abs(residuals) > thres*ma.std(NoiseData.data, 0)
    # Mask the data.
    Data.data[mask] = ma.masked
    NoiseData.data[mask] = ma.masked


# If this file is run from the command line, execute the main function.
if __name__=="__main__":
    import sys
    ReFlag(str(sys.argv[1])).execute()
