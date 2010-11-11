#! /usr/bin/python
"""This module flags rfi and other forms of bad data.
"""

import os

import numpy.ma as ma
import scipy as sp
from scipy import weave

import kiyopy.custom_exceptions as ce
import base_single

class FlagData(base_single.BaseSingle) :
    """Pipeline module that flags rfi and other forms of bad data.

    For lots of information look at the doc-string for 
    flag_data_aravind.apply_cuts.
    """

    prefix = 'fd_'
    params_init = {
                   # In multiples of the standard deviation of the whole block
                   # once normalized to the time median.
                   'sigma_thres' : 5,
                   # In multiples of the measured standard deviation.
                   'pol_thres' : 3
                   }
    feedback_title = 'New flags each data Block: '
    def action(self, Data) :
        already_flagged = ma.count_masked(Data.data)
        apply_cuts(Data, self.params['sigma_thres'], self.params['pol_thres'])
        new_flags = ma.count_masked(Data.data) - already_flagged
        self.block_feedback = str(new_flags) + ', '

        Data.add_history('Flagged Bad Data.', ('Sigma threshold: ' +
                    str(self.params['sigma_thres']), 'Polarization threshold: '
                    + str(self.params['pol_thres'])))
        return Data


def apply_cuts(Data, sig_thres=5, pol_thres=3.0, width=0, flatten=True,
               der_flags=10, der_width=0) :
    """Flags bad data from RFI and far outliers..
    
    Masks any data that is further than 'sig_thres' sigmas from the mean of the
    entire data block (all frequencies, all times).  Each cal, fequency and
    polarization is first normalized to the time mean.  Cut only checked for XX
    and YY polarizations but all polarizations are masked.  This cut is
    deactivated by setting 'sig_thres' < 0.

    Aravind rfi cut is also performed, cutting data with polarization lines, as
    well as lines in XX and YY.

    Arguments:
        Data -      A DataBlock object containing data to be cleaned.
        sig_thres - (float) Any XX or YY data that deviates by this more than 
                    this many sigmas is flagged.  Data first normalized to 
                    time median.
        pol_thres - (float) Any data cross polarized by more than this 
                    many sigmas is flagged.
        width -     (int) In the polarization cut, flag data within this many
                    frequency bins of offending data.
        flatten -   (Bool) Preflatten the polarization spectrum.
        der_flags - (int) Find RFI in XX and YY by looking for spikes in the
                    derivative.  Flag this many spikes.
        der_width - (int) Same as width but for the derivative cut.
    """
    
    # Here we check the polarizations and cal indicies
    xx_ind = 0
    yy_ind = 3
    xy_inds = [1,2]
    if (Data.field['CRVAL4'][xx_ind] != -5 or
        Data.field['CRVAL4'][yy_ind] != -6 or
        Data.field['CRVAL4'][xy_inds[0]] != -7 or
        Data.field['CRVAL4'][xy_inds[1]] != -8) :
            raise ce.DataError('Polarization types not as expected,'
                               ' function needs to be generalized.')
    on_ind = 0
    off_ind = 1
    if (Data.field['CAL'][on_ind] != 'T' or
        Data.field['CAL'][off_ind] != 'F') :
            raise ce.DataError('Cal states not in expected order.')
    
    if pol_thres > 0 :
        Data.calc_freq()
        freq = Data.freq
        dims = Data.dims
        data = Data.data

        if dims[3] != 2048 :
            raise ce.DataError('C code expects 2048 frequency channels.')
        
        derivative_cut = 0
        if der_flags > 0:
            derivative_cut = 1

        # Cast variables to a definate type for C code.
        pol_thres = float(pol_thres)
        width = int(width)
        flatten = int(flatten)
        der_flags = int(der_flags)
        der_width = int(der_width)
        
        # Get Aravind's rfi C code.
        this_dir = os.path.dirname(__file__)
        flag_code = open(this_dir + '/rfi_flag.c', 'r').read()
        # Driver code for Aravind's code.  Note that the C code is only
        # recompiled if the md5sum of the below string changes.  If you change
        # flag_code, you have to force a recompile.
        code = """
        #line 106 "flag_data_aravind.py"
        clean(pol_thres, width, flatten, derivative_cut, der_flags, 
              der_width, cross, thisdata, freq, thismask);
        """
        # Variables that should be carried from python to C.
        variables = ['pol_thres', 'width', 'flatten', 'derivative_cut', 
                     'der_flags', 'der_width', 'freq', 'thisdata', 'thismask', 
                     'cross']
        
        # Outer loops performed in python.
        for tii in range(dims[0]) :
            for cjj in range(dims[2]) :
                
                # Polarization cross correlation coefficient.
                cross = ma.sum(data[tii,xy_inds,cjj,:]**2, 0)
                cross /= data[tii,xx_ind,cjj,:]*data[tii,yy_ind,cjj,:]
                cross = ma.filled(cross, 100.)
                print sp.where(cross),
                
                # XX then YY
                # This may be confusing: Data is a DataBlock object which has
                # an attribute data which is a masked array.  Masked arrays
                # have attributes data (an array) and mask (a bool array).  So
                # thisdata and thismask are just normal arrays.
                thisdata = data.data[tii,xx_ind,cjj,:]
                thismask = sp.array(data.mask[tii,xx_ind,cjj,:], dtype=sp.int32)
                print sp.where(thismask),
                # Compile (if required) and execute C code.
                weave.inline(code, variables, support_code=flag_code)
                data[tii,:,cjj,thismask] = ma.masked
                print sp.where(thismask)

                thisdata = data.data[tii,yy_ind,cjj,:]
                thismask = sp.array(data.mask[tii,yy_ind,cjj,:], dtype=sp.int32)
                weave.inline(code, variables, support_code=flag_code)
                data[tii,:,cjj,thismask] = ma.masked

    if sig_thres > 0 :
        # Will work with squared data so no square roots needed.
        data = Data.data[:,[xx_ind, yy_ind], :, :]
        # Median is crucial, don't want to mess up a whole channel normalization
        # with one bad data point.
        norm_data = (data/ma.median(data, 0) - 1)**2
        var = ma.mean(norm_data)

        bad_mask = ma.where(ma.logical_or(norm_data[:,0,:,:] > var*sig_thres**2, 
                                          norm_data[:,1,:,:] > var*sig_thres**2))
        for ii in range(4) :
            data_this_pol = Data.data[:,ii,:,:]
            data_this_pol[bad_mask] = ma.masked

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    FlagData(str(sys.argv[1])).execute()

