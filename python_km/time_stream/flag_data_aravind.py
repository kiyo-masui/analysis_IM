#! /usr/bin/python
"""This module flags rfi and other forms of bad data.
"""

import os

import numpy.ma as ma
import scipy as sp
from scipy import weave

import kiyopy.custom_exceptions as ce
import base_single
import an.rfi as rfi

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
                   'pol_thres' : 5
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


def apply_cuts(Data, sig_thres=5.0, pol_thres=5.0, width=2, flatten=True,
               der_flags=10, der_width=2) :
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

        # Allocate memory in SWIG array types.
        freq_array = rfi.new_doublearray(dims[3])
        data_array = rfi.new_doublearray(dims[3])
        cross_array = rfi.new_doublearray(dims[3])
        fit_array = rfi.new_doublearray(dims[3])
        mask_array = rfi.new_intarray(dims[3])
        for fii in xrange(dims[3]) :
            rfi.doublearray_setitem(freq_array,fii,freq[fii])
        
        # Outer loops performed in python.
        for tii in range(dims[0]) :
            for cjj in range(dims[2]) :
                # Polarization cross correlation coefficient.
                cross = ma.sum(data[tii,xy_inds,cjj,:]**2, 0)
                cross /= data[tii,xx_ind,cjj,:]*data[tii,yy_ind,cjj,:]
                cross = ma.filled(cross, 1000.)

                # Copy data to the SWIG arrays.
                # This may be confusing: Data is a DataBlock object which has
                # an attribute data which is a masked array.  Masked arrays
                # have attributes data (an array) and mask (a bool array).  So
                # thisdata and thismask are just normal arrays.
                for fkk in xrange(dims[3]) :
                    rfi.doublearray_setitem(data_array,fkk, 
                                            data.data[tii,xx_ind,cjj,fkk]);
                    rfi.doublearray_setitem(cross_array,fkk,cross[fkk]);
                rfi.get_fit(cross_array,freq_array,fit_array)

                #print (pol_thres, width, int(flatten), derivative_cut, 
                #       der_flags, der_width)
                rfi.clean(pol_thres, width, int(flatten), derivative_cut, 
                          der_flags, der_width, fit_array, cross_array,
                          data_array, freq_array, mask_array)
                # Copy mask the flagged points and set up for YY.
                for fkk in xrange(dims[3]) :
                    if rfi.intarray_getitem(mask_array, fkk) :
                        data[tii,:,cjj,fkk] = ma.masked
                    rfi.doublearray_setitem(data_array,fkk, 
                                            data.data[tii,yy_ind,cjj,fkk]);
                # Clean The YY pol.
                rfi.clean(pol_thres, width, int(flatten), derivative_cut, 
                          der_flags, der_width, fit_array, cross_array,
                          data_array, freq_array, mask_array)
                # Mask flagged YY data.
                for fkk in xrange(dims[3]) :
                    if rfi.intarray_getitem(mask_array, fkk) :
                        data[tii,:,cjj,fkk] = ma.masked
        
        # Not sure if the desctors already do this, but free up memory.
        rfi.delete_intarray(mask_array)
        rfi.delete_doublearray(freq_array)
        rfi.delete_doublearray(data_array)
        rfi.delete_doublearray(cross_array)
        rfi.delete_doublearray(fit_array)
    
    if sig_thres > 0 :
        # Will work with squared data so no square roots needed.
        nt = Data.dims[0]
        data = Data.data[:,[xx_ind, yy_ind], :, :]
        norm_data = (data/ma.mean(data, 0) - 1)**2
        var = ma.mean(norm_data, 0)
        # Use an iteratively flagged mean instead of a median as medians are
        # very computationally costly.
        for jj in range(3) :
            # First check for any outliers that could throw the initial mean off
            # more than sqrt(3*nt/4) sigma.  This is the 'weak' cut.
            bad_mask = ma.where(ma.logical_or(
                                norm_data[:,0,:,:] > 3.*nt*var[0,:,:]/4., 
                                norm_data[:,1,:,:] > 3.*nt*var[1,:,:]/4.))
            if len(bad_mask[0]) == 0:
                break
            for ii in range(4) :
                data_this_pol = Data.data[:,ii,:,:]
                data_this_pol[bad_mask] = ma.masked
            # Recalculate the mean and varience for the next iteration or step.
            norm_data = (data/ma.mean(data, 0) - 1)**2
            var = ma.mean(norm_data, 0)

        # Strong cut, flag data deviating from the t and f var.
        # Want the varience over t,f. mean_f(var_t(data)) = var_t,f(data)
        var = ma.mean(var, -1)
        bad_mask = ma.where(ma.logical_or(
                        norm_data[:,0,:,:] > var[0,:,sp.newaxis]*sig_thres**2, 
                        norm_data[:,1,:,:] > var[0,:,sp.newaxis]*sig_thres**2))
        for ii in range(4) :
            data_this_pol = Data.data[:,ii,:,:]
            data_this_pol[bad_mask] = ma.masked



# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    FlagData(str(sys.argv[1])).execute()

