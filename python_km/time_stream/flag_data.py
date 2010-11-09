#! /usr/bin/python
"""This module flags rfi and other forms of bad data.
"""

import numpy.ma as ma

import kiyopy.custom_exceptions as ce
import base_single

class FlagData(base_single.BaseSingle) :
    """Pipeline module that flags rfi and other forms of bad data.

    For lots of information look at the doc-string for flag_data.apply_cuts.
    """

    prefix = 'fd_'
    params_init = {
                   # In multiples of the standard deviation of the whole block
                   # once normalized to the time median.
                   'sigma_thres' : 5,
                   # In multiples of the *theoredical* standard deviation
                   'pol_thres' : 15
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


def apply_cuts(Data, sig_thres=5, pol_thres=15) :
    """Flags bad data from RFI and far outliers..
    
    Masks any data that is further than 'sig_thres' sigmas from the mean of the
    entire data block (all frequencies, all times).  Each cal, fequency and
    polarization is first normalized to the time mean.  Cut only checked for XX
    and YY polarizations but all polarizations are masked.  This cut is
    deactivated by setting 'sig_thres' < 0.

    Also masks data that is more polarized than expected assuming it's RFI.
    the theoredical expectations value of the cross correlation coefficient
    (q^2) (in the absence of a polarized signal) is 1/dnudt the standard
    deviation is ~3/dnudt. In practice the expectation value can be up to 20
    times too high do to leakages.  This cut cuts any data with q^2 >
    'pol_thres'*3/dnudt.  It is reccomended that 'pol_thres' be set in the 10-20
    range.  This cut is deactivated by setting 'pol_thres' < 0.
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
        # Stdev of cross correlation known from theory.
        cutoff = pol_thres*3/abs(Data.field['CDELT1']*Data.field['EXPOSURE'])
        # Deal with cal off data first.
        data = Data.data
        cross = ma.sum(data[:,xy_inds,off_ind,:]**2, 1)
        cross /= data[:,xx_ind,off_ind,:]*data[:,yy_ind,off_ind,:]
        tinds, finds = ma.where(cross > cutoff)
        data[tinds,:,off_ind,finds] = ma.masked
        
        # In the cal on data, the cal is polarized, so need to subtract the
        # time mean of cal_on - cal_off out of the cross polarizations first.

        cal = ma.median(data[:,xy_inds,on_ind,:] - data[:,xy_inds,off_ind,:],
                        0) 
        cross = ma.sum((data[:,xy_inds,on_ind,:] - cal)**2, 1)
        cross /= data[:,xx_ind,on_ind,:]*data[:,yy_ind,on_ind,:]
        tinds, finds = ma.where(cross > cutoff)
        data[tinds,:,on_ind,finds] = ma.masked

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

