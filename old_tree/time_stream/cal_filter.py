#!/usr/bin/python
"""Module filters the cal time stream and applys filter to data."""

import scipy as sp
import numpy.ma as ma
from scipy import signal
import matplotlib.pyplot as plt

import kiyopy.custom_exceptions as ce
import base_single


class CalFilter(base_single.BaseSingle) :
    """Module filters the noise cal data and then scales the data.
    
    The idea here is to elimiate achromatic gain variations above a given time
    scale.  The band average noise cal data is low pass filtered to give an
    estimate of the achromatic gain variation.  The data is then divided by
    this extimate to get rid of the variations.
    
    Since the estimate is noisy, it is certainly possible to acctually 
    introduce gain variations this way, so it's all a question of time scales
    and where you want these variations to be.  Any variations you introduce
    will be perfectly achromatic, so it shouldn't kill any signal.
    """

    prefix = 'cf_'
    params_init = {'filter_type' : 'rectangular',
                   'filter_size' : 10.,
                   # Options are 'seconds' and 'bins'.
                   'filter_size_units' : 'bins'
                   }

    def action(self, Data):
        filter_cal_scale(Data, self.params['filter_size'],
                     filter_type=self.params['filter_type'],
                     filter_size_units=self.params['filter_size_units'])
        Data.add_history('Scaled by achromatic filtered cal.')
        return Data


def filter_cal_scale(Data, size, filter_type='rectangular',
                     filter_size_units='bins') :
    """Estimates Achromatic gain fluctuations and scales by them.
    """

    on_ind = 0
    off_ind = 1
    n_time = Data.data.shape[0]
    if (Data.field['CAL'][on_ind] != 'T' or
        Data.field['CAL'][off_ind] != 'F') :
            raise ce.DataError('Cal states not in expected order.')
    
    if tuple(Data.field['CRVAL4']) == (-5, -7, -8, -6) :
        # Here we check the polarizations and cal indicies
        xx_ind = 0
        yy_ind = 3
        xy_inds = [1,2]
        
        diff_xx = Data.data[:,xx_ind,on_ind,:] - Data.data[:,xx_ind,off_ind,:]
        diff_yy = Data.data[:,yy_ind,on_ind,:] - Data.data[:,yy_ind,off_ind,:]
        cal_xx = ma.mean(diff_xx, -1)
        cal_yy = ma.mean(diff_yy, -1)
        mask = ma.getmaskarray(Data.data)
        #mask = sp.any(sp.any(mask, 1), 1)
        # XXX Wrong, Just needs to be factorizable.  Need to devellop and
        # algorithem for ensuring this, but for now, just use this. This
        # shouldn't be too bad for the current RFI flagging algorithm (Apr
        # 2012).
        time_mask = sp.any(sp.any(sp.any(mask, 1), 1), 1)
        # Sanity check the masks.
        if not (sp.all(time_mask[ma.getmaskarray(cal_xx)])
                and sp.all(time_mask[ma.getmaskarray(cal_yy)])):
            msg = "Doesn't make since, this should always be true."
            raise RuntimeError(msg)
        # Convert to normal arrays.
        cal_xx = cal_xx.filled(0)
        cal_yy = cal_yy.filled(0)
        # XXX
        #Data.calc_time()
        #time_unmask = sp.logical_not(time_mask)
        #plt.plot(Data.time[time_unmask], cal_xx[time_unmask])
        # Now set up the filter.
        if filter_type == 'rectangular':
            if filter_size_units=='bins':
                n_bins_filter = int(size)
                if n_bins_filter % 2 == 0:
                    raise ValueError("Rectangular filter should have an odd" 
                                     "number of bins.")
            elif filter_size_units == 'seconds':
                Data.calc_time()
                dt = abs(sp.mean(sp.diff(Data.time)))
                n_bins_filter = size / dt
                # Round to the nearest odd number.
                n_bins_filter = 2 * int(round(n_bins_filter / 2. - 0.5)) + 1
            else:
                msg = 'Filter unit type unsupported.'
                raise ValueError(msg)
            kernal = sp.ones(n_bins_filter) / n_bins_filter
        else:
            msg = 'Filter type unsupported.'
            raise ValueError(msg)
        # Now that we know the kernal size, figure out what elements will
        # be newly masked by the smoothing.
        half_width = n_bins_filter // 2
        old_mask = time_mask.copy()
        for ii in range(n_time):
            if old_mask[ii]:
                time_mask[ii - half_width:ii + half_width + 1] = True
        # Also mask the edges.
        time_mask[:half_width] = True
        time_mask[-half_width:] = True
        # Now acctually do the convolution.
        cal_xx = signal.convolve(cal_xx, kernal, mode='same')
        cal_yy = signal.convolve(cal_yy, kernal, mode='same')
        # XXX
        #Data.calc_time()
        #time_unmask = sp.logical_not(time_mask)
        #plt.plot(Data.time[time_unmask], cal_xx[time_unmask])
        #plt.plot(Data.time, time_mask)
        plt.show()
        # Replace invalid entries with unity (They get masked later anyway).
        cal_xx[time_mask] = 1.
        cal_yy[time_mask] = 1.
        # Calibrate and apply mask.
        Data.data[:,xx_ind,:,:] /= cal_xx[:,None,None]
        Data.data[:,yy_ind,:,:] /= cal_yy[:,None,None]
        cross_cal = sp.sqrt(cal_xx * cal_yy)
        Data.data[:,xy_inds,:,:] /= cross_cal[:,None,None,None]
        # Apply the mask.
        Data.data[time_mask,...] = ma.masked
    #elif tuple(Data.field['CRVAL4']) == (1, 2, 3, 4) :
    else :
        raise ce.DataError("Unsupported polarization states.")

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    CalScale(str(sys.argv[1])).execute()

