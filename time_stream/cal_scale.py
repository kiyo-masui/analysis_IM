#!/usr/bin/python
"""Module that put data in units of cal temperture and subtracts median."""

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
import base_single


class CalScale(base_single.BaseSingle) :
    """Pipeline module that performs scales data by the cal on the data.
    
    See the main funciton of this module: cal_scale.scale_by_cal for a detailed
    doc string.
    """

    prefix = 'cs_'
    params_init = {
                   'scale_time_average' : True,
                   'scale_freq_average' : False,
                   'subtract_time_median' : False,
                   'scale_freq_average_mod' : False,
                   'rotate_to_cal_phase' : False,
                   'invert_cal' : False,
                   }

    def action(self, Data):
        scale_by_cal(Data, self.params['scale_time_average'],
                     self.params['scale_freq_average'], 
                     self.params['subtract_time_median'],
                     self.params['scale_freq_average_mod'],
                     self.params['rotate_to_cal_phase'],
                     self.params['invert_cal']
                     )
        Data.add_history('Converted to units of noise cal temperture.')
        return Data


def scale_by_cal(Data, scale_t_ave=True, scale_f_ave=False, sub_med=False,
                 scale_f_ave_mod=False, rotate=False,invert=False) :
    """Puts all data in units of the cal temperature.
    
    Data is put into units of the cal temperature, thus removing dependence on
    the gain.  This can be done by dividing by the time average of the cal
    (scale_t_ave=True, Default) thus removing dependence on the frequency-
    dependant gain.  Alternatively, you can scale by the frequency average to
    remove the time-dependent gain (scale_f_ave=True). Data is then in units of
    the frequency averaged cal temperture. You can also do both (recommended).
    After some scaling the data ends up in units of the cal temperture as a
    funciton of frequency.

    Optionally you can also subtract the time average of the data off here
    (subtract_time_median), since you might be done with the cal information at
    this point.
    """
    
    on_ind = 0
    off_ind = 1
    if invert:
        on_ind = 1
        off_ind = 0
    else: 
        if (Data.field['CAL'][on_ind] != 'T' or
            Data.field['CAL'][off_ind] != 'F') :
                raise ce.DataError('Cal states not in expected order.')
    
    if tuple(Data.field['CRVAL4']) == (-5, -7, -8, -6) :
        # Here we check the polarizations and cal indicies
        xx_ind = 0
        yy_ind = 3
        xy_inds = [1,2]
        
        # A bunch of calculations used to test phase closure.  Not acctually
        # relevant to what is being done here.
        #a = (Data.data[5, xy_inds, on_ind, 15:20]
        #     - Data.data[5, xy_inds, off_ind, 15:20])
        #a /= sp.sqrt( Data.data[5, xx_ind, on_ind, 15:20] 
        #              - Data.data[5, xx_ind, off_ind, 15:20])
        #a /= sp.sqrt( Data.data[5, yy_ind, on_ind, 15:20] 
        #              - Data.data[5, yy_ind, off_ind, 15:20])
        #print a[0,:]**2 + a[1,:]**2
        
        diff_xx = Data.data[:,xx_ind,on_ind,:] - Data.data[:,xx_ind,off_ind,:]
        diff_yy = Data.data[:,yy_ind,on_ind,:] - Data.data[:,yy_ind,off_ind,:]
        
        if scale_t_ave :
            # Find the cal means (in time) and scale by them.
            # Means work much better than medians.  Medians seems to bias the
            # result by up to 10%.  This seems to be discretization noise.  Cal
            # switches fast enough that we shouldn't need this anyway.
            cal_tmed_xx = ma.mean(diff_xx, 0)
            cal_tmed_yy = ma.mean(diff_yy, 0)
            cal_tmed_xx[sp.logical_or(cal_tmed_xx<=0, cal_tmed_yy<=0)] = ma.masked
            cal_tmed_yy[cal_tmed_xx.mask] = ma.masked

            Data.data[:,xx_ind,:,:] /= cal_tmed_xx
            Data.data[:,yy_ind,:,:] /= cal_tmed_yy
            Data.data[:,xy_inds,:,:] /= ma.sqrt(cal_tmed_yy*cal_tmed_xx)

        if scale_f_ave :
            # The frequency gains have have systematic structure to them, 
            # they are not by any approximation gaussian distributed.  Use
            # means, not medians across frequency.
            operation = ma.mean
            cal_fmea_xx = operation(diff_xx, -1)
            cal_fmea_yy = operation(diff_yy, -1)
            
            # Flag data with wierd cal power.  Still Experimental.
            cal_fmea_xx[sp.logical_or(cal_fmea_xx<=0,cal_fmea_yy<=0)] = ma.masked
            cal_fmea_yy[cal_fmea_xx.mask] = ma.masked
            cal_xx = ma.mean(cal_fmea_xx)
            cal_yy = ma.mean(cal_fmea_yy)
            cal_fmea_xx[sp.logical_or(abs(cal_fmea_xx.anom()) >= 0.1*cal_xx,
                            abs(cal_fmea_yy.anom()) >= 0.1*cal_yy)] = ma.masked
            cal_fmea_yy[cal_fmea_xx.mask] = ma.masked
            
            ntime = len(cal_fmea_xx)
            cal_fmea_xx.shape = (ntime, 1, 1)
            cal_fmea_yy.shape = (ntime, 1, 1)
            Data.data[:,xx_ind,:,:] /= cal_fmea_xx
            Data.data[:,yy_ind,:,:] /= cal_fmea_yy
            cal_fmea_xx.shape = (ntime, 1, 1, 1)
            cal_fmea_yy.shape = (ntime, 1, 1, 1)
            Data.data[:,xy_inds,:,:] /= ma.sqrt(cal_fmea_yy*cal_fmea_xx)

        if scale_f_ave_mod :
            # The frequency gains have have systematic structure to them, 
            # they are not by any approximation gaussian distributed.  Use
            # means, not medians across frequency.
            operation = ma.mean
            cal_fmea_xx = operation(diff_xx, -1)
            cal_fmea_yy = operation(diff_yy, -1)
            cal_fmea_xx_off = operation(Data.data[:,xx_ind,off_ind,:], -1)
            cal_fmea_yy_off = operation(Data.data[:,yy_ind,off_ind,:], -1)


            sys_xx = cal_fmea_xx_off/cal_fmea_xx
            sys_yy = cal_fmea_yy_off/cal_fmea_yy
            percent_ok = 0.03
            sys_xx_tmed = ma.median(sys_xx)
            sys_yy_tmed = ma.median(sys_yy)

            maskbad_xx = (sys_xx > sys_xx_tmed + sys_xx_tmed*percent_ok)|(sys_xx < sys_xx_tmed - sys_xx_tmed*percent_ok)
            maskbad_yy = (sys_yy > sys_yy_tmed + sys_yy_tmed*percent_ok)|(sys_yy < sys_yy_tmed - sys_yy_tmed*percent_ok)

            cal_fmea_xx[sp.logical_or(cal_fmea_xx<=0,cal_fmea_yy<=0)] = ma.masked
            cal_fmea_yy[cal_fmea_xx.mask] = ma.masked
            cal_fmea_xx[maskbad_xx] = ma.masked
            cal_fmea_yy[maskbad_yy] = ma.masked  
            cal_xx = ma.mean(cal_fmea_xx)
            cal_yy = ma.mean(cal_fmea_yy)

            ntime = len(cal_fmea_xx)
            cal_fmea_xx.shape = (ntime, 1, 1)
            cal_fmea_yy.shape = (ntime, 1, 1)
            Data.data[:,xx_ind,:,:] /= cal_fmea_xx
            Data.data[:,yy_ind,:,:] /= cal_fmea_yy
            cal_fmea_xx.shape = (ntime, 1, 1, 1)
            cal_fmea_yy.shape = (ntime, 1, 1, 1)
            Data.data[:,xy_inds,:,:] /= ma.sqrt(cal_fmea_yy*cal_fmea_xx) 

        if scale_f_ave and scale_t_ave :
            # We have devided out t_cal twice so we need to put one factor back
            # in.
            cal_xx = operation(cal_tmed_xx)
            cal_yy = operation(cal_tmed_yy)
            Data.data[:,xx_ind,:,:] *= cal_xx
            Data.data[:,yy_ind,:,:] *= cal_yy
            Data.data[:,xy_inds,:,:] *= ma.sqrt(cal_yy*cal_xx)

        if scale_f_ave_mod and scale_t_ave :
            #Same divide out twice problem.
            cal_xx = operation(cal_tmed_xx)
            cal_yy = operation(cal_tmed_yy)
            Data.data[:,xx_ind,:,:] *= cal_xxcal_imag_mean
            Data.data[:,yy_ind,:,:] *= cal_yy
            Data.data[:,xy_inds,:,:] *= ma.sqrt(cal_yy*cal_xx)
           
        if scale_f_ave and scale_f_ave_mod :
            raise ce.DataError("time averaging twice") 

        if rotate:
            # Define the differential cal phase to be zero and rotate all data
            # such that this is true.
            cal_real_mean = ma.mean(Data.data[:,1,0,:] - Data.data[:,1,1,:], 0)
            cal_imag_mean = ma.mean(Data.data[:,2,0,:] - Data.data[:,2,1,:], 0)
            # Get the cal phase angle as a function of frequency.
            cal_phase = -ma.arctan2(cal_imag_mean, cal_real_mean)

            # Rotate such that the cal phase is zero. Imperative to have a
            # temporary variable.
            New_data_real = (ma.cos(cal_phase) * Data.data[:,1,:,:]
                             - ma.sin(cal_phase) * Data.data[:,2,:,:])
            New_data_imag = (ma.sin(cal_phase) * Data.data[:,1,:,:]
                             + ma.cos(cal_phase) * Data.data[:,2,:,:])
            Data.data[:,1,:,:] = New_data_real
            Data.data[:,2,:,:] = New_data_imag

    elif tuple(Data.field['CRVAL4']) == (1, 2, 3, 4) :
        # For the shot term, just devide everything by on-off in I.
        I_ind = 0
        cal_I_t = Data.data[:,I_ind,on_ind,:] - Data.data[:,I_ind,off_ind,:]
        cal_I = ma.mean(cal_I_t, 0)

        Data.data /= cal_I
    else :
        raise ce.DataError("Unsupported polarization states.")

    # Subtract the time median if desired.
    if sub_med :
        Data.data -= ma.median(Data.data, 0)

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    CalScale(str(sys.argv[1])).execute()

