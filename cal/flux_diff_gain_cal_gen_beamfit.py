"""Procedure to calculate only the flux and differential gain between XX and YY for each frequency from on-off scans of a calibrator such as 3C286.

Currently setup so that it uses all the onoff scans from a session for a particular calibrator and gets the best fit. 
Run in analysis_IM: python cal/flux_diff_gain_cal_gen.py input/tcv/diff_gain_gen_guppi.ini
Note that the .ini file should indicate which session(s) and sourse you want to use. Script is run using data from a single source. The output is saved in my data directory under the folder diff_gain_params as a .txt file with three columns (freq, XXGain, YYGain).
Need to change script to use correct source data when changing sources.
Also, script currently also includes rotation measure correction to the rotation matrix (which only matters when generating using a polarized source such as 3C286).
 """
import os

from scipy.optimize import *
import scipy as sp
import numpy.ma as ma
import numpy as np

from utils import misc
import source
import matplotlib.pyplot as plt
from scipy import optimize, linalg

class FormattedData(object):
    """All data that goes into a fit.

    This class holds GBT SDfits data in a format that makes it easy for
    fitting.
    
    Parameters
    ----------
    Blocks : list of DataBlock objects.
        Data to be used for fit.
    """

    def __init__(self, Blocks):
        self._data = self.preprocess_data(Blocks)
        # Figure out all dimensions.
        self.n_chan = self._data[0]['data'].shape[0]
        self.n_pol = self._data[0]['data'].shape[1]
        n_time = 0
        for data in self._data:
            s = data['data'].shape
            n_time += s[2]
            if s[0] != self.n_chan or s[1] != self.n_pol:
                msg = "Data blocks have conflicting dimensions."
                raise ValueError(msg)
        self.n_time = n_time
        self.n_scan = len(self._data)

    def get_data_weight_chan(self, chan_ind):
        """Get all data and weights for a single channel."""

        data_arr = np.empty((self.n_pol, self.n_time), dtype=np.float64)
        weights_arr = np.empty((self.n_pol, self.n_time), dtype=np.float64)
        for data, start, end in self.iterate_scans():
            this_chan_data = data['data'][chan_ind,:,:]
            this_chan_weights = data['weight'][chan_ind,:,:]
            data_arr[:,start:end] = this_chan_data
            weights_arr[:,start:end] = this_chan_weights
        return data_arr, weights_arr

    def get_all_data_weight(self):
        data_arr = np.empty((self.n_chan, self.n_pol, self.n_time),
                            dtype=np.float64)
        weights_arr = np.empty((self.n_chan, self.n_pol, self.n_time),
                               dtype=np.float64)
        for data, start, end in self.iterate_scans():
            this_chan_data = data['data']
            this_chan_weights = data['weight']
            data_arr[:,:,start:end] = this_chan_data
            weights_arr[:,:,start:end] = this_chan_weights
        return data_arr, weights_arr

    def get_all_field(self, field):
        """Get all of a specific data field."""

        format = self._data[0][field].dtype
        out = np.empty(self.n_time, dtype=format)
        for data, start, end in self.iterate_scans():
            this_field = data[field]
            out[start:end] = this_field
        return out

    def generate_scan_basis_polys(self, order):
        """Generate basis time polynomials for each scan up to order.
        
        Generate basis polynomials as a function of time for each scan 
        up to the passed maximum order.

        Parameters
        ----------
        order : integer
            Maximum order of the polynomials.

        Returns
        -------
        polys : ndarray of shape (n_scan, `order`, n_time)
        """

        polys = np.zeros((self.n_scan, order, self.n_time),
                         dtype=np.float64)
        scan_ind = 0
        for data, start, end in self.iterate_scans():
            time = data['time']
            this_poly = misc.ortho_poly(time, order)
            polys[scan_ind,:,start:end] = this_poly
            scan_ind +=1
        return polys

    def iterate_scans(self):
        """Returns an iterator over the data scans, with start and end time
        indeces."""

        class iterator(object):

            def __init__(self_itr):
                self_itr.scan = 0
                self_itr.time_ind = 0
                self_itr._data = self._data

            def __iter__(self_itr):
                return self_itr

            def next(self_itr):
                if len(self_itr._data) <= self_itr.scan:
                    raise StopIteration()
                start = self_itr.time_ind
                data = self._data[self_itr.scan]
                self_itr.scan += 1
                end = start + data['data'].shape[2]
                self_itr.time_ind = end
                return data, start, end

        return iterator()

    def preprocess_data(self, Blocks):
        """Puts data into convienient format and calculates weights."""

        # Where we will store formated data.
        data_list = []
        for Data in Blocks:
            # First do some basic checks on data format.
            if not tuple(Data.field['CRVAL4']) == (-5, -7, -8, -6):
                msg = "Expected data polarizations to be XX, XY, YX, YY."
                raise NotImplementedError(msg)
            if not Data.data.shape[2] == 1:
                msg = "Multiple cal states, ignoring all but first."
                raise Warning(msg)

            # Rearrange data.
            data = np.swapaxes(Data.data[:,:,0,:], 0, 2)
            mask = ma.getmaskarray(data).copy()
            data = data.filled(0)
            # XX and YY data should all be positive.
            mask[data[:,[0],:] < 0] = True
            mask[data[:,[3],:] < 0] = True
            # If any polarizations are masked, all polarizations should be.
            mask[np.any(mask, 1)[:,None,:]] = True
            data[mask] = 0

            # Make the noise weights.
            weight = data.copy()
            cross_norm = np.sqrt(weight[:,0,:] * weight[:,3,:])
            weight[:,[1,2],:] = cross_norm[:,None,:]
            weight[mask] = 1
            weight = 1 / weight
            weight[mask] = 0

            # Get the pointing.
            Data.calc_pointing()
            ra = Data.ra
            dec = Data.dec
            # Telescope frame.
            az = Data.field['CRVAL2']
            el = Data.field['CRVAL3']
            UT = Data.field['DATE-OBS']
            # Get the time array.
            Data.calc_time()
            time = Data.time

            # Store it all and add to list.
            this_data = {}
            this_data['data'] = data
            this_data['weight'] = weight
            this_data['ra'] = ra
            this_data['dec'] = dec
            this_data['az'] = az
            this_data['el'] = el
            this_data['time'] = time
            this_data['UT'] = UT
            data_list.append(this_data)
        return data_list

def calcGain(OnData,OffData,file_num,freq_len,src,beamwidth):
    """ Perform gain calibration on dataset.

    """
    def peval(p,data):
        d = data
        XG = p[0]
        YG = p[1]
        act = sp.zeros(len(data)*4)
        for i in range(0,len(act),4):
            act[i] = XG*d[i/4,0]
            act[i+1] = 0
            act[i+2] = 0
            act[i+3] =YG*d[(i+3)/4,3]
        return act

    def residuals(p,errors,freq_val,src,theta,data,width,file_num):

        wavelength = 300.0/freq_val
        BW = width*180./sp.pi
        JtoK = (sp.pi*wavelength**2)/(8*1380.648*BW**2)
	Jsrc_name = ['3C286','3C48','3C67','3C147','3C295']
	Jsrc_val = [19.74748409*pow((750.0/freq_val),0.49899785),
		    25.15445092*pow((750.0/freq_val),0.75578842),
		    4.56303633*pow((750.0/freq_val),0.59237327),
		    31.32846821*pow((750.0/freq_val),0.52113534),
                    34.11187767*pow((750.0/freq_val),0.62009421)]
	for i in range(0,len(Jsrc_name)):
	    if Jsrc_name[i]==src:
		src_ind = i
	PAsrc = [33.*sp.pi/180.,0.,0.,0.,0.,0.]
	Psrc = [0.07,0.,0.,0.,0.]
	Isrc = Jsrc_val[src_ind]*JtoK
	Qsrc = Isrc*Psrc[src_ind]*sp.cos(2*PAsrc[src_ind])
	Usrc = Isrc*Psrc[src_ind]*sp.sin(2*PAsrc[src_ind])
        Vsrc = 0
        XXsrc0 = Isrc-Qsrc
        YYsrc0 = Isrc+Qsrc
        expec =sp.zeros(4*file_num)
        for i in range(0,len(source),4):
            expec[i] = (0.5*(1+sp.cos(2*theta[i]))*XXsrc0-sp.sin(2*theta[i])*Usrc+0.5*(1-sp.cos(2*theta[i]))*YYsrc0)
            expec[i+1] = 0
            expec[i+2] = 0
            expec[i+3] = (0.5*(1-sp.cos(2*theta[i]))*XXsrc0+sp.sin(2*theta[i])*Usrc+0.5*(1+sp.cos(2*theta[i]))*YYsrc0)
        err = (expec-peval(p,data))/errors
        return err
    
###################################################
# Setting labels for indices for later
    XX_ind = 0
    YY_ind = 3
    XY_ind = 1
    YX_ind = 2
 
    S_med_src = sp.zeros((file_num,4,freq_len))
    S_med = sp.zeros((file_num,4,freq_len))

    PA_on = []
    m=0
    for Data in OnBlocks:
        S_med_src[m,0,:] = ma.median(Data.data[:,XX_ind,:],axis=0)
        S_med_src[m,1,:] = ma.median(Data.data[:,XY_ind,:],axis=0)
        S_med_src[m,2,:] = ma.median(Data.data[:,YX_ind,:],axis=0)
        S_med_src[m,3,:] = ma.median(Data.data[:,YY_ind,:],axis=0)
        Data.calc_PA()
        for i in range(0,4):
  	    PA_on.append(ma.mean(Data.PA))
	Data.calc_freq()
	freq_val = Data.freq/1e6
        m+=1

    PA_off = []
    m=0
    for Data in OffBlocks:
        S_med[m,0,:] = ma.median(Data.data[:,XX_ind,:],axis=0)
        S_med[m,1,:] = ma.median(Data.data[:,XY_ind,:],axis=0)
        S_med[m,2,:] = ma.median(Data.data[:,YX_ind,:],axis=0)
        S_med[m,3,:] = ma.median(Data.data[:,YY_ind,:],axis=0)
        Data.calc_PA()
        for i in range(0,4):
	    PA_off.append(ma.mean(Data.PA))
        m+=1
 
    S_data = sp.zeros((file_num,4,freq_len))
    for i in range(0,len(S_med)):
        S_data[i,0,:] = S_med_src[i,0,:]-S_med[i,0,:]
        S_data[i,1,:] = S_med_src[i,1,:]-S_med[i,1,:]
	S_data[i,2,:] = S_med_src[i,2,:]-S_med[i,2,:]
	S_data[i,3,:] = S_med_src[i,3,:]-S_med[i,3,:]
#There are 2 parameters for this version p[0] is XX gain and p[1] is YY gain. 
    p0 = [1,1] # guessed preliminary values
    error = sp.ones(4*file_num)
    #Note that error can be used to weight the equations if not all set to one.

    p_val_out = sp.zeros((freq_len, 3))
    for f in range(0,freq_len):   
        plsq = leastsq(residuals,p0,args=(error,freq_val[f],src,PA_on,S_data[:,:,f],beamwidth[f],file_num),full_output=0, maxfev=5000)
        pval = plsq[0] # this is the 1-d array of results0

        p_val_out[f,0] = freq_val[f]
        p_val_out[f,1] = pval[0]
        p_val_out[f,2] = pval[1]

#    out_path = output_root+sess+'_diff_gain_calc'+output_end
#    np.savetxt(out_path,p_val_out,delimiter = ' ')
    return p_val_out

#If this file is run from the command line, execute the main function.
#if __name__ == "__main__":
#    import sys
#    DiffGainGen(str(sys.argv[1])).execute()

