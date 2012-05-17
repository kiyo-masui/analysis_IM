"""Time stream module reprocesses data after comparing with map.

Right now just reflags data but will eventually do a fine calibration maybe?"""

import scipy as sp
import numpy.ma as ma

from utils import misc
from kiyopy import utils
from core import fitsGBT
import base_single
import kiyopy.custom_exceptions as ce

import matplotlib.pyplot as plt

params_init = {
               'thres' : 3.,
               'max_noise_factor' : 3.,
               'smooth_modes_subtract' : 1,
               'filter_type' : 'edge',
               'subtracted_input_root' : './testdata/',
               'subtracted_output_root' : './subtracted_'
              }

prefix = 'sf_'

class ReFlag(base_single.BaseSingle) :

    prefix = prefix
    params_init = params_init

    def execute(self, n_processes=1) :
        utils.mkparents(self.params['subtracted_output_root'])
        base_single.BaseSingle.execute(self, n_processes)
    
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
        if self.feedback > 1 :
            print "New flags each block:",
        # Loop over scans and IFs
        for thisscan in scan_inds :
            for thisIF in if_inds :
                Data = Reader.read(thisscan, thisIF)
                SubData = SubReader.read(thisscan, thisIF)
                # Make sure they have agreeing masks to start.
                SubData.data[ma.getmaskarray(Data.data)] = ma.masked
                Data.data[ma.getmaskarray(SubData.data)] = ma.masked
                # Get initial number of flags.
                n_flags = ma.count_masked(Data.data)
                # Now do the flagging.
                flag(Data, SubData, params['thres'],
                     params['max_noise_factor'], 
                     params['smooth_modes_subtract'],
                     params['filter_type'])
                Data.add_history("Reflaged for outliers.", ("Used file: "
                    + utils.abbreviate_file_path(sub_input_fname),))
                SubData.add_history("Reflaged for outliers.")
                Writer.add_data(Data)
                SubWriter.add_data(SubData)
                # Report the number of new flags.
                n_flags = ma.count_masked(Data.data) - n_flags
                if self.feedback > 1 :
                    print n_flags,
        if self.feedback > 1 :
            print ''
        # Finally write the data back to file.
        utils.mkparents(output_fname)
        utils.mkparents(sub_output_fname)
        Writer.write(output_fname)
        SubWriter.write(sub_output_fname)

    
def flag(Data, NoiseData, thres=3.0, max_noise_factor=-1, modes_subtract=1,
         filter_type='edge'):
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
    modes_subtract : int
        How many modes to remove for high pass filtering.
    filter_type : {'edge', 'gaussian', 'gaussian/edge'}
        Type of high pass filtering to use.
    """
    
    # Get the mask and the data as normal arrays.
    # Copy seems to be nessisary if the mask is None.
    data = NoiseData.data.filled(0).copy()
    mask = ma.getmaskarray(NoiseData.data)
    ## High pass filter the data to make outliers stand out.
    un_mask = sp.logical_not(mask)
    NoiseData.calc_time()
    time = NoiseData.time
    n_time = len(time)
    # How many basis polynomials we need and with what fraction of each mode
    # gets subtracted out..
    if filter_type == 'edge':
        n_polys = modes_subtract
        subtract_weights = sp.ones(n_polys)
    elif filter_type == 'gaussian' or filter_type == 'gaussian/edge':
        n_polys = 4 * modes_subtract
        subtract_weights = sp.exp(-(sp.arange(n_polys, dtype=float)
                                     / modes_subtract)**2 / 2.)
        if filter_type == 'gaussian/edge':
            subtract_weights[0:2] = 1.
    # Test if the mask is the same for all slices.  If it is, that greatly
    # reduces the work as we only have to generate one set of polynomials.
    all_masks_same = True
    for jj in range(n_time):
        if sp.all(un_mask[jj,...] == un_mask[jj,0,0,0]):
            continue
        else:
            all_masks_same = False
            break
    if all_masks_same:
        polys = misc.ortho_poly(time, n_polys, un_mask[:,0,0,0], 0)
        polys.shape = (n_polys, len(time), 1, 1, 1)
    else:
        polys = misc.ortho_poly(time[:,None,None,None], n_polys, un_mask, 0)
    # Subtract the slope mode (1th mode) out of the NoiseData.
    amps = sp.sum(data * un_mask * polys, 1)
    amps *= subtract_weights[:,None,None,None]
    data -= sp.sum(amps[:,None,:,:,:] * un_mask[None,:,:,:,:] * polys, 0)
    ## Do the main outlier flagging.
    # Iteratively flag on sliding scale to get closer and closer to desired
    # threshold.
    max_thres = sp.sqrt(n_time)/2.
    n_iter = 3
    thresholds = (max_thres ** (n_iter - 1 - sp.arange(n_iter))
                 * thres ** sp.arange(n_iter)) ** (1./(n_iter - 1))
    for threshold in thresholds:
        # Subtract the mean from every channel.
        this_data = masked_subtract_mean(data, mask, 0)
        # Calculate the variance.
        un_mask = sp.logical_not(mask)
        counts = sp.sum(un_mask, 0)
        counts[counts == 0] = 1
        std = sp.sqrt(sp.sum(this_data**2 * un_mask, 0) / counts)
        bad_inds = abs(this_data) > threshold * std
        # If any polarization or cal state is masked, they all should be.
        bad_inds = sp.any(sp.any(bad_inds, 1), 1)
        mask[bad_inds[:,None,None,:]] = True
    ## Now look for times with excusion frequency average
    # (achromatic out-liers).
    # Compute the frequency mean.
    un_mask = sp.logical_not(mask)
    counts = sp.sum(un_mask, -1)
    fmean_un_mask = counts >= 1
    counts[counts == 0] = 1
    fmean = sp.sum(data * un_mask, -1) / counts
    # Subtract the time mean.
    fmean = masked_subtract_mean(fmean, sp.logical_not(fmean_un_mask), 0)
    # Get the variance.
    counts = sp.sum(fmean_un_mask, 0)
    counts[counts == 0] = 1
    fmean_std = sp.sqrt(sp.sum(fmean**2 * fmean_un_mask, 0) / counts)
    # Flag any time that is an outlier (for any polarization or cal state).
    bad_times = sp.any(sp.any(abs(fmean) > thres * fmean_std, 1), 1)
    mask[bad_times,:,:,:] = True
    ## Flag for very noisy channels.
    if max_noise_factor > 0:
        # Do this a few times to make sure we get everything.
        for ii in range(3):
            this_data = masked_subtract_mean(data, mask, 0)
            # Compute varience accounting for the mask.
            un_mask = sp.logical_not(mask)
            counts = sp.sum(un_mask, 0)
            vars_un_mask = counts >= 1
            counts[counts == 0] = 1
            vars = sp.sum(this_data**2 * un_mask, 0) / counts
            # Find the mean of the variences.
            counts = sp.sum(vars_un_mask, -1)
            counts[counts == 0] = 1
            mean_vars = sp.sum(vars * vars_un_mask, -1) / counts
            # Find channels that stand out (for any polarization or cal state).
            bad_chans = sp.any(sp.any(vars > max_noise_factor *
                                      mean_vars[:,:,None], 0), 0)
            mask[:,:,:,bad_chans] = True
    ## Transfer the mask to the DataBlock objects.
    Data.data[mask] = ma.masked
    NoiseData.data[mask] = ma.masked

def masked_subtract_mean(data, mask, axis):
    """Subtracts the mean from data along axis given a mask."""

    # Slice for broadcasting to the origional shape.
    up_broad = [slice(None)] * data.ndim
    up_broad[axis] = None
    up_broad = tuple(up_broad)
    # Calculate the mean.
    un_mask = sp.logical_not(mask)
    counts = sp.sum(un_mask, axis)
    counts[counts == 0] = 1
    mean = sp.sum(data * un_mask, axis) / counts
    # Subtract the mean.
    data = data - mean[up_broad] * un_mask
    return data








def nothing(noth):
    # If requested, remove the time gradient from all channels.
    if remove_slope:
        un_mask = sp.logical_not(ma.getmaskarray(NoiseData.data))
        NoiseData.calc_time()
        time = NoiseData.time
        n_time = len(time)
        # Test if the mask is the same for all slices.  If it is, that greatly
        # reduces the work as we only have to generate one set of polynomials.
        all_masks_same = True
        for jj in range(n_time):
            if sp.all(un_mask[jj,...] == un_mask[jj,0,0,0]):
                continue
            else:
                all_masks_same = False
                break
        if all_masks_same:
            polys = misc.ortho_poly(time, 2, un_mask[:,0,0,0], 0)
            polys.shape = (2, len(time), 1, 1, 1)
        else:
            polys = misc.ortho_poly(time[:,None,None,None], 2, un_mask, 0)
        # Subtract the slope mode (1th mode) out of the NoiseData.
        slope_amps = sp.sum(polys[1,...] * un_mask * NoiseData.data.filled(0),
                            0)
        NoiseData.data -= polys[1,...] * slope_amps
    # Iteratively flag on sliding scale to get closer and closer to desired
    # threshold.
    n_time = Data.data.shape[0]
    max_thres = sp.sqrt(n_time)/2.
    n_iter = 3
    thresholds = (max_thres ** (n_iter - 1 - sp.arange(n_iter))
                 * thres ** sp.arange(n_iter)) ** (1./(n_iter - 1))
    for threshold in thresholds:
        # Get the deviation from the mean.
        residuals = ma.anom(NoiseData.data, 0).filled(0)
        # Get indices above the threshold.
        mask = abs(residuals) > threshold * ma.std(NoiseData.data, 0)
        # Mask the data.
        Data.data[mask] = ma.masked
        NoiseData.data[mask] = ma.masked
    
    # Now flag for very noisey channels.
    if max_noise_factor > 0:
        vars = ma.var(NoiseData.data, 0)
        mean_vars = ma.mean(vars, -1).filled(0)
        bad_chans = vars.filled(0) > max_noise_factor * mean_vars[:,:,None]
        Data.data[:,bad_chans] = ma.masked
        NoiseData.data[:,bad_chans] = ma.masked


# If this file is run from the command line, execute the main function.
if __name__=="__main__":
    import sys
    ReFlag(str(sys.argv[1])).execute()
