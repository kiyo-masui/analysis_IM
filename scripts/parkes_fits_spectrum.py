from astropy.io import fits
import os
import numpy as np
import scipy.stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core import fitsGBT
import numpy.ma as ma
import datetime
from matplotlib.colors import LogNorm

#data_dir_orig = '/scratch2/p/pen/andersoc/second_parkes_pipe/rebinned/SDFITS/'
data_dir_rot = '/scratch2/p/pen/andersoc/second_parkes_pipe/rebinned/SDFITS_bpdivide/SVD_beams_split/ra165rot_to_I/'
data_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/rebinned/SDFITS_bpdivide/SVD_beams_split/ra165_I_map_subtracted_reflagged/'
#data_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/rebinned/SDFITS_bpdivide/SVD_beams_split/ra165_I_map_subtracted/'
out = '/scratch2/p/pen/andersoc/second_parkes_pipe/bp_estimates/all_data_spectrum/map_subtracted_I/'
out_orig = '/scratch2/p/pen/andersoc/second_parkes_pipe/bp_estimates/all_data_spectrum/'
#out = '/scratch2/p/pen/andersoc/second_parkes_pipe/bp_estimates/all_data_spectrum/'
#out = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/cal_nobp_jansky/'

def load_data(fname):
    hdu = fits.open(fname)
    return hdu[1].data

def load_datablock(fname):
    Reader = fitsGBT.Reader(fname)
    Blocks = Reader.read((), 0,force_tuple=True)
    return Blocks

def make_data_list(dir, mean=True, time=False, multiple_beams = True, which_beam = 0):
    #If mean is false, will calculate variance of each scan
    files = fname_list(dir, beam_choose = which_beam)
    count = 0
    time_bins = 0
    if multiple_beams:
        averages = ma.zeros((len(files),13,2,1,64))
        for file in files:
            blocks = load_datablock(file)
            for num in range(13):
                if mean==True:
                     averages[count,num,:,:,:] = ma.mean(blocks[num].data, axis=0)
                elif time == False:
                     averages[count,num,:,:,:] = ma.var(blocks[num].data, axis=0)
                if time == True:
                    time_bins += blocks[num].data.shape[0]
                    averages = time_bins
            count += 1
    else:
        averages = ma.zeros((len(files),1,2,1,64))
        for file in files: 
            blocks = load_datablock(file)
            if mean==True:
                averages[count,:,:,:,:] = ma.mean(blocks[0].data, axis=0)
            elif time == False:
                averages[count,:,:,:,:] = ma.var(blocks[0].data, axis=0)
            if time == True:
                time_bins += blocks[0].data.shape[0]
                averages = time_bins
            count += 1
    return averages
        
                

def data_list(dir):
    hdu_list = []
    for path, subdirs, files in os.walk(dir):
        for name in files:
            hdu_list.append(load_data(os.path.join(path, name)))
    return hdu_list

def fname_list(dir, sort ='True', beam_choose = 0):
    list = []
    for path, subdirs, files in os.walk(dir):
        for name in files:
            if beam_choose == 0:
                list.append(os.path.join(path, name))
            else:
                if 'beam' + str(beam_choose) + '.' in name:
                    list.append(os.path.join(path, name))
    if sort:
        list.sort(key = lambda x: datetime.datetime.strptime(x[x.find('/2014-') + 1:x.find('/2014-') + 16], '%Y-%m-%d_%H%M')) 
    return list

def combine_hdu(hdu_list, type):
    array_list = []
    for data in hdu_list:
        array_list.append(data[type])
    return np.concatenate(array_list)

def bp_estimate(beam, pol, dir):
    files = fname_list(dir)
    bp = avg_data(load_data(files[0]), beam, pol)
    bp = bp[None,:]
    for file in files[1:]:
        avg = avg_data(load_data(file), beam, pol)
        bp = np.concatenate((bp, avg[None,:]))
    bp = scipy.stats.nanmean(bp)
    return bp

def bp_sd(beam, pol, dir):
    files = fname_list(dir)
    bp = avg_data(load_data(files[0]), beam, pol)
    bp = bp[None,:]
    for file in files[1:]:
        avg = avg_data(load_data(file), beam, pol)
        bp = np.concatenate((bp, avg[None,:]))
    bp -= scipy.stats.nanmean(bp, axis=0)
    bp = bp**2
    bp = scipy.stats.nanmean(bp, axis=0)
    bp = bp**0.5
    return bp

def avg_data(data, beam, pol):
    #Will produce time average of data at particular beam, pol.
    #Could be a decent bandpass estimate, if Tspill isn't too freq dependent.
    #pol should be pol code, either -5 or -6 for XX, YY.
    beam_bool = data['BEAM'] == beam
    pol_bool = data['CRVAL4'] == pol
    bool = np.logical_and(beam_bool, pol_bool)
    data = np.array(data['DATA'][bool])
    data = scipy.stats.nanmean(data)
    return data

def plot_waterfall(data, beam, pol, file_name, log = False):
    #Takes in [scans, beams, pol, cal, freq] array.  Could be means, var, std.
    #Makes color plot with scans and freq along the two axes.
    #Beam is from 0 to 12, pol can be 0 (xx) or 1 (yy).
    data = data[:, beam, pol, 0, :]
    plt.axes()
    max = np.max(data)
    min = np.min(data)
    print min
    print max
    if log == True:
        plt.imshow(np.transpose(data), extent=[0, data.shape[0], 1347.5,1283.5] , interpolation = 'nearest', norm=LogNorm(vmin=np.min(data), vmax=np.max(data)))
    else:
        plt.imshow(np.transpose(data), vmax = max, vmin = min, interpolation = 'nearest')
    plt.axes().set_aspect(35.)
    plt.xlabel(r'Block number', fontsize=14, color='black')
    plt.ylabel(r'Frequency (MHz)', fontsize=14, color='black')
    plt.colorbar()
    plt.savefig(file_name)
    plt.clf()

def plot_bp_estimates(data_dir, file_name, norm = 'False'):
    fig, ax = plt.subplots(nrows=13, ncols=2, figsize=(12,18))
    #pols = [-5,-6]
    pols = [0,1]
    beams = range(1,14)
    col = 0
    #bps = ma.mean(make_data_list(data_dir), axis=0)
    #bpsds = ma.std(make_data_list(data_dir), axis=0)
    bps = ma.median(make_data_list(data_dir), axis=0)
    #bps = make_data_list(data_dir)
    #bps -= ma.median(bps, axis=0)
    #bps = (bps**2)**0.5
    #bps = ma.median(bps, axis = 0)
    #bps = ma.mean(make_data_list(data_dir, mean=False), axis=0)**0.5
    bps = np.transpose(bps)
    bps /= ma.mean(bps, axis=0)
    bps = np.transpose(bps)
    for pol in pols:
        for beam in beams:
            #bp = bp_estimate(beam, pol, data_dir)
            bp = bps[beam -1, pol, 0,:]
            #bp = bpsds[beam -1, pol, 0,:]
            #bp = bp_sd(beam, pol, data_dir)
            ax[beam - 1][col].plot(bp, 'b-')
        col += 1
    plt.savefig(file_name)
    plt.clf()

if __name__ == '__main__':
    #plot_bp_estimates(data_dir, out + 'masked_median_meanscan_normalized')
    #bps = ma.median(make_data_list(data_dir), axis=0)
    #Normalize bandpass estimates
    #bps = np.transpose(bps)
    #bps /= ma.mean(bps, axis=0)
    #bps = np.transpose(bps)
    #Save bandpass estimates
    #bps = ma.filled(bps)
    #np.save(out + 'masked_median_meanscan_bp', bps)
    #std_rot = make_data_list(data_dir_rot, mean=False, multiple_beams = False, which_beam=1)**0.5
    #mean_rot = make_data_list(data_dir_rot, mean=True, multiple_beams = False, which_beam=1)
    #np.ma.dump(std_rot, out + 'standard_dev_byscan_nomapsubtraction')
    #np.ma.dump(mean_rot, out + 'mean_byscan_nomapsubtraction')
    std_rot = np.load(out + 'standard_dev_byscan_nomapsubtraction')
    mean_rot = np.load(out + 'mean_byscan_nomapsubtraction')
    #print std.shape
    mean = np.load(out + 'mean_byscan')
    std = np.load(out + 'standard_dev_byscan')
    #mean_orig = np.load(out_orig + 'mean_byscan')
    #std_orig = np.load(out_orig + 'standard_dev_byscan')
    list = fname_list(data_dir)
    np.save(out + 'fits_list', np.array(list, dtype=object))
    #print len(orig_list)
    #file_cuts = [data_dir_orig + '20140428/2014-04-28_1343-P641_2df2_west1_1315_P641.fits', data_dir_orig + '20140505/2014-05-05_1242-P641_2df2_west2_1315_P641.fits']
    #cut_indeces = [orig_list.index(file_cuts[0]),orig_list.index(file_cuts[1])]
    #print cut_indeces
    #Get beam 1 and cut the 2 files.
    #std_orig = np.delete(std_orig, cut_indeces, axis=0)
    #mean_orig = np.delete(mean_orig, cut_indeces, axis=0) 
    #std = ma.log10(std)
    #plot_waterfall(std, 0, 0, out + 'beam1xx_stdevperscan_log')
    plot_waterfall(mean/mean_rot, 0, 0, out + 'beam1xx_mean_ratio_to_original', log=False)
    plot_waterfall(std/std_rot, 0, 0, out + 'beam1xx_std_ratio_to_original', log=False)
    print np.sum(std**-2)/np.sum(std_rot**-2)
    #plot_waterfall(mean/33., 0, 0, out + 'beam1xx_meanperscan_divby33Jy_log', log=True)
    #plot_waterfall(std/33., 0, 0, out + 'beam1xx_stdperscan_divby33Jy_log', log=True)
    #plot_waterfall(std, 0, 0, out + 'beam1xx_stdperscan_Jy_log', log=True)
    #plot_waterfall(std/mean, 0, 0, out + 'beam1xx_stddivbymean_perscan_Jy_log', log=True)
    #plot_waterfall(mean/33.,9,0, out + 'beam10xx_meanperscan_divby33Jy_log', log=True)
    #plot_waterfall(std/33.,9,0, out + 'beam10xx_stdperscan_divby33Jy_log', log=True)
    #plot_waterfall(std, 9, 0, out + 'beam10xx_stdperscan_Jy_log', log=True)
    #plot_waterfall(std/mean, 9, 0, out + 'beam10xx_stddivbymean_perscan_Jy_log', log=True)
    #time_bins = make_data_list(data_dir, mean=False, time = True, multiple_beams = False, which_beam = 1)
    #np.ma.dump(time_bins, out + 'time_bins')
    #print time_bins
    print np.sum(std[np.logical_not(std.mask)]**-2)/np.sum(std_rot[np.logical_not(std.mask)]**-2)
    print np.max((std[np.logical_not(std.mask)]**-2)/(std_rot[np.logical_not(std.mask)]**-2))
    print np.min((std[np.logical_not(std.mask)]**-2)/(std_rot[np.logical_not(std.mask)]**-2))
