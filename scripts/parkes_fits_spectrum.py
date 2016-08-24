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

data_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/rebinned/SDFITS/'
out = '/scratch2/p/pen/andersoc/second_parkes_pipe/bp_estimates/all_data_spectrum/'

def load_data(fname):
    hdu = fits.open(fname)
    return hdu[1].data

def load_datablock(fname):
    Reader = fitsGBT.Reader(fname)
    Blocks = Reader.read((), 0,force_tuple=True)
    return Blocks

def make_data_list(dir, mean=True):
    #If mean is false, will calculate variance of each scan
    files = fname_list(dir)
    count = 0
    averages = ma.zeros((len(files),13,2,1,64))
    for file in files:
        blocks = load_datablock(file)
        for num in range(13):
            if mean==True:
                averages[count,num,:,:,:] = ma.mean(blocks[num].data, axis=0)
            else:
                averages[count,num,:,:,:] = ma.var(blocks[num].data, axis=0)
        count += 1
    return averages
        
                

def data_list(dir):
    hdu_list = []
    for path, subdirs, files in os.walk(dir):
        for name in files:
            hdu_list.append(load_data(os.path.join(path, name)))
    return hdu_list

def fname_list(dir, sort ='True'):
    list = []
    for path, subdirs, files in os.walk(dir):
        for name in files:
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

def plot_waterfall(data, beam, pol, file_name):
    #Takes in [scans, beams, pol, cal, freq] array.  Could be means, var, std.
    #Makes color plot with scans and freq along the two axes.
    #Beam is from 0 to 12, pol can be 0 (xx) or 1 (yy).
    data = data[:, beam, pol, 0, :]
    plt.axes()
    max = np.max(data)
    min = np.min(data)
    print min
    print max
    plt.imshow(np.transpose(data), vmax = max, vmin = min, interpolation = 'nearest')
    plt.axes().set_aspect(35.)
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
    #bps = ma.median(make_data_list(data_dir), axis=0)
    #bps = make_data_list(data_dir)
    #bps -= ma.median(bps, axis=0)
    #bps = (bps**2)**0.5
    #bps = ma.median(bps, axis = 0)
    bps = ma.mean(make_data_list(data_dir, mean=False), axis=0)**0.5
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
    #plot_bp_estimates(data_dir, out + 'masked_sqrt_meanscanvar')
    std = make_data_list(data_dir, mean=False)**0.5
    std = ma.log10(std)
    plot_waterfall(std, 0, 0, out + 'beam1xx_stdevperscan_log')
