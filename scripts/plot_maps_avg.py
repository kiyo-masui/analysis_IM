import core.algebra as al
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import h5py
from itertools import combinations
import glob
import os
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker

def plot_avg(file_list, av_ax, filename):
    count = 0
    rows = float(len(file_list))
    two = float(2)
    rows = int(math.ceil(rows/two))
    fig, ax = plt.subplots(nrows=rows, ncols=2, figsize=(18,18))
    for file in file_list:
        data = al.load(file)
        data = np.transpose(np.mean(data, axis=av_ax))
        print data.shape
        print np.max(data)
        ra_step = data.info['ra_delta']
        dec_step = data.info['dec_delta']
        ra_range = [data.info['ra_centre'] + ra_step*data.shape[1]/2, data.info['ra_centre'] - ra_step*data.shape[1]/2]
        dec_range = [data.info['dec_centre'] + dec_step*data.shape[0]/2, data.info['dec_centre'] - dec_step*data.shape[0]/2]
        data = np.fliplr(data)
        data = np.flipud(data)
        im = ax[count%rows][count/rows].imshow(data, extent = [ra_range[0], ra_range[1], dec_range[0], dec_range[1]])
        #ax[count/rows].imshow(data, extent = [ra_range[0], ra_range[1], dec_range[0], dec_range[1]])
        cbar_ax = fig.add_axes([0.85, 0.3, 0.05, 0.4])
        fig.colorbar(im, cbar_ax)
        count += 1
        print ra_range
        print dec_range
    plt.savefig(filename)
    plt.clf()

def plot_2df(file_list, noise_list, freq, fname, subtract_mean=True, avg_freqs = False, jansky=False, nweight = False):
    rows = int(len(file_list))
    #fig, ax = plt.subplots(figsize = (35,20), nrows=rows, ncols=1)
    #fig, ax = plt.subplots(figsize = (50,35), nrows=rows, ncols=1)
    fig, ax = plt.subplots(nrows=rows, ncols=1)
    #fig, ax = plt.subplot2grid((2,1),(0,0))
    #plt.subplots_adjust(wspace=0.01,hspace=0.01)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=None, hspace=None)
    #fig.subplots_adjust(hspace = 0.02, wspace = 0.02)
    #gs1 = gridspec.GridSpec(rows, 1)
    #gs1.update(wspace=0.025, hspace=0.05)
    count = 0
    maps = []
    high = 0
    low = 0
    for files in zip(file_list, noise_list):
        print files[0]
        print files[1]
        data = al.load(files[0])
        #Use line below only to correct calibration.
        #data = 0.5*al.load(files[0])
        noise = al.load(files[1])
        if nweight:
            noise = np.transpose(noise)
            print noise.shape
            norm = np.mean(np.mean(noise, axis = 0), axis = 0)
            #norm = np.max(np.max(noise, axis = 0), axis = 0)
            print norm.shape
            noise /= norm
            noise = np.transpose(noise)
            data = np.multiply(data,noise)
            noise[:] = 1.0
        if not avg_freqs:
            data = np.transpose(data[freq,:,:])
            noise = np.transpose(noise[freq,:,:])
        else:
            data = np.transpose(np.mean(data, axis=0))
            noise = np.transpose(np.mean(noise, axis=0))
        #Subtract noise weighted mean from slice
        if subtract_mean:
            data -= np.mean(np.multiply(data,noise))/(np.mean(noise)*1.0)
        maps.append(data)
        high = max(np.max(data),high)
        low = min(np.min(data), low)
        print data.shape
        s_d = (np.mean(data**2*noise)/np.mean(noise))**0.5
        print 'Map standard deviation in this slice is ' + str(s_d)
        print np.max(data)
    im = 0
    for data in maps:
        ra_step = data.info['ra_delta']
        dec_step = data.info['dec_delta']
        ra_range = [data.info['ra_centre'] + ra_step*data.shape[1]/2, data.info['ra_centre'] - ra_step*data.shape[1]/2]
        dec_range = [data.info['dec_centre'] + dec_step*data.shape[0]/2, data.info['dec_centre'] - dec_step*data.shape[0]/2]
        data = np.fliplr(data)
        data = np.flipud(data)
        #im = ax[count].imshow(data, cmap = cm.Greys_r, extent = [ra_range[0], ra_range[1], dec_range[1], dec_range[0]], vmax = high, vmin = low)
        im = ax[count].imshow(data, extent = [ra_range[0], ra_range[1], dec_range[1], dec_range[0]], vmax = high, vmin = low)
        if count != 0:
        #if count == 2:
            ax[count].set_xlabel('Right Ascension')
        #ax[count].set_ylabel('Declination')
        ax[count].yaxis.set_major_locator(ticker.MultipleLocator(2))
        count += 1
        print ra_range
        print dec_range
    #fig.subplots_adjust(right=0.7)
    #plt.tight_layout()
    cbar_ax = fig.add_axes([0.85, 0.3, 0.05, 0.4])
    #cbar_ax = fig.add_axes([0.8, 0.3, 0.05, 0.4])
    cbar = fig.colorbar(im, cbar_ax)
    if not jansky:
        cbar.set_label('Temperature fluctuations in K', rotation = 270, labelpad=13)
    else:
        cbar.set_label('Flux in Jk', rotation = 270, labelpad=14)
    #fig.tight_layout()
    #plt.savefig(fname,bbox_inches='tight')
    #plt.savefig(fname, pad_inches = 0.0)
    #plt.gca().set_axis_off()
    plt.subplots_adjust(top = 0.99, bottom = 0, right = 0.835, left = 0, 
            hspace = -0.68, wspace = 0.01)
    #plt.subplots_adjust(top = 0.99, bottom = 0, right = 0.835, left = 0,
    #        hspace = 0.2, wspace = 0.01)
    #plt.margins(0,0)
    #plt.gca().xaxis.set_major_locator(plt.NullLocator())
    #plt.gca().yaxis.set_major_locator(plt.NullLocator())
    #fig.text(0.5, 0.01, 'Right Ascension', ha='center')
    fig.text(-0.08, 0.5, 'Declination', va='center', rotation='vertical')
    #fig.text(0.00, 0.5, 'Declination', va='center', rotation='vertical')
    plt.savefig(fname, bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

def plot_wide(file_list, clean_list, noise_list, freq, fname, subtract_mean=True, avg_freqs = False, jansky=False, nweight = False):
    cols = int(len(file_list))
    #fig, ax = plt.subplots(nrows=3, ncols=cols, gridspec_kw = {'width_ratios':[486, 568]}, figsize = [12.0,9.0])
    fig, ax = plt.subplots(nrows=3, ncols=cols, gridspec_kw = {'width_ratios':[486, 568]}, figsize = [18.0,13.5])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                wspace=None, hspace=None)
    count = 0
    maps = []
    noises = []
    cleans = []
    high = 0
    low = 0
    for files in zip(file_list, clean_list, noise_list):
        print files[0]
        print files[1]
        data = al.load(files[0])
        noise = al.load(files[2])
        clean = al.load(files[1])
        if nweight:
            noise = np.transpose(noise)
            print noise.shape
            norm = np.mean(np.mean(noise, axis = 0), axis = 0)
            print norm.shape
            noise /= norm
            noise = np.transpose(noise)
            data = np.multiply(data,noise)
            clean = np.multiply(data,noise)
            noise[:] = 1.0
        if not avg_freqs:
            data = np.transpose(data[freq,:,:])
            clean = np.transpose(clean[freq,:,:])
            noise = np.transpose(noise[freq,:,:])
        else:
            data = np.transpose(np.mean(data, axis=0))
            clean = np.transpose(np.mean(clean, axis=0))
            noise = np.transpose(np.mean(noise, axis=0))
        #Subtract noise weighted mean from slice
        if subtract_mean:
            print np.mean(data)
            print np.mean(noise)
            #data -= np.mean(np.multiply(data,noise))/(np.mean(noise)*1.0)
            data -= np.mean(data)
            #clean -= np.mean(np.multiply(clean,noise))/(np.mean(noise)*1.0)
            clean -= np.mean(clean)
        #Convert clean map to mK
        sn = np.mean(np.multiply(np.abs(clean), noise**0.5))
        print 'Typical signal to noise is ' + str(sn)
        clean *= 1000.
        maps.append(data)
        noises.append(noise)
        cleans.append(clean)
        high = [0,0,0]
        low = [0,0,0]
        high[0] = max(np.max(data),high[0])
        low[0] = min(np.min(data), low[0])
        high[1] = max(np.max(clean),high[1])/4.
        low[1] = min(np.min(clean), low[1])/4.
        high[2] = max(np.max(noise),high[2])
        low[2] = min(np.min(noise), low[2])
        print data.shape
        s_d = (np.mean(data**2*noise)/np.mean(noise))**0.5
        print 'Map standard deviation in this slice is ' + str(s_d)
        print np.max(data)
    im = 0
    col = 0
    for datas in zip(maps, cleans, noises):
        count = 0
        for data in datas:
            ra_step = data.info['ra_delta']
            dec_step = data.info['dec_delta']
            ra_range = [data.info['ra_centre'] + ra_step*data.shape[1]/2, data.info['ra_centre'] - ra_step*data.shape[1]/2]
            dec_range = [data.info['dec_centre'] + dec_step*data.shape[0]/2, data.info['dec_centre'] - dec_step*data.shape[0]/2]
            data = np.fliplr(data)
            data = np.flipud(data)
            #im = ax[count, col].imshow(data, extent = [ra_range[0], ra_range[1], dec_range[1], dec_range[0]], vmax = high[count], vmin = low[count], cmap = plt.get_cmap('plasma'))
            im = ax[count, col].imshow(data, extent = [ra_range[0], ra_range[1], dec_range[1], dec_range[0]], vmax = high[count], vmin = low[count])
            if col != 0:
                 #ax[count, col].set_xlabel('Right Ascension')
                 ax[count, col].set_yticklabels([])
            if count != 2:
                 ax[count, col].set_xticklabels([])
            ax[count, col].yaxis.set_major_locator(ticker.MultipleLocator(2))
            count += 1
            print ra_range
            print dec_range
            if count == 1:
                cbar_ax = fig.add_axes([0,0.2, 0.25,0.05])
                cbar = fig.colorbar(im, cbar_ax, orientation = 'horizontal')
                cbar.set_label('Temperature fluctuations in K', rotation = 0, labelpad=13)
            if count == 2:
                cbar_ax = fig.add_axes([0.29 ,0.2, 0.25,0.05])
                cbar = fig.colorbar(im, cbar_ax, orientation = 'horizontal')
                cbar.set_label('Temperature fluctuations in mK', rotation = 0, labelpad=13)
            if count == 3:
                cbar_ax = fig.add_axes([0.58 ,0.2, 0.25,0.05])
                cbar = fig.colorbar(im, cbar_ax, orientation = 'horizontal')
                cbar.set_label(r'Inverse noise weight in $\mathrm{K^{-2}}$', rotation = 0, labelpad=13)
                tick_locator = ticker.MaxNLocator(nbins=5)
                cbar.locator = tick_locator
                cbar.update_ticks()
        col += 1
    plt.subplots_adjust(top = 0.99, bottom = 0, right = 0.835, left = 0, 
            hspace = -0.87, wspace = 0.01)
    fig.text(-0.06, 0.5, 'Declination', va='center', rotation='vertical')
    fig.text(0.38, 0.29, 'Right Ascension', va='center', rotation='horizontal')
    #plt.savefig(fname, bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(fname, bbox_inches = 'tight', pad_inches = 0.1, format='eps', dpi=300)
    plt.clf()



def save_2df_map_avg(file_list, noise_list, path_out, fname):
    print len(file_list)
    print len(noise_list)
    count = 0
    map_sum = 0
    noise_sum = 0
    for files in zip(file_list, noise_list):
        print files[0]
        print files[1]
        data = al.load(files[0])
        noise = al.load(files[1])
        data = np.multiply(data,noise)
        if count == 0:
            map_sum = data
            noise_sum = noise
        else:
            map_sum += data
            noise_sum += noise
        count += 1
    map_sum = np.divide(map_sum, noise_sum)
    al.save(path_out + fname + 'map' + '.npy', map_sum)
    al.save(path_out + fname + 'weights' + '.npy', noise_sum)



if __name__ == '__main__':
    base_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/conv_fixed/renormalized/'
    noise_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_to_share/'
    output_dir = '/home/p/pen/andersoc/plots/'
    #Plot beam 1 xx maps for 2df1, at freq 32.  
    combine_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/conv_fixed/renormalized/beams_combined/'
    #cleaned_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_bp_divide/beams_removed/flux_cal/proper_conv/XX/'
    #cleaned_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_bp_divide/hitconv_sync27_mbcal/I/'
    cleaned_dir = '/scratch2/p/pen/nluciw/parkes/maps/correct_effic/I/'

#for el in os.walk(base_dir):
#    print el[0]
#noise_files = glob.glob(base_dir + '*' + 'ra165dec0' + '*noise_inv*' + '.npy')
#map_files = glob.glob(base_dir + '*' + 'ra165dec0' + '*clean*' + '.npy')
#noise_files = glob.glob(base_dir + '*' + 'ra33decn30' + '*noise_inv*' + '.npy')
    #map_files = glob.glob(base_dir + '*' + 'ran18decn30' + '*beam*' + '*clean*'  + '*XX*' + '.npy')
    #map_files = glob.glob(combine_dir + 'ran18decn30*' + '*map*' + '.npy')
    #map_files = glob.glob(base_dir + '*' + 'ran18decn30' + '*beam*' + '*clean*'  + '*XX*' + '.npy')
    #noise_files = glob.glob(noise_dir + '*' + 'ran18decn30' + '*beam*' + '*noise_inv_diag*'  + '*XX*' + '.npy')
    #noise_files = glob.glob(combine_dir + 'ran18decn30*' + '*weights*' + '.npy')

    #map_files.append(glob.glob(combine_dir + 'ra33decn30*' + '*map*' + '.npy')[0])
    #noise_files.append(glob.glob(combine_dir + 'ra33decn30*' + '*weights*' + '.npy')[0])

    #map_files = [cleaned_dir + 'ran18/' + 'combined_clean_map_10modes.npy', cleaned_dir + 'ra33/' + 'combined_clean_map_10modes.npy']
    #noise_files = [cleaned_dir + 'ran18/' + 'combined_clean_weight_10modes.npy', cleaned_dir + 'ra33/' + 'combined_clean_weight_10modes.npy']
    clean_files = [cleaned_dir + 'ran18/' + 'combined_clean_map_10modes.npy', cleaned_dir + 'ra33/' + 'combined_clean_map_10modes.npy']
    dir_noise_test = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/noise_reestimated/correct_effic_sync27/'
    clean_files = [dir_noise_test + 'parkes_parallel_thread_ra165dec0_p08_440by136_beam1,2,3_dirty_map_I_1316.npy']

    #map_files = [cleaned_dir + 'ra165/' + 'combined_clean_map_10modes.npy', cleaned_dir + 'ra199/' + 'combined_clean_map_10modes.npy', cleaned_dir + 'ra216/' + 'combined_clean_map_10modes.npy']
    #noise_files = [cleaned_dir + 'ra165/' + 'combined_clean_weight_10modes.npy', cleaned_dir + 'ra199/' + 'combined_clean_weight_10modes.npy', cleaned_dir + 'ra216/' + 'combined_clean_weight_10modes.npy']

    #map_files = [cleaned_dir + 'ra165/' + 'combined_clean_map_0modes.npy', cleaned_dir + 'ra199/' + 'combined_clean_map_0modes.npy', cleaned_dir + 'ra216/' + 'combined_clean_map_0modes.npy']
    #noise_files = [cleaned_dir + 'ra165/' + 'combined_clean_weight_0modes.npy', cleaned_dir + 'ra199/' + 'combined_clean_weight_0modes.npy', cleaned_dir + 'ra216/' + 'combined_clean_weight_0modes.npy']

    #map_files = [cleaned_dir + 'ran18/' + 'combined_clean_map_20modes.npy', cleaned_dir + 'ra33/' + 'combined_clean_map_20modes.npy']
    #noise_files = [cleaned_dir + 'ran18/' + 'combined_clean_weight_20modes.npy', cleaned_dir + 'ra33/' + 'combined_clean_weight_20modes.npy']
    orig_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps/'

    map_files = [cleaned_dir + 'ran18/' + 'combined_clean_map_0modes.npy', cleaned_dir + 'ra33/' + 'combined_clean_map_0modes.npy']
    noise_files = [cleaned_dir + 'ran18/' + 'combined_clean_weight_0modes.npy', cleaned_dir + 'ra33/' + 'combined_clean_weight_0modes.npy']
    map_files = [dir_noise_test + 'parkes_parallel_thread_ra165dec0_p08_440by136_beam1,2,3_dirty_map_I_1316.npy']
    noise_files = [dir_noise_test + 'parkes_parallel_thread_ra165dec0_p08_440by136_beam1,2,3_noise_inv_diag_I_1316.npy']
    orig_map_files = [orig_dir + 'parkes_parallel_thread_ra165dec0_p08_440by136_beam1_dirty_map_XX_1316.npy']
    orig_noise_files = [orig_dir + 'parkes_parallel_thread_ra165dec0_p08_440by136_beam1_noise_inv_diag_XX_1316.npy']
    orig_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/correct_beam_effic/sync27/combined/I/'
    orig_map_files = [orig_dir + 'parkes_parallel_thread_ra165dec0_p08_440by136_beam123_clean_map_bp_div_I_1316.npy']
    orig_noise_files = [orig_dir + 'parkes_parallel_thread_ra165dec0_p08_440by136_beam123_noise_inv_diag_bp_div_I_1316.npy']



    jk_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/cal_nobp_jansky/combined/I/'
    #map_files = [jk_dir + 'parkes_parallel_thread_ran18decn30_p08_488by106_beam_all__clean_map_bp_div_I_1316.npy', jk_dir + 'parkes_parallel_thread_ra33decn30_p08_568by106_beam_all__clean_map_bp_div_I_1316.npy']
    #map_files = [jk_dir + 'parkes_parallel_thread_ra165dec0_p08_440by136_beam_all__clean_map_bp_div_I_1316.npy', jk_dir + 'parkes_parallel_thread_ra182dec0_p08_440by136_beam_all__clean_map_bp_div_I_1316.npy', jk_dir + 'parkes_parallel_thread_ra199dec0_p08_440by136_beam_all__clean_map_bp_div_I_1316.npy', jk_dir + 'parkes_parallel_thread_ra216dec0_p08_440by136_beam_all__clean_map_bp_div_I_1316.npy']
    #noise_files = [jk_dir + 'parkes_parallel_thread_ran18decn30_p08_488by106_beam_all__noise_inv_diag_bp_div_I_1316.npy', jk_dir + 'parkes_parallel_thread_ra33decn30_p08_568by106_beam_all__noise_inv_diag_bp_div_I_1316.npy']
    #noise_files = [jk_dir + 'parkes_parallel_thread_ra165dec0_p08_440by136_beam_all__noise_inv_diag_bp_div_I_1316.npy', jk_dir + 'parkes_parallel_thread_ra182dec0_p08_440by136_beam_all__noise_inv_diag_bp_div_I_1316.npy', jk_dir + 'parkes_parallel_thread_ra199dec0_p08_440by136_beam_all__noise_inv_diag_bp_div_I_1316.npy', jk_dir + 'parkes_parallel_thread_ra216dec0_p08_440by136_beam_all__noise_inv_diag_bp_div_I_1316.npy']
    #noise_files = []
    #for file in map_files:
    #    noise_files.append(glob.glob(noise_dir + '*' + 'ran18decn30*' + file[file.index('beam'):file.index('clean')] + '*noise_inv_diag*'  + '*XX*' + '.npy')[0])

    #map_files = [map_files[0],map_files[2],map_files[3]]
    #noise_files = [noise_files[0],noise_files[2],noise_files[3]]
    #plot_2df(map_files, noise_files, 48, output_dir + '2df2_I_freqavg_allgoodbeams_calibratedJy_fixdec', subtract_mean=True, avg_freqs = True, jansky=True)
    #plot_2df(map_files, noise_files, 48, output_dir + '2df1_I_freqavg_allgoodbeams_calibratedJy_fixdec', subtract_mean=True, avg_freqs = True, jansky=True)
    #plot_2df(map_files, noise_files, 32, output_dir + '2df2_I_f32_allgoodbeams_calibrated10modes_fixdec', subtract_mean=True)
    #plot_2df(map_files, noise_files, 32, output_dir + '2df1_I_f32_allgoodbeams_calibrated10modes_fixdec', subtract_mean=True)

    #plot_wide(map_files, clean_files, noise_files, 32, output_dir + '2df1_I_f32_allgoodbeams_calibrated_wide_fixdec', subtract_mean=True)
    #plot_wide(map_files, clean_files, noise_files, 32, output_dir + '2df1_I_f32_allgoodbeams_calibrated_wide_fixdec_300dpi.eps', subtract_mean=True)
    #plot_wide(map_files, clean_files, noise_files, 32, output_dir + 'ra165_test_2ndround_mapmaker.eps', subtract_mean=True)
    plot_avg(map_files + noise_files + orig_map_files + orig_noise_files, 0, output_dir + 'ra165_test_2dnround_mapmaker')

    #plot_2df(map_files, noise_files, 32, output_dir + '2df1_I_f32_allgoodbeams_calibrated10modes_fixdec_nweighted', subtract_mean=True, nweight=True)
    #plot_2df(map_files, noise_files, 32, output_dir + '2df2_I_f32_allgoodbeams_calibrated10modes_fixdec_nweighted', subtract_mean=True, nweight=True)
    #save_2df_map_avg(map_files, noise_files, combine_dir, 'ran18decn30_XX_')
#noise_files.sort()
#map_files.sort()

#print map_files
#print len(map_files)

#plot_dir = '/home/p/pen/andersoc/plots/'

#plot_avg(noise_files, 0, plot_dir + 'noise_diags') 
#plot_avg(map_files, 0, plot_dir + 'freq_avg_maps_ra33decn30')
