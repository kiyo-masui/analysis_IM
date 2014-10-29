import os
import fnmatch
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy import optimize

from core import fitsGBT

plt.ion()

#fnames = ['2012-10-24_0711-P641_east1_1392_P641.sdfits']
#fnames = ['01_wigglez1hr_azel_121-128.fits']
#fnames = ['15_wigglez1hr_centre_ralongmap_122-131.fits', '15_wigglez1hr_centre_ralongmap_12-21.fits', '15_wigglez1hr_centre_ralongmap_132-141.fits', '15_wigglez1hr_centre_ralongmap_144-153.fits', '15_wigglez1hr_centre_ralongmap_154-163.fits', '15_wigglez1hr_centre_ralongmap_164-173.fits', '15_wigglez1hr_centre_ralongmap_176-185.fits', '15_wigglez1hr_centre_ralongmap_22-31.fits', '15_wigglez1hr_centre_ralongmap_32-41.fits', '15_wigglez1hr_centre_ralongmap_48-57.fits', '15_wigglez1hr_centre_ralongmap_58-67.fits', '15_wigglez1hr_centre_ralongmap_68-77.fits']
#fnames = ['02_3C295_track_104.fits'] #, '02_3C295_track_48.fits']
#root = os.getenv('GBT_DATA') + 'GBT12A_418/'

def find_pattern(pattern,root_dir):
    matches = []
    for root, dirnames, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))

    return matches

#root = '/mnt/raid-project/gmrt/kiyo/gbt_out/flagged/GBT10B_036/'
root = '/mnt/raid-project/gmrt/kiyo/gbt_out/flagged/GBT10B_036/'

f_names = find_pattern("**.fits",root)

#print f_names

fnames = []
for file in f_names:
    if int(file[55:57]) >= 42:
        #print file[55:57]
        fnames.append(file)


#Blocks = []

#for fname in fnames:
     #Read.
    #fpath = root + fname
    #Reader = fitsGBT.Reader(fpath)
    #Reader = fitsGBT.Reader(fname)
    #Data_list = Reader.read((), 0, force_tuple=True)
    #Data_list.shape
    #Blocks += list(Data_list)

plt.figure()
#a = ma.count_masked(Blocks[2].data[:,0,0,:], 0)
#b = np.size(Blocks[2].data[:,0,0,:], 0)
#Blocks[2].calc_freq()
#print Blocks[2].field.keys()
#c = Blocks[2].freq
#print len(c)
#print a
#test=0
#for file in a:
#    test = test + file
#print test
#plt.plot(Blocks[2].freq, a)
#plt.show()

def make_plot(files):
    i = 1
    Reader = fitsGBT.Reader(files[0])
    Data1 = list(Reader.read(0, 0, force_tuple = True))
    Data1[0].calc_freq()
    freq = Data1[0].freq
    num_bins = len(Data1[0].freq)
    a = np.arange(num_bins)
    flagged=np.zeros(num_bins)
    number = flagged
    for file in files:
        if i <= 30:
            Reader = fitsGBT.Reader(file)
            Data_list = list(Reader.read((), 0, force_tuple = True))
            for data in Data_list:
                data.calc_freq()
                current_spect = data.freq
                if current_spect.all() == freq.all():
                    flags = data.data[:,0,0,:]
                    flag_count = ma.count_masked(flags, 0)
                    #print 'flags'
                    #print np.sum(flag_count)
                    fbins = np.size(flags, 0)
                    #print 'total'
                    #print fbins
                    flagged = [x + y for x, y in zip(np.histogram(a, bins=num_bins, weights = flag_count)[0], flagged)]
                    number = [x + y for x, y in zip(fbins*np.histogram(a,bins=num_bins, weights = None)[0], number)]
                    i = i + 1
    flagged_n = [x/y for x, y in zip(flagged, number)]
    return {'domain':flagged_n, 'range':freq}

plot_lib = make_plot(fnames)
plt.bar(plot_lib['range'], plot_lib['domain'])
plt.show()
