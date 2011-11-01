"""Simple script to read in a GBT fits file and trim it down to make a test
case.

A note on versioning of the test file: I have set up the git repository to
intensionally ignore the output fits file of this script.  That is why it
doesn't show up on git status.  However the testfile is only 2MB, so if could
be versioned.  Alternately, we could version MD5sums of all data files used for
testing, or regression testing.
"""

import os
import copy
import shelve

import scipy as sp
import pyfits

from core import fitsGBT
from time_stream import rebin_freq, rebin_time, split_bands, combine_cal
from time_stream import rotate_pol
from noise import measure_noise


#### The origional spectrometer gbt test data file.
data_file_name  = (os.getenv('GBT10B_SPEC')  + 
                  '/01_wigglez22hr_azel_17-24.raw.acs.fits')
#                  '04_wigglez1hr_azel_113-120.raw.acs.fits')
test_file_name = './testdata/testfile_GBTfits.fits'
# Some parameters of GBT spectrometer data.  These are assumed to be
# correct with only minimal checking.
npol = 4 # number of polarizations
ncal = 2 # number of noise cal states (on, off)

# Read data and copy it
datahdulist = pyfits.open(data_file_name, 'readonly')
testhdulist = pyfits.HDUList(datahdulist)
fitsdata = datahdulist[1].data
# Get the scans and IF of all records and the unique elements
scans_all = fitsdata.field('SCAN') # Again by reference
scan_set = sp.unique(scans_all)
IFs_all = fitsdata.field('CRVAL1')/1E6 # MHz
IFs_all = IFs_all.round(0)
IF_set = sp.unique(IFs_all)

# Now select find the indicies we want 
scans = [0,1]
IFs = [0,1]
times = range(10*npol*ncal)
inds_total = []
for thescan in scan_set[scans] :
    for theIF in IF_set[IFs] :
        (inds_sif,) = sp.where(sp.logical_and(IFs_all==theIF, 
                                              scans_all==thescan))
        inds_sif = list(inds_sif[times])
        inds_total += inds_sif

inds_total.sort()
# Seems to be nessisary for fitsdata[inds] to be the right type
inds = sp.array(inds_total)

testhdulist[1].data = fitsdata[inds]
testhdulist.writeto(test_file_name)


#### A series of test data files created from guppi data.
guppi_file_name  = (os.getenv('GBT_DATA')  + 
                    '/GBT10B_036/42_wigglez15hrst_ralongmap_230-237.fits')
Reader = fitsGBT.Reader(guppi_file_name)
Blocks = Reader.read((0,1), None)
for Data in Blocks:
    rebin_freq.rebin(Data, 32, True, True)
    rebin_time.rebin(Data, 2)

split_Blocks = ()
for Data in Blocks:
    split_Blocks += split_bands.split(Data, 2, 32, 25)

comb_Blocks = copy.deepcopy(split_Blocks)
for Data in comb_Blocks:
    combine_cal.combine(Data, sub_mean=False)

rot_Blocks = copy.deepcopy(comb_Blocks)
for Data in rot_Blocks:
    rotate_pol.rotate(Data)

# Measure some parameters from the noise.
#out_db = shelve.open('./testdata/testfile_guppi_noise_parameters.shelve')
#parameters = get_correlated_overf(rot_Blocks, 1.0, window='hanning')

Writer = fitsGBT.Writer(Blocks)
Writer.write('./testdata/testfile_guppi_rebinned.fits')
Writer = fitsGBT.Writer(split_Blocks)
Writer.write('./testdata/testfile_guppi_split.fits')
Writer = fitsGBT.Writer(comb_Blocks)
Writer.write('./testdata/testfile_guppi_combined.fits')
Writer = fitsGBT.Writer(rot_Blocks)
Writer.write('./testdata/testfile_guppi_rotated.fits')

