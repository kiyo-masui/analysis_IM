# Input parameter file for noise model code.

# Construct a list of filenames by looping over scans.
#input_root = "/mnt/raid-cita/tchang/wiggleZ/GBT10B_036/04_3c48_onoff_"
input_root = "/mnt/raid-cita/tchang/wiggleZ/GBT10B_036/00_3c48_onoff_"
# String
input_middles = ["161-164","165-168"] # List of strings
#input_middles = ["297-300","304-307"] # List of strings
input_ends = ".raw.acs.fits" # String

# How the files should be grouped:
files_scans_per_mean = 1
# Positive integer: take mean every this many files.
# negitive integer: take mean every -ev this many scans.
# 0: all files.

# Outputs
output_dir = "/mnt/raid-project/gmrt/kiyo/wiggleZ/noisemodel/04/1hr/" # String
output_root = "testing"

# Choose what IFs and scans to pick out of each file (applied uniformly accross
# files).
scans = [1,3] # empty for all
#scans = [0,2] # empty for all
IFs = [3] # empty for all

# Set up time lag bins
import scipy as sp
lags = sp.arange(0.01, 440.0, 30.0)

# Frequency bin settings
first_freq = 192
nfreqs = 1664
# Only the first nfreqs of all the frequencies read will be processed.  Must be
# compatible with IFs.
# Frequencies rebinned in blocks of freqs_per_bin
freqs_per_fbin = 64

