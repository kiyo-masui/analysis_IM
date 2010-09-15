# Input parameter file for noise model code.

# Construct a list of filenames by looping over scans.
input_root = "/mnt/raid-cita/tchang/wiggleZ/GBT10B_036/04_wigglez1hr_azel_" 
input_middles = [#"113-120","121-128","129-136","137-144"#,
                 "241-248","249-256","257-264","265-272",
                 "273-280","281-288"
                 #,"308-315","316-323","324-331","332-339"
                 ] # List of strings
input_ends = ".raw.acs.fits"

# Outputs
write_outputs_file = False # Doesn't work yet (except saving figures).
output_dir = "/mnt/raid-project/gmrt/kiyo/wiggleZ/noisemodel/04/1hr/" 
output_root = "caloff_nofreq_IF3_"

# How the files should be grouped.  All the data in a single block is processed
# at the same time, with the time median subtracted.  This parameter sets the
# size of the block:
# Results from all blocks are essentially averaged at the end.
# Positive integer: block is this many files (all scans per file).
# negitive integer: block is -ev this many scans.
# 0: block all files.
files_scans_per_block = -1

# Choose what IFs and scans to pick out of each file (applied uniformly accross
# files).
scans = range(8) # Empty for all but files_scans_per_block >= 0 for this.
IFs = [3] # empty for all
# Frequency bin settings
first_freq = 400
nfreqs = 1200
# Only the first nfreqs of all the frequencies read will be processed.  Must be
# compatible with IFs.
# Frequencies rebinned in blocks of freqs_per_bin
freqs_per_fbin = 80
# If the following is set to true, time fluctuations are scaled by the time RMS
# and then the freqency median is subtracted.

# Polarizations (XX, XY, YX, YY) and Cal States (On, Off)
pol_weights = [1., 0., 0., 1.]
cal_weights = [0., 1.]
# Divide time stream by the time mean of the below weighted time stream.
# All zeros will matches pol_weights and/or cal_weights.
norm_pol_weights = [1., 0., 0., 1.]
norm_cal_weights = [0., 1.]

calculate_covariance = True
calculate_power_spectrum = True
subtract_freq_average = False

# Set up time lag bins
import scipy as sp
lags = [0.9]
for ii in range(10) :
    lags.append(1.5*lags[-1])
                                                    
make_figures = True
# Figure options
plot_save_file = True
plot_fbins = range(1,nfreqs//freqs_per_fbin,3)
plot_norm_to_first_lag = False
plot_show = False
plot_output_root = output_root
plot_output_dir = output_dir
