"""Dirty map making module.

Module converts data in the time domain into noise weighted data in the map
domain, i.e. it creats the dirty map.  Module also contains many utilities like
the pointing operator (`Pointing`) and the time domain noise operator
(`Noise`).
"""
import os
os.environ['PYTHON_EGG_CACHE'] = '/scratch/p/pen/andersoc/.python-eggs'
import sys
from memory_profiler import profile

from mpi4py import MPI

import math
import threading
from Queue import Queue
import shelve
import sys
import time as time_mod
import cPickle
import h5py

import scipy as sp
import numpy.ma as ma
import scipy.fftpack as fft
from scipy import linalg
from scipy import interpolate
import numpy as np
#import kiyopy.pickle_method
import warnings
#import matplotlib.pyplot as plt

import core.algebra as al
from core import fitsGBT
import utils.misc as utils
import tools
from noise import noise_power
from foreground import ts_measure
import kiyopy.custom_exceptions as ce
import kiyopy.utils
from kiyopy import parse_ini
import _mapmaker as _mapmaker_c
from constants import T_infinity, T_huge, T_large, T_medium, T_small, T_sys
from constants import f_medium, f_large
from utils import misc
from map import dirty_map as dm

prefix ='dm_'
params_init = {# IO:
               'input_root' : './',
               # The unique part of every fname
               'file_middles' : ("testfile_GBTfits",),
               'input_end' : ".fits",
               'output_root' : "./testoutput_",
               # Shelve file that holds measurements of the noise.
               # What data to include from each file.
               'scans' : (),
               'IFs' : (0,),
               'polarizations' : ('I',),
               # Map parameters (Ra (deg), Dec (deg)).
               'field_centre' : (325.0, 0.0),
               # In pixels.
               'map_shape' : (5, 5),
               'pixel_spacing' : 0.5, # degrees
               # Interpolation between pixel points.  Options are 'nearest',
               # 'linear' and 'cubic'.
               'interpolation' : 'linear',
               # How to treat the data.
               # How much data to include at a time (scan by scan or file by
               # file)
               'time_block' : 'file',
               # Number of files to process at a time.  Increasing this number
               # reduces IO but takes more memory.
               'n_files_group' : 0,
               # What kind of frequency correlations to use.  Options are
               # 'None', 'mean' and 'measured'.
               'frequency_correlations' : 'None',
               # Ignored unless 'frequency_correlations' is 'measured'.
               'number_frequency_modes' : 1,
               # Of the above modes, how many to completely discard instead of
               # noise weighting.
               'number_frequency_modes_discard' : 0,
               # Where to get noise parameters from.
               'noise_parameter_file' : '',
               'deweight_time_mean' : True,
               'deweight_time_slope' : False,
               # If there where any foregrounds subtracted in the time stream,
               # let the noise know about it.
               'n_ts_foreground_modes' : 0,
               'ts_foreground_mode_file' : '',
               'thread_divide' : False,
               # Which beam to use to make map.  Choose zero to use all beams.
               'beam' : 0
               }

comm = MPI.COMM_WORLD

class DirtyMapMaker(object):
    """Dirty map maker.
    """

    def __init__(self, parameter_file_or_dict=None, feedback=2):
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix=prefix, feedback=feedback)
        print self.params['file_middles']
        self.feedback = feedback
        self.rank = comm.Get_rank()
        self.nproc = comm.Get_size()
    
    def iterate_data(self, file_middles):
        """An iterator over the input data.
        
        This can either iterate over the file names, processing the whole file
        at once, or it can iterate over both the file and the scans.

        Returns the DataBlock objects as well as the file middle for that data
        (which is used as a key in several parameter data bases).
        """
    
        params = self.params
        for middle in file_middles:
            self.file_middle = middle
            fname = params["input_root"] + middle + params["input_end"]
            Reader = fitsGBT.Reader(fname, feedback=self.feedback)
            Blocks = Reader.read(self.params['scans'], self.band_ind,
                                 force_tuple=True)
            if params['time_block'] == 'scan':
                for Data in Blocks:
                    if params["beam"] == 0 or params["beam"] == Data.field['BEAM']: 
                        yield (Data,), middle
            elif params['time_block'] == 'file':
                if params["beam"] == 0:
                    yield Blocks, middle
                else:
                    yield tuple([data for data in Blocks if data.field['BEAM']
== params["beam"]]), middle
            else:
                msg = "time_block parameter must be 'scan' or 'file'."
                raise dm.ValueError(msg)
        
        
    def execute(self, n_processes):
        """Driver method."""
        
        # n_processes matters when using mpirun from MPI.COMM_WORLD and threading
        self.n_processes = n_processes
        params = self.params
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root']
                               + 'dirty_map_params.ini',
                               prefix='dm_')
        if self.feedback > 0:
            print "Input root is: " + params["input_root"]
            print "Output root is: " + params["output_root"]
        # Set flag if we are doing independant channels.
        if params['frequency_correlations'] == 'None':
            self.uncorrelated_channels = True
            if self.feedback > 0:
                print "Treating frequency channels as being independant."
        else:
            self.uncorrelated_channels = False
        # Open the file that stores noise parameters.
        if params['noise_parameter_file']:
            self.noise_params = shelve.open(params['noise_parameter_file'], 'r')
        else:
            self.noise_params = None
        # If we are subtracting out foreground modes, open that file.
        if params['ts_foreground_mode_file']:
            self.foreground_modes = shelve.open(
                params['ts_foreground_mode_file'], 'r')
        else:
            self.foreground_modes = None
        # Set the map dimensioning parameters.
        self.n_ra = params['map_shape'][0]
        self.n_dec = params['map_shape'][1]
        spacing = params["pixel_spacing"]
        # Negative sign because RA increases from right to left.
        ra_spacing = -spacing/sp.cos(params['field_centre'][1]*sp.pi/180.)
        self.ra_spacing = ra_spacing
        # To set some parameters, we have to read the first data file.
        first_file_name = (params["input_root"] + params["file_middles"][0]
                           + params["input_end"])
        Reader = fitsGBT.Reader(first_file_name, feedback=self.feedback)
        Blocks = Reader.read(0, (), force_tuple=True)
        self.n_chan = Blocks[0].dims[3]
        # Figure out which polarization indices to process.
        n_pols = Blocks[0].dims[1]
        pols = list(Blocks[0].field["CRVAL4"])
        self.pols = pols
        if not params['polarizations']:
            pol_inds = range(n_pols)
        else:
            pol_inds = []
            for stokes_par in params['polarizations']:
                for ii in range(n_pols):
                    if stokes_par == utils.polint2str(pols[ii]):
                        pol_inds.append(ii)
                        break
                else :
                    msg = "Polarization " + stokes_par + " not in Data files."
                    raise ce.ValueError(msg)
        # Get the centre frequncies of all the bands and figure out which bands
        # to process.
        n_bands = len(Reader.IF_set)
        band_inds = params["IFs"]
        if not band_inds:
            band_inds = range(n_bands)
        delta_freq = Blocks[0].field['CDELT1']
        self.delta_freq = delta_freq
        band_centres = []
        for Data in Blocks:
            Data.calc_freq()
            band_centres.append(Data.freq[self.n_chan//2])
        self.band_centres = band_centres
        del Blocks
        del Reader
        map_shape = (self.n_chan, self.n_ra, self.n_dec)
        for ii in pol_inds:
            for jj in band_inds:
                self.band_ind = jj
                band_centre = band_centres[jj]
                self.pol_ind = ii
                pol = pols[ii]
                if self.feedback > 1:
                    print ("Making map for polarization "
                           + utils.polint2str(pol) + " and band centred at "
                           +repr(round(band_centre/1e6)) + "MHz.")
                # Work out the file names.
                map_filename = (params["output_root"] + "dirty_map_"
                                + utils.polint2str(pol)
                                + "_" + str(int(round(band_centre/1e6)))
                                + '.npy')
                cov_filename = (params["output_root"] + "noise_inv_"
                                + utils.polint2str(pol)
                                + "_" + str(int(round((band_centre/1e6))))
                                + '.hdf5')
                self.map_filename = map_filename
                self.cov_filename = cov_filename
                # Initialization of the outputs.
                map = sp.zeros(map_shape, dtype=float)
                map = al.make_vect(map, axis_names=('freq', 'ra', 'dec'))
                map.set_axis_info('freq', band_centre, delta_freq)
                map.set_axis_info('ra', params['field_centre'][0], ra_spacing)
                map.set_axis_info('dec', params['field_centre'][1], 
                                   params['pixel_spacing'])
                '''if self.uncorrelated_channels:
                    cov_inv = al.open_memmap(cov_filename, mode='w+',
                                             shape=map_shape + map_shape[1:],
                                             dtype=float)
                    cov_inv = al.make_mat(cov_inv,
                                axis_names=('freq', 'ra', 'dec', 'ra', 'dec'),
                                row_axes=(0, 1, 2), col_axes=(0, 3, 4))
                else:
                    cov_inv = al.open_memmap(cov_filename, mode='w+',
                                             shape=map_shape + map_shape,
                                             dtype=float)
                    cov_inv = al.make_mat(cov_inv,
                                      axis_names=('freq', 'ra', 'dec',
                                                  'freq', 'ra', 'dec'),
                                      row_axes=(0, 1, 2), col_axes=(3, 4, 5))
                cov_inv.copy_axis_info(map)'''
                # The zeroing takes too long. Do it in memory, not hard disk.
                #print 'zeroing cov_inv'
                #cov_inv[...] = 0
                #print 'zeroed'
                self.map = map
                #self.cov_inv = cov_inv
                # Do work.
                try:
                    print "self.make_map() now"
                    self.make_map()
                except:
                    # I think yo need to do this to get a sensible error
                    # message, otherwise it will print cov_inv and fill the
                    # screen.
                    #del self.cov_inv, cov_inv
                    raise
                # IO.
                # To write the noise_inverse, all we have to do is delete the
                # memeory map object.
                #del self.cov_inv, cov_inv
                if self.rank == 0:
                    al.save(map_filename, map)
        if not self.noise_params is None:
            self.noise_params.close()
    
    def make_map(self):
        """Makes map for current polarization and band.
        
        This worker function has been split off for testing reasons.
        """
        
        params = self.params
        map = self.map
        #cov_inv = self.cov_inv
        # The current polarization index and value.
        pol_ind = self.pol_ind
        pol_val = self.pols[pol_ind]
        # The current band, index and value.
        band_ind = self.band_ind
        band_centre = self.band_centres[band_ind]
        # This outer loop is over groups of files.  This is so we don't need to
        # hold the noise matrices for all the data in memory at once.
        file_middles = params['file_middles']
        n_files = len(file_middles)
        n_files_group = params['n_files_group']
        if n_files_group == 0:
            n_files_group = n_files
        for start_file_ind in range(0, n_files, n_files_group):
            all_this_file_middles = file_middles[start_file_ind
                                             :start_file_ind + n_files_group]
            # Each process will handle an equal number of input files.
            # (Numbers that do not divide are taken care of as evenly as
            # possible). This way is most efficient in memory and efficiency.
            # This also allows a smaller number of mpi calls when passing
            # DataSets around later to fill cov_inv.

            #this_file_middles,junk = split_elems(all_this_file_middles,
                                             #self.nproc)
            #this_file_middles = this_file_middles[self.rank]

            this_file_middles = all_this_file_middles

            # Don't need the 'junk' from above right now. That was added
            # at a later point for something else.
            if self.feedback > 1:
                # Will have to change for each process saying something
                # unique.
                print ("Processing data files %d to %d of %d."
                       % (start_file_ind, start_file_ind + n_files_group,
                          n_files))
            # Initialize lists for the data and the noise.
            time_stream_list = []
            data_set_list = []
            # Loop an iterator that reads and preprocesses the data.
            # This loop can't be threaded without a lot of restructuring,
            # because iterate_data, preprocess_data and get_noise_parameter
            # all talk to each other through class attributes.
            for this_Data, middle in self.iterate_data(this_file_middles):
                try:
                    this_DataSet = dm.DataSet(params, this_Data, map, 
                                       self.n_chan, self.delta_freq, self.pols,
                                       self.pol_ind, self.band_centres,
                                       self.band_ind, middle)
                except dm.DataSetError:
                    continue
                # Check that the polarization for this data is the correct
                # one.
                if this_DataSet.pol_val != pol_val:
                    raise RuntimeError("Data polarizations not consistant.")
                # Making the Noise and adding the mask is done in the
                # constructor for DataSet.
                #
                # The thermal part.
                thermal_noise = this_DataSet.get_noise_parameter("thermal")
                try:
                    this_DataSet.Noise.add_thermal(thermal_noise)
                except dm.NoiseError:
                    continue
                # Get the noise parameter database for this data.
                if params['noise_parameter_file']:
#                    noise_entry = (self.noise_params[middle]
                    noise_entry = (this_DataSet.noise_params[middle]
                        [int(round(band_centre/1e6))][pol_val])
                    # Should save this into the Data_Set so it doesn't have
                    # to be kept again there.
                else:
                    noise_entry = None

                # Frequency correlations.
                if params['frequency_correlations'] == 'mean':
                    mean_overf = this_DataSet.get_noise_parameter("mean_over_f")
                    this_DataSet.Noise.add_correlated_over_f(*mean_overf)
                elif params['frequency_correlations'] == 'measured':
                    # Correlated channel mode noise.
                    # The first frequency modes are completly deweighted. This
                    # was seen to be nessisary from looking at the noise
                    # spectra.
                    n_modes_discard = params['number_frequency_modes_discard']
                    try:
                        for ii in range(n_modes_discard):
                            mode_noise_params = this_DataSet.get_noise_parameter(
                                    "over_f_mode_" + repr(ii))
                            # Instead of completly discarding, multiply 
                            # amplitude by 2.
                            #N.deweight_freq_mode(mode_noise_params[4])
                            mode_noise_params = ((mode_noise_params[0] * 2,)
                                                 + mode_noise_params[1:])
                            this_DataSet.Noise.add_over_f_freq_mode(*mode_noise_params)
                        # The remaining frequency modes are noise weighed
                        # normally.
                        n_modes = params['number_frequency_modes']
                        for ii in range(n_modes_discard, n_modes):
                            mode_noise_params = this_DataSet.get_noise_parameter(
                                "over_f_mode_" + repr(ii))
                            this_DataSet.Noise.add_over_f_freq_mode(*mode_noise_params)
                        # All channel low frequency noise.
                        all_chan_noise_params = this_DataSet.get_noise_parameter(
                                'all_channel')
                        # In some cases the corner frequncy will be so high that
                        # this scan will contain no information.
                        this_DataSet.Noise.add_all_chan_low(*all_chan_noise_params)
                    except dm.NoiseError:
                        if self.feedback > 1:
                            print ("Noise parameters risk numerical"
                                   " instability. Data skipped.")
                        continue
                elif params['frequency_correlations'] == 'None':
                    pass
                else:
                    raise dm.ValueError("Invalid frequency correlations.")


#                # From input parameters decide what the noise model should be.
#                if params['frequency_correlations'] == 'None':
#                    model = 'thermal_only'
#                else:
#                    # This must be set to something so code doesn't say
#                    # it doesn't exist.
#                    model = None
#                # TODO Other options.
#                this_DataSet.setup_noise_from_model(model,
#                        params['deweight_time_mean'],
#                        params['deweight_time_slope'], noise_entry)
                if params['deweight_time_mean']:
                    this_DataSet.Noise.deweight_time_mean(T_huge**2)
                if params['deweight_time_slope']:
                    this_DataSet.Noise.deweight_time_slope(T_huge**2)
                # Set the time stream foregrounds for this data.
                if (params['ts_foreground_mode_file']
                    and params['n_ts_foreground_modes']) :
                    foreground_modes_db = self.foreground_modes
                    key = ts_measure.get_key(file_middle)
                    # From the database, read the data for this file.
                    db_freq = foreground_modes_db[key + '.freq']
                    db_vects = foreground_modes_db[key + '.vects']
                    ts_foregrounds_entry = (db_freq, db_vects)
                    this_DataSet.setup_ts_foregrounds(
                        params['n_ts_foreground_modes'], 
                        ts_foregrounds_entry)
                # No need for time_stream_list since data in DataSet. 
                #time_stream_list.append(this_DataSet.time_stream)
                data_set_list.append(this_DataSet)
            if self.feedback > 1:
                print "Finalizing %d noise blocks: " % len(data_set_list)
            noise_queue = Queue()
            def thread_work():
                while True:
                    thread_N = noise_queue.get()
                    #None will be the flag that there is no more work to do.
                    if thread_N is None:
                        return
                    else:
                        thread_N.Noise.finalize(frequency_correlations=(
                                     not self.uncorrelated_channels),
                                               preserve_matrices=False)
                        if self.feedback > 1:
                            sys.stdout.write('.')
                            sys.stdout.flush()
            #Start the worker threads.
            thread_list = []
            for ii in range(self.n_processes):
                T = threading.Thread(target=thread_work)
                T.start()
                thread_list.append(T)
            #Now put work on the queue for the threads to do.
            for a_DataSet in data_set_list:
                noise_queue.put(a_DataSet)
            #At the end of the queue, tell the threads that they are done.
            for ii in range(self.n_processes):
                noise_queue.put(None)
            #Wait for the threads.
            for T in thread_list:
                T.join()
            if not noise_queue.empty():
                msg = "A thread had an error in Noise finalization."
                raise RuntimeError(msg)

            '''# Finalize the noises. Each process has a subset of the total
            # noises, so no threading, just a loop.
            for a_DataSet in data_set_list:
                a_DataSet.Noise.finalize(frequency_correlations=
                                       (not self.uncorrelated_channels),
                                   preserve_matrices=False)
                if self.feedback > 1:
                    sys.stdout.write('.')
                    sys.stdout.flush()'''

            if self.feedback > 1:
                print
                print "Noise finalized."

            #temporary code to be sure all noise has time_mode_update attribute
            for data in data_set_list:
                time_mode=data.Noise.time_mode_update

            # Each process holds a subset of the DataSets. They will be
            # passed around "in a circle" with mpi so that each process
            # sees all the sets. Right now, find the "previous" and "next"
            # process that will make the above statement true.
            # Note: This ordering will make others see the package
            #       in "numerical" order.
            if ((self.rank+1) == self.nproc):
                self.prev_guy = 0
            else:
                self.prev_guy = self.rank + 1

            if (self.rank == 0):
                self.next_guy = self.nproc - 1
            else:
                self.next_guy = self.rank - 1

            # Each process will handle an equal subset of indices of
            # the cov_inv to fill.
            # This does not have to be rediscovered for every DataSet list
            # passed in, so get the indices first.
            chan_index_list = range(self.n_chan)
            dtype=np.float64
            dsize=MPI.DOUBLE.Get_size()
            if self.uncorrelated_channels:
                # Using self.rank after gets the info for
                # the approproate process.
                # 'junk' for now.
                index_list,junk = split_elems(chan_index_list,self.nproc)
                index_list = index_list[self.rank]
                # Allocate hdf5 file
                if self.rank == 0:
                    print 'rank zero' + '\n'
                    print self.n_chan
                    print self.n_ra
                    print self.n_dec
                if start_file_ind == 0:
                    data_offset, file_size = allocate_hdf5_dataset(self.cov_filename, 'inv_cov', (self.n_chan, self.n_ra, self.n_dec,                                                                self.n_ra, self.n_dec), dtype)
            else:
                ra_index_list = range(self.n_ra)
                # cross product = A x B = (a,b) for all a in A, for all b in B
                chan_ra_index_list = cross([chan_index_list,ra_index_list])
                index_list,start_list = split_elems(chan_ra_index_list,
                                                             self.nproc)
                index_list = index_list[self.rank]
                # This is the point in n_chan x n_ra that the first index
                # is at. This is needed for knowing which part of the
                # total matrix a process has.
                f_ra_start_ind = start_list[self.rank]
                # Allocate hdf5 file
                if start_file_ind == 0:
                    data_offset, file_size = allocate_hdf5_dataset(self.cov_filename, 'inv_cov', (self.n_chan, self.n_ra, self.n_dec,                                                                self.n_chan, self.n_ra, self.n_dec), dtype)
            # Wait for file to be written before continuing.
            comm.Barrier()
            print '\n' + 'Process ' + str(self.rank) + ' Passed first barrier.' + '\n'
            # Since the DataSets are split evenly over the processes,
            # have to put in the pointing/noise at each index for
            # each DataSet list that the processes hold.
            
            #for run in range(self.nproc):
            # The commented loop above is for passing subsections of the data
            # between nodes.  It currently doesn't work due to pickling issues.
            # I replaced it with a loop of just one pass to avoid re-indenting.
            for run in range(1):
                print '\n' + 'Process' + ' ' + str(self.rank) + ' ' +  ',run' + ' ' + str(run) + '\n'
                # Each DataSet has to be applied to the dirty map, too.
                # Since the '_mpi'dirty map is much smaller than the noise,
                # We will just let processor 0 do all the work
                # for the dirty map.
                # Must be done here to not pass the DataSets around twice.
                if self.rank == 0:
                    for ii in xrange(len(data_set_list)):
                        print '\n' + 'Process ' + str(self.rank) + ' adding to map data piece' + ' ' + str(ii) + ' of ' + str(len(data_set_list)) +' during run ' + str(run) + '\n'
                        #time_stream = time_stream_list[ii]
                        time_stream = data_set_list[ii].time_stream
                        N = data_set_list[ii].Noise
                        P = data_set_list[ii].Pointing
                        w_time_stream = N.weight_time_stream(time_stream)
                        map += P.apply_to_time_axis(w_time_stream)

                #If first pass, prepare thread_cov_inv_chunk to be worked on by threads
                if run == 0 and start_file_ind == 0:
                    if self.uncorrelated_channels:
                        thread_cov_inv_chunk = sp.zeros((len(index_list),                                                                                                                                self.n_ra, self.n_dec,                                                                                                                            self.n_ra, self.n_dec),                                                                                                                            dtype=dtype)
                        #Adding prior to the diagonal.
                        for ii in xrange(len(index_list)):
                            thread_cov_inv_chunk[ii,...].flat[::self.n_ra * self.n_dec + 1] += \
                                                            1.0 / T_large**2
                    else:
                       thread_cov_inv_chunk = sp.zeros((len(index_list),                                                                                                                             self.n_dec, self.n_chan,                                                                                                                          self.n_ra, self.n_dec),                                                                                                                           dtype=dtype)
                index_queue = Queue()
                def thread_work():
                    while True:
                        thread_inds = index_queue.get()
                        # None will be the flag that there is no more work to do.
                        if thread_inds is None:
                            return
                        if self.uncorrelated_channels:
                            if params['thread_divide'] == False:
                                thread_f_ind = thread_inds
                                #thread_cov_inv_block = sp.zeros((self.n_ra,                                                                                                                          self.n_dec, self.n_ra, self.n_dec), dtype=dtype)
                                for thread_D in data_set_list:
                                    thread_D.Pointing.noise_channel_to_map(thread_D.Noise,                                                                                                                 thread_f_ind + index_list[0], thread_cov_inv_chunk[thread_f_ind,...])
                                #No lock below currently.  I don't think it is necessary without the memmap.
                                #thread_cov_inv_chunk[thread_f_ind,...] += thread_cov_inv_block
                            else:
                                thread_f_ind = thread_inds[0]
                                ra_vals = thread_inds[1]
                                ra_min = ra_vals[0]
                                ra_cap = ra_vals[-1]+1
                                ra_size = ra_cap - ra_min
                                #thread_cov_inv_block = sp.zeros((ra_size,                                                                                                                          self.n_dec, self.n_ra, self.n_dec), dtype=float)
                                for thread_D in data_set_list:
                                   thread_D.Pointing.noise_channel_to_map(thread_D.Noise,                                                                                                                   thread_f_ind + index_list[0], thread_cov_inv_chunk[thread_f_ind,ra_min:ra_cap,...], ra0_range = [ra_min, ra_cap])
                                   #No lock below currently. I don't think it is necessary without the memmap.
                                #thread_cov_inv_chunk[thread_f_ind,ra_min:ra_cap,...] += thread_cov_inv_block
                        else:
                            ii = thread_inds
                            #thread_cov_inv_row = sp.zeros((self.n_dec,                                                                                                                                       self.n_chan, self.n_ra,                                                                                                                           self.n_dec),dtype=dtype)
                            thread_f_ind = index_list[ii][0]
                            thread_ra_ind = index_list[ii][1]
                            for thread_D in data_set_list:
                                thread_D.Pointing.noise_to_map_domain(thread_D.Noise,                                                                                                                                   thread_f_ind,                                                                                                                                     thread_ra_ind, thread_cov_inv_chunk[ii,...])
                            if run == 0 and start_file_ind == 0:
                                #Adding prior to the diagonal.
                                #thread_cov_inv_row_adder = thread_cov_inv_row[:, thread_f_ind, thread_ra_ind, :]
                                thread_cov_inv_chunk[ii,:, thread_f_ind, thread_ra_ind,:].flat[:: self.n_dec + 1] += \
                                                           1.0 / T_large**2
                            if (self.feedback > 1 and thread_ra_ind == self.n_ra - 1):
                                print thread_f_ind,
                                sys.stdout.flush()
                            #thread_cov_inv_chunk[ii,...] += \
                            #                      thread_cov_inv_row
                #Now start the worker threads.
                thread_list = []
                for ii in range(self.n_processes):
                    T = threading.Thread(target=thread_work)
                    T.start()
                    thread_list.append(T)
                    print '\n' + 'Process ' + str(self.rank) + ' is using ' + str(self.n_processes) + 'threads.' + '\n'
                #Now put work on the queue for the threads to do.
                for ii in xrange(len(index_list)):
                    if params['thread_divide'] == False:
                        index_queue.put(ii)
                    else:
                        if self.uncorrelated_channels: 
                            split_ra = split_elems(xrange(self.n_ra),self.n_processes)[0]
                            for split in split_ra:
                                index_queue.put((ii, split))
                        else:
                            index_queue.put(ii)
                for ii in range(self.n_processes):
                    index_queue.put(None)
                for T in thread_list:
                    T.join()
                if not index_queue.empty():
                    msg = "A thread had an error while building map covariance."
                    raise RuntimeError(msg)
                                
                # Commented section below is without threading.
                '''if self.uncorrelated_channels:
                    # No actual threading going on.
                    if run == 0 and start_file_ind == 0:
                        thread_cov_inv_chunk = sp.zeros((len(index_list),
                                                        self.n_ra, self.n_dec,
                                                        self.n_ra, self.n_dec),                                                                                                                             dtype=float)
                        #for thread_f_ind in index_list:
                        print '\n' + 'Process ' + str(self.rank) + ' just created inv_chunk' + ' during run ' + str(run) + '\n'
                    for thread_f_ind in index_list:
                    #for ii in xrange(len(index_list)):
                        thread_cov_inv_block = sp.zeros((self.n_ra, self.n_dec,
                                    self.n_ra, self.n_dec), dtype=float)
                        #if run == 0 and start_file_ind == 0:
                        #    thread_cov_inv_block.flat[::self.n_ra * self.n_dec + 1] += \
                        #                              1.0 / T_large**2

                        for thread_D in data_set_list: 
                            thread_D.Pointing.noise_channel_to_map(
                                  thread_D.Noise, thread_f_ind, thread_cov_inv_block)
                        
                        if run == 0 and start_file_ind == 0:
                            thread_cov_inv_block.flat[::self.n_ra * self.n_dec + 1] += \
                                                      1.0 / T_large**2

                        #if start_file_ind == 0:
                            # The first time through the matrix. We 'zero'
                            # everything by just assigning a value instead of
                            # setting to zero then += value.
                            # TODO: writing to same place might hiccup.
                            #cov_inv[thread_f_ind,...] = thread_cov_inv_block
                            # Add stuff to diagonal now since it's not
                            # very fun later.
                            #thread_cov_inv_block.flat \
                            #        [::self.n_ra * self.n_dec + 1] += \
                            #        1.0 / T_large**2
                        #else:
                            # Not the first time through. Just add in vals.
                            #cov_inv[thread_f_ind,...] += thread_cov_inv_block
                        #thread_cov_inv_chunk[thread_f_ind,...] += thread_cov_inv_block
                        thread_cov_inv_chunk[thread_f_ind-index_list[0],...] += thread_cov_inv_block
                        print '\n' + 'Process ' + str(self.rank) + ' added to cov_inv block at freq ' + str(thread_f_ind) + ' , using data from pass ' + str(run) + ' of ' + str(self.nproc) + '\n'  
                        if self.feedback > 1:
                            print thread_f_ind,
                            sys.stdout.flush()
                else:
                    # Keep all of a process' rows in memory, then use
                    # MPI and file views to write it out fast.
                    if run == 0 and start_file_ind == 0:
                        thread_cov_inv_chunk = sp.zeros((len(index_list),
                                                     self.n_dec, self.n_chan,
                                                     self.n_ra, self.n_dec),
                                                     dtype=float)
                    #for thread_f_ind,thread_ra_ind in index_list:
                    for ii in xrange(len(index_list)):
                        thread_cov_inv_row = sp.zeros((self.n_dec, self.n_chan,
                                                       self.n_ra, self.n_dec),
                                                      dtype=float)
                        thread_f_ind = index_list[ii][0]
                        thread_ra_ind = index_list[ii][1]
                        if run == 0 and start_file_ind == 0:
                            thread_cov_inv_row_adder = thread_cov_inv_row[:, thread_f_ind, thread_ra_ind, :]
                            thread_cov_inv_row_adder.flat[:: self.n_dec + 1] += \
                                                1.0 / T_large**2
                        for thread_D in data_set_list:
                            thread_D.Pointing.noise_to_map_domain(
                                             thread_D.Noise, thread_f_ind,
                                             thread_ra_ind, thread_cov_inv_row)


                        #writeout_filename = self.cov_filename \
                        #                  + repr(thread_inds).zfill(20)
                        #f = open(writeout_filename, 'w')
                        #np.save(f,thread_cov_inv_row)
                        #f.close()
                        if (self.feedback > 1
                            and thread_ra_ind == self.n_ra - 1):
                            print thread_f_ind,
                            sys.stdout.flush()
                        thread_cov_inv_chunk[ii,...] += \
                                                 thread_cov_inv_row
                # Once done with one list of DataSets, do next.
                # Note: No need to pass anything the last time since
                #       that data was the process' original data and
                #       has already been processed.
                '''
                #if (run != (self.nproc - 1)):

                #Replace if statement below with if statement above to pass data between nodes.
                if run == 1:
                    # Send out the DataSet list I have to next guy.
                    print '\n' + 'Process ' + str(self.rank) + ' is passing data_set_list to process ' + str(self.next_guy) + ', at the end of run ' + str(run) + '\n'
                    #Pickling DataSet objects in data_set_list 
                    def pickle_list(list):
                        pickled_list = []
                        for item in list:
                            pickled_list.append(cPickle.dumps(item))
                        return pickled_list
                    data_set_list_pickled = pickle_list(data_set_list)
                    comm.send(data_set_list_pickled,dest=self.next_guy)
                    # Receive the DataSet list being sent to me by prev guy.
                    data_set_list_pickled = comm.recv(source=self.prev_guy)
                    #Unpickling DataSet objects
                    def unpickle_list(list):
                        unpickled_list = []
                        for item in list:
                            unpickled_list.append(cPickle.loads(item))
                        return unpickled_list
                    data_set_list=unpickle_list(data_set_list_pickled)
                    ## Same for time_stream_list.
                    #comm.send(time_stream_list,dest=self.next_guy)
                    #time_stream_list = comm.recv(source=self.prev_guy)

                    # Processes will be relatively in sync here since they
                    # will block on read until the "previous" process
                    # sends over its data.

            # Now that thread_cov_inv_chunk is all filled with
            # data from all DataSets, write it out.
            # Only dealing with correlated channels now.
        if not self.uncorrelated_channels:
                #total_shape = (self.n_chan*self.n_ra, self.n_dec,
                             #  self.n_chan, self.n_ra, self.n_dec)
                #start_ind = (f_ra_start_ind,0,0,0,0)
            lock_and_write_buffer(thread_cov_inv_chunk, self.cov_filename, data_offset + dsize*f_ra_start_ind*self.n_dec*self.n_chan*self.n_ra*self.n_dec, dsize*thread_cov_inv_chunk.size, self.rank)
            if self.rank == 0:
                f = h5py.File(self.cov_filename, 'r+')
                cov_inv_dset = f['inv_cov']
                attrs = cov_inv_dset.attrs
                attrs.__setitem__('rows', str((0, 1, 2)))
                attrs.__setitem__('cols', str((3, 4, 5)))
                attrs.__setitem__('dec_centre', str(self.params['field_centre'][1]))
                attrs.__setitem__('ra_centre', str(self.params['field_centre'][0]))
                attrs.__setitem__('type', "'mat'")
                attrs.__setitem__('axes', "('freq', 'ra', 'dec', 'freq', 'ra', 'dec')")
                attrs.__setitem__('freq_centre', str(band_centre))
                attrs.__setitem__('freq_delta', str(self.delta_freq))
                attrs.__setitem__('ra_delta', str(self.ra_spacing))
                attrs.__setitem__('dec_delta', str(self.params['pixel_spacing']))
                #cov_inv_shape = sp.zeros(shape = (self.n_chan, self.n_ra, self.n_dec, self.n_chan, self.n_ra, self.n_dec), dtype = dtype)
                #cov_inv_info = al.make_mat(cov_inv_shape, axis_names=('freq', 'ra', 'dec','freq', 'ra', 'dec'),                                                                                                           row_axes=(0, 1, 2), col_axes=(3, 4, 5))
                #cov_inv_info.copy_axis_info(map)
                #for key, value in cov_inv_info.info.iteritems():
                #    cov_inv_dset.attrs[key] = repr(value)
                # NOTE: using 'float' is not supprted in the saving because
                # it has to know if it is 32 or 64 bits.
                #dtype = thread_cov_inv_chunk.dtype
                #dtype = np.float64 # default for now
                #np.save(self.cov_filename+'_mpi_corr_proc'+str(self.rank),thread_cov_inv_chunk)
                # Save array.
                #mpi_writearray(self.cov_filename+'_mpi', thread_cov_inv_chunk,
                               #comm, total_shape, start_ind, dtype,
                               #order='C', displacement=0)
   

        if self.uncorrelated_channels:
                #total_shape = (self.n_chan, self.n_ra,
                             #  self.n_dec, self.n_ra, self.n_dec)
                #start_ind = (index_list[0],0,0,0,0)
            lock_and_write_buffer(thread_cov_inv_chunk, self.cov_filename, data_offset + dsize*index_list[0]*self.n_dec*self.n_ra*self.n_dec*self.n_ra, dsize*thread_cov_inv_chunk.size, self.rank)
            if self.rank == 0:
                f = h5py.File(self.cov_filename, 'r+')
                cov_inv_dset = f['inv_cov']
                attrs = cov_inv_dset.attrs
                attrs.__setitem__('rows', str((0, 1, 2)))
                attrs.__setitem__('cols', str((0, 3, 4)))
                attrs.__setitem__('dec_centre', str(self.params['field_centre'][1]))
                attrs.__setitem__('ra_centre', str(self.params['field_centre'][0]))
                attrs.__setitem__('type', "'mat'")
                attrs.__setitem__('axes', "('freq', 'ra', 'dec', 'ra', 'dec')")
                attrs.__setitem__('freq_centre', str(band_centre))
                attrs.__setitem__('freq_delta', str(self.delta_freq))
                attrs.__setitem__('ra_delta', str(self.ra_spacing))
                attrs.__setitem__('dec_delta', str(self.params['pixel_spacing']))
                #cov_inv = al.make_mat(cov_inv_dset,                                                                                                                 axis_names=('freq', 'ra', 'dec', 'ra', 'dec'),                                                                                  row_axes=(0, 1, 2), col_axes=(0, 3, 4))
                #cov_inv_shape = sp.zeros(shape = (self.n_chan, self.n_ra, self.n_dec, self.n_ra, self.n_dec), dtype = dtype)
                #cov_inv_shape = _mapmaker_c.large_empty((self.n_chan, self.n_ra, self.n_dec, self.n_ra, self.n_dec))
                #cov_inv_info = al.make_mat(cov_inv_shape,                                                                                                                      axis_names=('freq', 'ra', 'dec', 'ra', 'dec'),                                                                                       row_axes=(0, 1, 2), col_axes=(0, 3, 4))
                #cov_inv_info.copy_axis_info(map)
                #for key, value in cov_inv_info.info.iteritems():
                #    cov_inv_dset.attrs[key] = repr(value)
                # NOTE: using 'float' is not supprted in the saving because
                # it has to know if it is 32 or 64 bits.
                #dtype = thread_cov_inv_chunk.dtype
                #dtype = np.float64 # default for now
                #np.save(self.cov_filename+'_mpi_uncorr_proc'+str(self.rank),thread_cov_inv_chunk)
                # Save array.
                #mpi_writearray(self.cov_filename+'_mpi', thread_cov_inv_chunk,
                               #comm, total_shape, start_ind, dtype,                                                                                                              #order='C', displacement=0)                    



# Close self.noise_params in each DataSet.
# Along those lines, might want to open and close self.foreground_modes, too.
        #for a_data in data_set_list:
        #    if not a_data.noise_params is None:
        #        a_data.noise_params.close()
        #    if not a_data.foreground_modes is None:
        #        a_data.foreground_modes.close(




 


#### Utilities ####

def split_elems(points, num_splits):
    '''Split a list of arbitrary type elements, points, into num_splits sets.
    Also, returns a list of the index that each set starts at.
    Only tested with 'points' being:
    list of str
    list of int
    list of (tuple of int)'''

    size = len(points)
    start_list = []
    big_list = []
    for i in range(num_splits):
        big_list.append([])
    # Number of elements per split/process.
    divdiv = size / num_splits
    # Number of processes that get one extra element.
    modmod = size % num_splits
    # Put in the right number of elements for each process. The 'extra'
    # elements from mod get put into the first processes, so the difference
    # in the number of elements per process is either 1 or 0.
    start = 0
    end = divdiv
    for i in range(num_splits):
        if i < modmod:
            end += 1
        for j in range(start,end):
            if type(points[j]) == str:
                big_list[i].append(points[j])
            else:
                try:
                    big_list[i].append(tuple(points[j]))
                except TypeError:
                    big_list[i].append(points[j])
        start_list.append(start)
        start = end
        end += divdiv
    return big_list, start_list


def cross(set_list):
        '''Given a list of sets, return the cross product.'''
        # By associativity of cross product, cross the first two sets together
        # then cross that with the rest. The big conditional in the list
        # comprehension is just to make sure that there are no nested lists
        # in the final answer.
        if len(set_list) == 1:
                ans = []
                for elem in set_list[0]:
                        # In the 1D case, these elements are not lists.
                        if type(elem) == list:
                                ans.append(np.array(elem))
                        else:
                                ans.append(np.array([elem]))
                return ans
        else:
                A = set_list[0]
                B = set_list[1]
                cross_2 = [a+b if ((type(a) == list) and (type(b) == list)) else \
                (a+[b] if type(a) == list else ([a]+b if type(b) == list else \
                [a]+[b])) for a in A for b in B]
                remaining = set_list[2:]
                remaining.insert(0,cross_2)
                return cross(remaining)

@profile
def lock_and_write_buffer(obj, fname, offset, size, proc):
    """Write the contents of a buffer to disk at a given offset, and explicitly
    lock the region of the file whilst doing so.

    Parameters
    ----------
    obj : buffer
        Data to write to disk.
    fname : string
        Filename to write.
    offset : integer
        Offset into the file to start writing at.
    size : integer
        Size of the region to write to (and lock).
    """
    import os
    #import os.fcntl as fcntl
    import fcntl
    import h5py

    buf = buffer(obj)

    if len(buf) > size:
        raise Exception("Size doesn't match array length.")

    fd = os.open(fname, os.O_RDWR | os.O_CREAT)
    #fd = open(fname, 'rw')
    #fd = h5py.File(fname, 'r+')

    try:
        fcntl.lockf(fd, fcntl.LOCK_EX, size, offset, os.SEEK_SET)
    except:
        print "Could not obtain lock"

    os.lseek(fd, offset, 0)
    
    a=0
    while(True):
    #while(len(buf)>0):
        #nb = os.write(fd,buf)
        nb = os.write(fd, buf[a:])
        print "The buffer save started at byte " + str(a) + " for process " + str(proc)
        if nb < 0:
            raise Exception("Failed write")

        if nb == len(buf[a:]):
            break
        else:
            a += nb
        #buf = buf[nb:]

    '''
    nb = os.write(fd, buf)

    if nb != len(buf):
        #raise Exception("Something funny happened with the reading.")
        raise Exception("Something funny happened with the reading." + "  The buffer is length " + str(len(buf)) + " but the number of bytes written was " + str(nb) + " for process " + str(proc) + "." ) '''
    try:
        fcntl.lockf(fd, fcntl.LOCK_UN)
    except:
        print "Could not unlock"

    os.close(fd)


def allocate_hdf5_dataset(fname, dsetname, shape, dtype, comm=MPI.COMM_WORLD):
    """Create a hdf5 dataset and return its offset and size.

    The dataset will be created contiguously and immediately allocated,
    however it will not be filled.

    Parameters
    ----------
    fname : string
        Name of the file to write.
    dsetname : string
        Name of the dataset to write (must be at root level).
    shape : tuple
        Shape of the dataset.
    dtype : numpy datatype
        Type of the dataset.
    comm : MPI communicator
        Communicator over which to broadcast results.

    Returns
    -------
    offset : integer
        Offset into the file at which the dataset starts (in bytes).
    size : integer
        Size of the dataset in bytes.

    """

    import h5py

    state = None

    if comm.rank == 0:

        # Create/open file
        f = h5py.File(fname, 'a')

        # Create dataspace and HDF5 datatype
        sp = h5py.h5s.create_simple(shape, shape)
        tp = h5py.h5t.py_create(dtype)

        # Create a new plist and tell it to allocate the space for dataset
        # immediately, but don't fill the file with zeros.
        plist = h5py.h5p.create(h5py.h5p.DATASET_CREATE)
        plist.set_alloc_time(h5py.h5d.ALLOC_TIME_EARLY)
        plist.set_fill_time(h5py.h5d.FILL_TIME_NEVER)

        # Create the dataset
        dset = h5py.h5d.create(f.id, dsetname, tp, sp, plist)

        # Get the offset of the dataset into the file.
        state = dset.get_offset(), dset.get_storage_size()

        f.close()

    state = comm.bcast(state, root=0)

    return state




# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        par_file = sys.argv[1]
        nproc = 1
        DirtyMapMaker(par_file).execute(nproc)
    elif len(sys.argv) == 3:
        par_file = sys.argv[1]
        nproc = int(sys.argv[2])
        DirtyMapMaker(par_file).execute(nproc)
    else:
        print "Usage: `python map/dirty_map.py parameter_file [num_threads]"

