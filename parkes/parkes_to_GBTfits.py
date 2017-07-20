"""
This program sdfits files from the Parkes format to the GBT format.

This program is pipeline compatible.
"""

import os
import datetime
import multiprocessing as mp

import numpy as np
import pyfits

from kiyopy import parse_ini, utils
import kiyopy.pickle_method
from core import data_block
from core import fitsGBT

params_init = {
               # IO:
               'input_root' : '/home/p/pen/ycli/scratch/data/pks/2012/RAWDATA_SDF/20121024/',
               # The unique part of every fname
               'file_middles' : ("testfile_GBTfits",),
               'input_files' : {"testfile_GBTfits":["2012-10-24_1505-P641_west1_1315_P641.sdfits",]},
               'output_root' : "./",
               'output_end' : ".testout.fits",
               # What data to process within each file.
               'beams' : (1,2,3,4,5,6,7),
               'table_jy2k' : "", 
               'table_ginv' : "",
               'regular_shape' : (1222, 1, 1, 2, 1024),
               }

prefix = 'p2g_'

class Converter(object):

    def __init__(self, parameter_file_or_dict, feedback=2):
        self.feedback = feedback
        
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init,
                                      prefix=prefix, feedback=feedback)

    def execute(self, n_processes=1):
        params = self.params
        # Make parent directories if need be.
        print params['output_root']
        utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)

        self.get_cal()

        n_new = n_processes - 1
        n_files = len(params['file_middles'])
        # Loop over the files to process.
        if n_new <= 0 :
            # Single process mode.
            for file_ind in range(n_files) :
                self.process_file(file_ind)
        elif n_new > 32 :
            raise ValueError("Asked for a rediculouse number of processes: " +
                             str(n_new) + ".  Limit is 32.")
        else :
            # Spawn a bunch of new processes each with a single file to
            # analyse.
            # Can't us an mp.Pool here because we don't want to reuse processes
            # due to pyfits memory leak.
            process_list = range(n_new)
            for ii in xrange(n_files + n_new) :
                if ii >= n_new :
                    process_list[ii%n_new].join()
                    if process_list[ii%n_new].exitcode != 0 : 
                        raise RuntimeError("A thread failed with exit code: "
                                        + str(process_list[ii%n_new].exitcode))
                if ii < n_files :
                    process_list[ii%n_new] = mp.Process(
                        target=self.process_file, args=(ii,))
                    process_list[ii%n_new].start()

    def process_file(self, file_ind):

        params = self.params
        file_middle = params['file_middles'][file_ind]
        output_fname = (params['output_root']
                        + file_middle + params['output_end'])

        input_files = params['input_files'][file_middle]
        Blocks = []
        for input_file in input_files:
            input_fname = (params['input_root'] + input_file)

            if self.feedback > 1:
                print "Reading Parkes fits file: " + input_fname
            hdulist = pyfits.open(input_fname)
            #if hdulist[1].data['DATA'].shape != params['regular_shape']:
            #    print 'DATA Ignore:\n' +\
            #          '\tdo not have the regular data shape\n' +\
            #          '\t',hdulist[1].data['DATA'].shape
            #    continue
            # Convert to GBT format.
            Blocks.append(self.blocks_from_parkes_hdulist(hdulist))
        # Write out.
        if Blocks == []:
            return
        beams = params['beams']
        if not beams:
            beams = range(len(Blocks)//len(input_files))
        for beam in range(len(beams)):
            Writer = fitsGBT.Writer(feedback=self.feedback)
            output_fname = (params['output_root'] + file_middle + 
                    '_beam_%02d'%beams[beam] + params['output_end'])
            Blocks_beam = []
            for Blocks_file in Blocks:
                Blocks_beam = Blocks_file[beam]
                print Blocks_beam
                Writer.add_data(Blocks_beam)
            utils.mkparents(output_fname)        
            Writer.write(output_fname)

    def blocks_from_parkes_hdulist(self, hdulist):

        hdudata = hdulist[1]
        fheader = hdudata.header
        fdata = hdudata.data
        # Get some basic info about the data.
        data_entry = hdudata.columns.names.index('DATA')
        data_dims = eval(hdudata.columns[data_entry].dim)
        n_chan = data_dims[0]
        n_pol = data_dims[1]
        # Assume that if there are 2 polarizations, they are XX and YY.
        if (n_pol == 2 and fheader['CRVAL2'] == -5 and fheader['CRPIX2'] == 1
            and fheader['CDELT2'] == -1):
            pol_codes = [-5, -6]
        else:
            raise NotImplementedError('polarizations not in expected order.')
        # First figure out the beams, make sure that the cycle regulairily etc.
        beam_axis = fdata.field('BEAM').copy()
        n_beam = 13 # Generalize later.
        n_time = len(fdata) // n_beam
        # Make sure the beam index cycles regularly.
        #beam_axis.shape = (n_time, n_beam)
        if not np.all(beam_axis.reshape(n_time, n_beam) == range(1, n_beam + 1)):
            raise NotImplementedError("Beams out of order.")
        # Figure out which beams to process.
        beams = self.params['beams']
        if not beams:
            beams = range(n_beam)

        # Loop thorugh the beams and get a data_block for each one.
        # Initialze list of output blocks.
        out_Blocks = []
        for beam_index in beams:
            this_beam_inds = np.arange(n_time) * n_beam + beam_index
            this_fdata = fdata[this_beam_inds]
            #this_fdata = fdata[:][beam_index]
            if not np.all(this_fdata['BEAM'] == beam_index + 1):
                print this_fdata['BEAM']
                raise RuntimeError('Something went horribly wrong.')
            data = this_fdata.field('DATA').copy()
            mask = this_fdata.field('FLAGGED').copy()
            # Reshape to the GBTfits order.
            data.shape = mask.shape = (n_time, n_pol, 1, n_chan)
            # Cal and convert to K
            if self.jy2k != None:
                data *= self.jy2k[beam_index][None, :, None, None]
            if self.ginv != None:
                data *= self.ginv[beam_index][None, :, None, None]
            # Create the GBT DataBlock object.
            this_Data = data_block.DataBlock(data)
            # Apply the mask. mask == 1 is nessisary to recast to bool.
            this_Data.data[mask == 1] = None
            # Fill in the time axis.
            this_Data.set_field('DATE-OBS', np.empty(n_time, dtype='S23'),
                                ('time',), '23A')
            for ii in xrange(n_time):
                date = this_fdata['DATE-OBS'][ii]
                date = datetime.datetime.strptime(date, "%Y-%m-%d")
                time = this_fdata['TIME'][ii]
                time = datetime.timedelta(seconds=time)
                date_obs_obj = date + time
                milliseconds = "%06d" % date_obs_obj.microsecond
                milliseconds = milliseconds[:3]
                date_obs = (date_obs_obj.strftime("%Y-%m-%dT%H:%M:%S.") 
                            + milliseconds)
                this_Data.field['DATE-OBS'][ii] = date_obs
            # Now fill in other fields that might be usefull.
            # These fields are functions of time.
            t_sys = np.array(this_fdata['TSYS'], dtype=np.float64)
            this_Data.set_field('TSYS', t_sys, ('time', 'pol'), '1D') 
            copy_time_field_float(this_fdata, this_Data, 'CRVAL3', 'RA')
            copy_time_field_float(this_fdata, this_Data, 'CRVAL4', 'DEC')
            copy_time_field_float(this_fdata, this_Data, 'AZIMUTH', 'AZIMUTH')
            copy_time_field_float(this_fdata, this_Data, 'ELEVATIO', 'ELEVATIO')
            copy_time_field_float(this_fdata, this_Data, 'PARANGLE', 'PARANGLE')
            # These fields are constants.
            this_Data.set_field('CRVAL4', pol_codes, ('pol',), '1I')
            this_Data.set_field('CAL', ['A'], ('cal',), '1A')
            # For some reason this is formated as a float in Parkes data.
            crpix1 = round(this_fdata['CRPIX1'][0])
            if not np.allclose(this_fdata['CRPIX1'], crpix1):
                raise RuntimeError('CRPIX1 issue.')
            # Set the freq increase
            if this_fdata['CDELT1'][0] < 0:
                print "Note: change the delta frequency to positave"
                this_Data.set_field('CDELT1', -this_fdata['CDELT1'][0], (), '1D')
                this_Data.set_field('CRVAL1',  this_fdata['CRVAL1'][0], (), '1D')
                this_Data.data = this_Data.data[:, :, :, ::-1]
            #if this_fdata['CDELT1'][0] < 0:
            #    this_Data.set_field('CDELT1', -this_fdata['CDELT1'][0], (), '1D')
            #    freq = (np.arange(data.shape[-1], dtype=float) + 1.0 -
            #            this_fdata['CRPIX1'][0]) * (-this_fdata['CDELT1'][0]) + \
            #            this_fdata['CRVAL1'][0]
            #    this_Data.data = this_Data.data[:, :, :, ::-1]
            #    this_Data.data[:, :, :, 1:] = this_Data.data[:, :, :, :-1]
            #    this_Data.data[:, :, :, 0] = 0.
            #    this_Data.set_field('CRVAL1', freq[this_fdata['CRPIX1'][0] - 1], 
            #            (), '1D')
            #    print "invserse frequency index: centre ", \
            #            freq[this_fdata['CRPIX1'][0] - 1]
            else:
                this_Data.set_field('CDELT1',  this_fdata['CDELT1'][0], (), '1D')
                this_Data.set_field('CRVAL1',  this_fdata['CRVAL1'][0], (), '1D')
            this_Data.set_field('CRPIX1', crpix1, (), '1I')

            copy_constant_field(this_fdata, this_Data, 'BEAM', 'BEAM', '1I')
            copy_constant_field(this_fdata, this_Data, 'SCAN', 'SCAN', '1I')
            #copy_constant_field(this_fdata, this_Data, 'CRVAL1', 'CRVAL1', '1D')
            #copy_constant_field(this_fdata, this_Data, 'CDELT1', 'CDELT1', '1D')
            copy_constant_field(this_fdata, this_Data, 'BANDWID', 'BANDWID',
                                '1D')
            copy_constant_field(this_fdata, this_Data, 'EXPOSURE', 'EXPOSURE',
                                '1D')
            # Add a history entry, should make this more descriptive.
            this_Data.add_history('Converted from parkes SDfits data.')

            this_Data.verify()
            out_Blocks.append(this_Data)

        return out_Blocks

    def get_cal(self):
        params = self.params

        if os.path.exists(params['table_jy2k']):
            self.jy2k = np.loadtxt(params['table_jy2k']).T
        else:
            self.jy2k = None
            print "no jy to k table"
        if os.path.exists(params['table_ginv']):
            self.ginv = np.loadtxt(params['table_ginv']).T
        else:
            self.ginv = None
            print "no gian, set to 1"

def copy_time_field_float(fdata, Data, old_key, new_key):
    '''
    field = np.array(fdata[old_key], dtype=np.float64)
    Data.set_field(new_key, field, ('time',), '1D')
    '''

    if old_key == 'CRVAL3':
        print "convert to [-180, 180]"
        ra_max = 180
        field = np.array(fdata[old_key], dtype=np.float64)
        ra_to_renorm = (field > ra_max)
        field[ra_to_renorm] = field[ra_to_renorm] -360
    else:
        field = np.array(fdata[old_key], dtype=np.float64)
    Data.set_field(new_key, field, ('time',), '1D')

def copy_constant_field(fdata, Data, old_key, new_key, format):
    field = fdata[old_key][0]
    if not np.allclose(fdata[old_key], field):
        raise RuntimeError('Expected %s (%s) field to be constant/'
                           % (old_key, new_key))
    Data.set_field(new_key, field, (), format)


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Converter(None).execute()



