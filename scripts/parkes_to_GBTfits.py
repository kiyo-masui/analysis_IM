"""
This program sdfits files from the Parkes format to the GBT format.

This program is pipeline compatible.
"""

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
               'input_root' : './testdata/',
               # The unique part of every fname
               'file_middles' : ("testfile_GBTfits",),
               'input_end' : ".fits",
               'output_root' : "./",
               'output_end' : ".testout.fits",
               # What data to process within each file.
               'beams' : ()
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
        utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
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
        input_fname = (params['input_root'] + file_middle +
                       params['input_end'])
        output_fname = (params['output_root']
                        + file_middle + params['output_end'])
        Writer = fitsGBT.Writer(feedback=self.feedback)

        if self.feedback > 1:
            print "Reading Parkes fits file: " + input_fname
        hdulist = pyfits.open(input_fname)
        # Convert to GBT format.
        Blocks = self.blocks_from_parkes_hdulist(hdulist)
        # Write out.
        Writer.add_data(Blocks)
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
        beam_axis = fdata.field('BEAM')
        n_beam = 13 # Generalize later.
        n_time = len(fdata) // n_beam
        # Make sure the beam index cycles regularly.
        beam_axis.shape = (n_time, n_beam)
        if not np.all(beam_axis == range(1, n_beam + 1)):
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
            if not np.all(this_fdata['BEAM'] == beam_index + 1):
                raise RuntimeError('Something went horribly wrong.')
            data = this_fdata.field('DATA').copy()
            mask = this_fdata.field('FLAGGED').copy()
            # Reshape to the GBTfits order.
            data.shape = mask.shape = (n_time, n_pol, 1, n_chan)
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
            copy_time_field_float(this_fdata, this_Data, 'AZIMUTH', 'CRVAL2')
            copy_time_field_float(this_fdata, this_Data, 'ELEVATIO', 'CRVAL3')
            copy_time_field_float(this_fdata, this_Data, 'PARANGLE', 'PARANGLE')
            # These fields are constants.
            this_Data.set_field('CRVAL4', pol_codes, ('pol',), '1I')
            this_Data.set_field('CAL', ['A'], ('cal',), '1A')
            # For some reason this is formated as a float in Parkes data.
            crpix1 = round(this_fdata['CRPIX1'][0])
            if not np.allclose(this_fdata['CRPIX1'], crpix1):
                raise RuntimeError('CRPIX1 issue.')
            this_Data.set_field('CRPIX1', crpix1, (), '1I')
            copy_constant_field(this_fdata, this_Data, 'BEAM', 'BEAM', '1I')
            copy_constant_field(this_fdata, this_Data, 'SCAN', 'SCAN', '1I')
            copy_constant_field(this_fdata, this_Data, 'CRVAL1', 'CRVAL1', '1D')
            copy_constant_field(this_fdata, this_Data, 'CDELT1', 'CDELT1', '1D')
            copy_constant_field(this_fdata, this_Data, 'BANDWID', 'BANDWID',
                                '1D')
            copy_constant_field(this_fdata, this_Data, 'EXPOSURE', 'EXPOSURE',
                                '1D')
            # Add a history entry, should make this more descriptive.
            this_Data.add_history('Converted from parkes SDfits data.')

            this_Data.verify()
            out_Blocks.append(this_Data)

        return out_Blocks

def copy_time_field_float(fdata, Data, old_key, new_key):
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
    Converter(str(sys.argv[1])).execute()



