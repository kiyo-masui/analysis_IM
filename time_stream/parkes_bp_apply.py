#Script for dividing Parkes timestream data by a bandpass estimate.
#The timestream data is in the form of GBT style data blocks.
#Each scan should have 13 blocks, for beams 1-13.
#The bandpass estimate is a .npy file.  It can also include 
#a calibration factor.

import numpy as np
from core import fitsGBT
from kiyopy import parse_ini
import os

params_init = {'input_root': '/scratch2/p/pen/andersoc/second_parkes_pipe/rebinned/SDFITS/',
            'file_middles': ['2df'],
            'input_end':  '.fits',
            'output_root': '/scratch2/p/pen/andersoc/second_parkes_pipe/rebinned/SDFITS/bpdivide/',
            'output_middles': ['2df'],
            'bp_root': '/scratch2/p/pen/andersoc/second_parkes_pipe/bp_estimates/svd_estimates/correct_effic/',
            'bp_yy_fname': 'bandpasses_[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12]_ra165_dec0_YY.npy',
            'bp_xx_fname': 'bandpasses_[1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13]_ra165_dec0_XX.npy',
             'split_beams': True}

prefix = 'ts_'

def bp_divide(data_blocks, x_bp_array, y_bp_array, x_beams = [1,2,3,4,5,6,7,8,9,11,12,13], y_beams = [1,2,3,4,5,6,7,8,9,10,12]):
    #bp_array files should have dimensions (beam, freq).
    #The beams arrays tell the order of the beams in the bp_array.
    #Pol is field['CRVAL4'].  -5 is XX, -6 is YY.
    #In data matrix, pol is the second index.  Index 0 is XX, 1 is YY.
    if len(data_blocks) != 13:
        raise Exception("Data file does not have 13 data blocks.")
    for i in range(len(x_beams)):
        if data_blocks[x_beams[i]-1].field['BEAM'][()] != x_beams[i]:
            raise Exception("Beam of data block does not match bandpass beam.")
        data_blocks[x_beams[i]-1].data[:,0,0,:] /= x_bp_array[i,:]
    for i in range(len(y_beams)):
        if data_blocks[y_beams[i]-1].field['BEAM'][()] != y_beams[i]:
            raise Exception("Beam of data block does not match bandpass beam.")
        data_blocks[y_beams[i]-1].data[:,1,0,:] /= y_bp_array[i,:]
    return data_blocks

def get_blocks(fname):
    Reader = fitsGBT.Reader(fname)
    Blocks = Reader.read((), 0)
    return Blocks

class bp_divider():

    def __init__(self, parameter_file_or_dict = None):
        self.params = params_init
        if parameter_file_or_dict:
            self.params = parse_ini.parse(parameter_file_or_dict, params_init,
                                          prefix=prefix)
        print self.params

    def execute(self):
        params = self.params
        xx_bp = np.load(params['bp_root'] + params['bp_xx_fname'])
        yy_bp = np.load(params['bp_root'] + params['bp_yy_fname'])
        for names in zip(params['file_middles'],params['output_middles']):
            blocks = get_blocks(params['input_root'] + names[0] + params['input_end'])
            blocks = bp_divide(blocks, xx_bp, yy_bp)
            if not params['split_beams']:
                writer = fitsGBT.Writer()
                writer.add_data(blocks)
                out_fname = params['output_root'] + names[1] + params['input_end']
                out_dir = out_fname[:out_fname.rfind("/")+1]
                print 'output directory is ' + out_dir
                try:
                    os.makedirs(out_dir)
                except OSError:
                    if not os.path.isdir(out_dir):
                        raise Exception("Output path doesn't exist and can't be made.")
                print 'output fname is ' + out_fname
                writer.write(out_fname)
            else:
                for block in blocks:
                    beam_num = str(block.field['BEAM'])
                    writer = fitsGBT.Writer()
                    writer.add_data(block)
                    out_fname = params['output_root'] + names[1] + '_beam' + beam_num + params['input_end']
                    out_dir = out_fname[:out_fname.rfind("/")+1]
                    try:
                        os.makedirs(out_dir)
                    except OSError:
                        if not os.path.isdir(out_dir):
                            raise Exception("Output path doesn't exist and can't be made.")
                    writer.write(out_fname)

if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        par_file = sys.argv[1]
        bp_divider(par_file).execute()
