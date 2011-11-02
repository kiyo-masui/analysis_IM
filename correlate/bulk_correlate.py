import sys
import scipy as sp
from utils import file_tools as ft
from kiyopy import parse_ini
import kiyopy.utils
from utils import data_paths as dp
# using the prescription for which noise to associate to the map sets form all
# pairs of maps for the correlation code.:s
# GBT_15hr_cleaned_Eric_0mode A_with_B;map

# calculate correlations in batch (designed to be used with PBS/scinet)
# problem: each correlation wants a different .ini file.
# can group correlations and use multiprocessing, but would rather have pbs
# spin one-off correlations to run on scinet
# standalone function that 

# want: A_with_B_x_B_with_A -- these six pairs
# split on cross_sym, _with_ A B

# need an intelligible way to summarize correlations:
# for each map: field name, {15hr, wigglez, etc} method it was cleaned with
# A_with_B, etc.

# map1_key_interiorname_noise1_key_interiorname
# GBT_15hr_cleaned_Eric_0mode_A_with_B_x_GBT_15hr_cleaned_Eric_0mode_B_with_A.ini

params_init = {
               'output_root': "./data_test/",
               'ini_root': "./data_test/",
               'output_filename': "corr.pkl",
               'pairlist': [{"map1": 'GBT_15hr_cleaned_Eric_0mode',
                              "map2": 'GBT_15hr_cleaned_Eric_0mode',
                              "noise_inv1": 'GBT_15hr_cleaned_Eric_0mode',
                              "noise_inv2": 'GBT_15hr_cleaned_Eric_0mode'}],
               'nfreq_bin': 256,
               'freq_list': (),
               # Angular lags at which to calculate the correlation.  Upper
               # edge bins in degrees.
               'lags': "tuple(sp.arange(0.002, 0.2, 0.12))"
               }
prefix = 'fs_'


class BulkCorrelate(ft.ClassPersistence):
    def __init__(self, parameter_file_or_dict=None):
        self.params = parse_ini.parse(parameter_file_or_dict, params_init,
                                      prefix=prefix)

        self.freq_list = sp.array(self.params['freq_list'], dtype=int)
        self.lags = self.params['lags']
        self.nfreq_bin = self.params['nfreq_bin']

        #self.output_root = self.datapath_db.fetch(self.params['output_root'],
        #                                          intend_write=True)
        self.output_root = self.params['output_root']
        self.ini_root = self.params['ini_root']

        # Write parameter file.
        kiyopy.utils.mkparents(self.ini_root)
        parse_ini.write_params(self.params, self.ini_root + 'params.ini',
                               prefix=prefix)

    def execute(self):
        for pairitem in self.params['pairlist']:
            print pairitem
            self.define_map_pairs(pairitem['map1'], pairitem['map2'],
                                  pairitem['noise_inv1'], pairitem['noise_inv2'])

    def define_map_pairs(self, map1_key, map2_key, noise1_key, noise2_key):
        par = self.params
        (self.pairlist, self.pairdict) = \
                dp.cross_maps(map1_key, map2_key,
                              noise1_key, noise2_key,
                              map_suffix=";map",
                              noise_inv_suffix=";noise_inv",
                              cross_sym="_x_",
                              pair_former="GBTauto_cross_pairs",
                              ignore=['param'],
                              tag1prefix = map1_key + "_",
                              tag2prefix = map2_key + "_",
                              verbose=False)

        freq_list_full = range(self.nfreq_bin)
        for pairname in self.pairlist:
            ininame = pairname + ".ini"
            fh = open(self.ini_root + ininame, 'w')
            print ininame
            fileinfo = self.pairdict[pairname]
            fh.write("import os\nimport scipy as sp\n\n")
            fh.write("fs_output_root = '%s'\n" % self.params['output_root'])
            fh.write("fs_output_filetag = '%s'\n" % pairname)
            fh.write("fs_map1 = '%s'\n" %fileinfo['map1'])
            fh.write("fs_noise_inv1 = '%s'\n" % fileinfo['noise_inv1'])
            fh.write("fs_map2 = '%s'\n" % fileinfo['map2'])
            fh.write("fs_noise_inv2 = '%s'\n" % fileinfo['noise_inv2'])
            fh.write("fs_lags = %s\n" % self.lags)
            fh.write("fs_freq_list = range(%d)\n" % self.nfreq_bin)
            cutlist = []
            for fbin in range(self.nfreq_bin):
                if fbin not in self.freq_list:
                    cutlist.append(fbin)
            fh.write("cutlist = %s\n" % repr(cutlist))
            comp = '[x for x in fs_freq_list if x not in cutlist]'
            fh.write("fs_freq_list = tuple(%s)\n" % comp)
            fh.close()

if __name__ == '__main__':
    if len(sys.argv) == 2:
        BulkCorrelate(str(sys.argv[1])).execute()
    else:
        print 'Need one argument: parameter file name.'

