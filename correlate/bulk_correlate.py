r"""prepare a correlation analysis of all relevant pairs of maps"""
import sys
import scipy as sp
from utils import file_tools as ft
from kiyopy import parse_ini
import kiyopy.utils
from utils import data_paths as dp

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
    r"""produce a set of map pairs implied by the database keys"""
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
        r"""perform the main function of this class (for pipeline call)"""
        for pairitem in self.params['pairlist']:
            print pairitem
            self.generate_ini_batch(pairitem['map1'], pairitem['map2'],
                                    pairitem['noise_inv1'],
                                    pairitem['noise_inv2'])

    def generate_ini_batch(self, map1_key, map2_key, noise1_key, noise2_key):
        r""" produce a batch of ini files which can be fed to PBS or run
        independently on prawn. The correlations are built from all pairs and
        then labelled by their database keys, e.g.
        [map1_key]_A_with_B_x_[map2_key]_B_with_A.ini
        the noise_inv keys are not annotated in the filename but on can read
        off the files in the .ini file
        """
        (pairlist, pairdict) = \
                dp.cross_maps(map1_key, map2_key,
                              noise1_key, noise2_key,
                              map_suffix=";map",
                              noise_inv_suffix=";noise_inv",
                              cross_sym="_x_",
                              pair_former="GBTauto_cross_pairs",
                              ignore=['param'],
                              tag1prefix=map1_key + "_",
                              tag2prefix=map2_key + "_",
                              verbose=False)

        for pairname in pairlist:
            ininame = pairname + ".ini"
            inifile = open(self.ini_root + ininame, 'w')
            print ininame
            fileinfo = pairdict[pairname]
            inifile.write("import os\nimport scipy as sp\n\n")
            inifile.write("fs_output_root = '%s'\n" % \
                          self.params['output_root'])
            inifile.write("fs_output_filetag = '%s'\n" % pairname)
            inifile.write("fs_map1 = '%s'\n" % fileinfo['map1'])
            inifile.write("fs_noise_inv1 = '%s'\n" % fileinfo['noise_inv1'])
            inifile.write("fs_map2 = '%s'\n" % fileinfo['map2'])
            inifile.write("fs_noise_inv2 = '%s'\n" % fileinfo['noise_inv2'])
            inifile.write("fs_lags = %s\n" % self.lags)
            inifile.write("fs_freq_list = range(%d)\n" % self.nfreq_bin)
            cutlist = []
            for fbin in range(self.nfreq_bin):
                if fbin not in self.freq_list:
                    cutlist.append(fbin)
            inifile.write("cutlist = %s\n" % repr(cutlist))
            comp = '[x for x in fs_freq_list if x not in cutlist]'
            inifile.write("fs_freq_list = tuple(%s)\n" % comp)
            inifile.close()

if __name__ == '__main__':
    if len(sys.argv) == 2:
        BulkCorrelate(str(sys.argv[1])).execute()
    else:
        print 'Need one argument: parameter file name.'
