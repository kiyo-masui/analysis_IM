from kiyopy import parse_ini
import numpy as np
from utils import data_paths as dp
import os
import shutil
from core import algebra

subtractmap_init = {
        "map_key_1": "test_map",
        "tack_on_1": None,
        "map_key_2": "test_map",
        "tack_on_2": None,
        "map_key_out": "test_map",
        "tack_on_out": None,
               }
subtractmap_prefix = 'sm_'


class SubtractMaps(object):
    r"""given two database key for maps, subtract
    """
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        self.datapath_db = dp.DataPath()

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          subtractmap_init,
                                          prefix=subtractmap_prefix)

    def execute(self, processes):
        file_list_1 = self.datapath_db.fetch(self.params['map_key_1'],
                                             tack_on=self.params["tack_on_1"],
                                             silent=True)

        file_list_2 = self.datapath_db.fetch(self.params['map_key_2'],
                                             tack_on=self.params["tack_on_2"],
                                             silent=True)

        file_list_out = self.datapath_db.fetch(self.params['map_key_out'],
                                           tack_on=self.params["tack_on_out"],
                                           silent=True)

        for file_key in file_list_1[0]:
            infile = file_list_1[1][file_key]
            subfile = file_list_2[1][file_key]
            outfile = file_list_out[1][file_key]

            rootdir = "/".join(outfile.split("/")[0:-1])
            if len(rootdir) > 0 and rootdir != ".":
                if not os.path.isdir(rootdir):
                    print "print_multicolumn: making dir " + rootdir
                    os.mkdir(rootdir)

            print "input: ", infile
            print "out: ", outfile

            if "map" in file_key:
                print "minus: ", subfile
                inmap = algebra.make_vect(algebra.load(infile))
                submap = algebra.make_vect(algebra.load(subfile))
                print inmap.shape, submap.shape
                algebra.save(outfile, inmap - submap)
            else:
                shutil.copy2(infile, outfile)
                shutil.copy2(infile + ".meta", outfile + ".meta")
