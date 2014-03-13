import numpy as np
import numpy.ma as ma
from quadratic_products import pwrspec_estimator as pe
from quadratic_products import power_spectrum as ps
import glob
import copy
from plotting import plot_slice
import shelve
import h5py
from kiyopy import parse_ini
from utils import file_tools
import shutil

# if do_transfer:
#     pipe_modules.append((aggregate_bulksim.SubtractMap, ('sub1_', 'sub_')))
# sub1_plussim_file = ast3_statfile_out
# sub1_mappower_file = sigpwr_base + basemap + ".shelve"
# sub1_output_file = "%s/%s_stat_modes.hd5" % (pwrout_root, output_tag)
subtractmap_init = {
        "plussim_file": "file",
        "mappower_file": "file",
        "output_file": "file"
    }
subtractmap_prefix = 'sub_'


class SubtractMap(object):
    """Input: clean(map+sim)xclean(map+sim), clean(map)xclean(map)
    subtract these for the purpose of finding the transfer funcion
    """

    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        np.seterr(invalid='raise')

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          subtractmap_init,
                                          prefix=subtractmap_prefix)

        print self.params["plussim_file"], "-", self.params["mappower_file"]
        print "writing to ", self.params["output_file"]

    def execute(self, processes):
        shutil.copyfile(self.params["plussim_file"], self.params["output_file"])
        outfile = h5py.File(self.params["output_file"], "r+")
        pwr_map = ps.PowerSpectrum(self.params["mappower_file"])

        signal2d_agg = pwr_map.agg_stat_2d_pwrspec()
        print signal2d_agg.keys()

        for treatment in outfile:
            plussim_power = \
                copy.deepcopy(outfile[treatment]['pk_2d_stat']['mean'].value)

            map_power = signal2d_agg[treatment]["mean"]
            del outfile[treatment]['pk_2d_stat']['mean']

            outfile[treatment]['pk_2d_stat']['mean'] = plussim_power - \
                                                       map_power

        outfile.close()



calculatedatalike_init = {
        "powerfile": "file",
        "powerdatalike_out": "file",
    }
calculatedatalike_prefix = 'cdl_'


class CalculateDatalike(object):
    """Reformat aggregated simulations to look like data (vestigial?)
    """

    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        np.seterr(invalid='raise')

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          calculatedatalike_init,
                                          prefix=calculatedatalike_prefix)

        print self.params["powerfile_in"], "->", self.params["powerdatalike_out"]

        self.stats_in = h5py.File(self.params["powerfile"], "r")
        self.treatments_in = self.stats_in["results"].keys()

        # maybe make this hd5?
        self.stats_dataout = shelve.open(self.params["powerdatalike_out"], "n")
        #file_tools.convert_numpytree_hdf5(out_dicttree, outfile)

    def execute(self, processes):
        """produce a list of files to combine and run"""
        for treatment in self.treatments_in:
            print treatment
            print self.stats_in["results"][treatment].keys()

            stat_in = self.stats_in["stats"]

            stats_1d_in = stat_in["0modes"]["pk_1d_stat"]["mean"]
            stats_2d_in = stat_in["0modes"]["pk_2d_stat"]["mean"]
            counts_1d_in = stat_in["0modes"]["pk_1d_counts"]["mean"]
            counts_2d_in = stat_in["0modes"]["pk_2d_counts"]["mean"]

            k_1d = self.stats_in["k_1d"]
            k_1d_from_2d = self.stats_in["k_1d_from_2d"]
            kx_2d = self.stats_in["kx_2d"]
            ky_2d = self.stats_in["ky_2d"]

            # make a package that looks like the data
            pwrspec2d_product = {}
            pwrspec2d_product["bin_x_left"] = kx_2d["left"]
            pwrspec2d_product["bin_x_center"] = kx_2d["center"]
            pwrspec2d_product["bin_x_right"] = kx_2d["right"]
            pwrspec2d_product["bin_y_left"] = ky_2d["left"]
            pwrspec2d_product["bin_y_center"] = ky_2d["center"]
            pwrspec2d_product["bin_y_right"] = ky_2d["right"]
            pwrspec2d_product["counts_histo"] = counts_2d_in
            pwrspec2d_product["binavg"] = stats_2d_in

            pwrspec1d_product = {}
            pwrspec1d_product["bin_left"] = k_1d["left"]
            pwrspec1d_product["bin_center"] = k_1d["center"]
            pwrspec1d_product["bin_right"] = k_1d["right"]
            pwrspec1d_product["counts_histo"] = counts_1d_in
            pwrspec1d_product["binavg"] = stats_1d_in

            datakey = "data:%s" % treatment
            self.stats_dataout[datakey] = (0, (pwrspec2d_product,
                                               pwrspec1d_product))

        self.stats_in.close()
        self.stats_dataout.close()


