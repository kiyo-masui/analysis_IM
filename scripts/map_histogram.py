import h5py
import numpy as np
import matplotlib.pyplot as plt
from core import algebra
from utils import binning
import bin_map

def make_hist(map, bins, name, params = ['X-axis','Y-axis','Title'],
              set_logx = False, set_logy = False, savefig = True,
              plotdir = "/cita/d/www/home/mufma/plots/"):
    '''
    Make histogram of the map
    '''
    num_bin = bins
    hist, bin_edges, patches = plt.hist(map, num_bin, histtype='step')
    plt.xlabel(params[0])
    plt.ylabel(params[1])
    plt.title(params[2])
    if set_logx == True:
        plt.xscale('log')
    if set_logy == True:
        plt.yscale('log')
    if savefig == True:
        plt.savefig(plotdir + name, format = "png")

def open_data_hdf5(catalog_name, group, subgroup):
    '''
    Returns array of data from the catalog[group][subgroup]
    '''
    catalog = h5py.File(catalog_name)
    data = catalog[group][subgroup][:]
    return data



binned_map = bin_map.data_binning('/cita/h/home-2/mufma/code/
                                  analysis_IM/Tryout.hdf5', 41.,
                                  'HI_masses', 'HI_masses')

map = '/cita/h/home-2/mufma/code/analysis_IM/Tryout.hdf5'
group = '21_brightness_probe'
subgroup = '21_probe'
params = ['mK', 'number of clusters', '21 brightness fluctuations']
