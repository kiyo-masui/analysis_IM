import h5py
from core import algebra
from utils import binning
import numpy as np

def converter_halo(catalog, num_bins, group, subgroup, directory, savefile):
    """
    Convert h5py catalog to binned algebra map with making deltas.
    """
    file = h5py.File(catalog, 'r')
    # reading data from h5py
    posx = file['Positions']['x'][:]
    posy = file['Positions']['y'][:]
    posz = file['Positions']['z'][:]
    weights = file[group][subgroup][:]
    file.close()
    # creating info for binning function
    xedges = np.linspace(0, posx.max(), num_bins)
    yedges = np.linspace(0, posy.max(), num_bins)
    zedges = np.linspace(0, posz.max(), num_bins)
    sample = np.zeros((len(posx), 4))
    sample[:, 0] = posx
    sample[:, 1] = posy
    sample[:, 2] = posz
    sample[:, 3] = weights
    # calling binning function
    binned_map = binning.histogram3d_weights(sample, xedges, yedges, zedges)
    # prepare deltas
    binned_map = (binned_map - binned_map.mean())/(binned_map.mean())
    # info for pipeline
    x_delta = posx.max()/num_bins
    y_delta = posy.max()/num_bins
    z_delta = posz.max()/num_bins
    x_centre = posx.max()/2.
    y_centre = posy.max()/2.
    z_centre = posz.max()/2.
    info = {'ra_delta': x_delta,
            'dec_delta': y_delta,
            'dec_centre': y_centre,
            'axes': ('freq', 'ra', 'dec'),
            'ra_centre': x_centre,
            'freq_centre': z_centre,
            'freq_delta': z_delta,
            'type': 'vect'}
    # save data for pipeline in manager.py
    map = algebra.make_vect(binned_map, axis_names=('ra', 'dec', 'freq'))
    map.info = info
    save_file = open(directory%(num_bins-1) + savefile, "w")
    algebra.save(save_file, map)
    save_file.close()

def bin_map(catalog, num_bins, group, subgroup):
    """
    Read h5py catalog, return binned algebra map without making deltas.
    Output: map,info.
    """
    file = h5py.File(catalog, 'r')
    # read data from h5py catalogs
    posx = file['Positions']['x'][:]
    posy = file['Positions']['y'][:]
    posz = file['Positions']['z'][:]
    weights = file[group][subgroup][:]
    file.close()
    # create info for binning function
    xedges = np.linspace(0, posx.max(), num_bins)
    yedges = np.linspace(0, posy.max(), num_bins)
    zedges = np.linspace(0, posz.max(), num_bins)
    sample = np.zeros((len(posx), 4))
    sample[:, 0] = posx
    sample[:, 1] = posy
    sample[:, 2] = posz
    sample[:, 3] = weights
    # create info for pipeline
    x_delta = posx.max()/num_bins
    y_delta = posy.max()/num_bins
    z_delta = posz.max()/num_bins
    x_centre = posx.max()/2.
    y_centre = posy.max()/2.
    z_centre = posz.max()/2.
    info = {'ra_delta': x_delta,
            'dec_delta': y_delta,
            'dec_centre': y_centre,
            'axes': ('freq', 'ra', 'dec'),
            'ra_centre': x_centre,
            'freq_centre': z_centre,
            'freq_delta': z_delta,
            'type': 'vect'}
    return binning.histogram3d_weights(sample, xedges, yedges, zedges), info

def join_index_algebra(catalog_list, num_bins, group, subgroup,
                       directory, savefile):
    binned_map, info = bin_map(catalog_list[0], num_bins, group, subgroup)
    for num in range(1,len(catalog_list)):
        binned_map = binned_map + bin_map(catalog_list[num], num_bins, group,
subgroup)[0]
    # save data for pipeline in manager.py
    map = algebra.make_vect(binned_map, axis_names=('ra', 'dec', 'freq'))
    map.info = info
    save_file = open(directory%(num_bins-1) + savefile, "w")
    algebra.save(save_file, map)
    save_file.close()

def convert_npy(directory):
    file = np.load(directory)
    x_delta = 0.62890625
    y_delta = 0.62890625
    z_delta = 0.62890625
    x_centre = 161.0
    y_centre = 161.0
    z_centre = 161.0
    info = {'ra_delta': x_delta,
            'dec_delta': y_delta,
            'dec_centre': y_centre,
            'axes': ('freq', 'ra', 'dec'),
            'ra_centre': x_centre,
            'freq_centre': z_centre,
            'freq_delta': z_delta,
            'type': 'vect'}
    box = algebra.make_vect(file, axis_names=('ra', 'dec', 'freq'))
    box.info = info
    directory = "/mnt/raid-project/gmrt/mufma/"
    savefile = "kappa512.npy"
    save_file = open(directory + savefile, "w")
    algebra.save(save_file, box)
    save_file.close()



