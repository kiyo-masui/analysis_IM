import h5py
from core import algebra
from utils import binning
import numpy as np


def data_binning(h5py_file,num_bins):
    file = h5py.File(h5py_file)
    # reading data from h5py
    posx = file['Positions']['x'][:]
    posy = file['Positions']['y'][:]
    posz = file['Positions']['z'][:]
    halo_mass = file['Halo_Masses']['Halo_Masses'][:]
    # creating info for binning function
    xedges = np.linspace(0,posx.max(),num_bins)
    yedges = np.linspace(0,posy.max(),num_bins)
    zedges = np.linspace(0,posz.max(),num_bins)
    sample = np.zeros((len(posx),4))
    sample[:,0] = posx
    sample[:,1] = posy
    sample[:,2] = posz
    sample[:,3] = halo_mass
    # calling binning function
    z = binning.histogram3d_weights(sample,xedges,yedges,zedges)
    # info for make_cube_movie.py
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
    # save data for make_cube_movie.py
    map = algebra.make_vect(z, axis_names=('ra','dec','freq'))
    map.info = info
    save_file = open("/tmp/mufma/data/halos_20_bins_weights.npy","w")
    algebra.save(save_file,map)
    save_file.close()

data_binning('scripts/Halo_Catalog.hdf5',20.)

