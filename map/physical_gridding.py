import numpy as np
import scipy as sp
import scipy.ndimage
from core import algebra
from utils import data_paths
from utils import units
from utils import cosmology as cosmo
from utils import batch_handler


#@batch_handler.memoize_persistent
def physical_grid(input_array, refinement=2, pad=5, order=2):
    r"""Project from freq, ra, dec into physical coordinates

    Parameters
    ----------
    input_array: np.ndarray
        The freq, ra, dec map

    Returns
    -------
    cube: np.ndarray
        The cube projected back into physical coordinates

    """
    freq_axis = input_array.get_axis('freq') / 1.e6
    ra_axis = input_array.get_axis('ra')
    dec_axis = input_array.get_axis('dec')

    nu_lower, nu_upper = freq_axis.min(), freq_axis.max()
    ra_fact = sp.cos(sp.pi * input_array.info['dec_centre'] / 180.0)
    thetax, thetay = np.ptp(ra_axis), np.ptp(dec_axis)
    thetax *= ra_fact
    (numz, numx, numy) = input_array.shape

    cosmology = cosmo.Cosmology()
    z1 = units.nu21 / nu_upper - 1.0
    z2 = units.nu21 / nu_lower - 1.0
    d1 = cosmology.proper_distance(z1)
    d2 = cosmology.proper_distance(z2)
    c1 = cosmology.comoving_distance(z1)
    c2 = cosmology.comoving_distance(z2)
    c_center = (c1 + c2) / 2.

    # Make cube pixelisation finer, such that angular cube will
    # have sufficient resolution on the closest face.
    phys_dim = np.array([c2 - c1,
                         thetax * d2 * units.degree,
                         thetay * d2 * units.degree])

    # Note that the ratio of deltas in Ra, Dec in degrees may
    # be different than the Ra, Dec in physical coordinates due to
    # rounding onto this grid
    n = np.array([numz, int(d2 / d1 * numx), int(d2 / d1 * numy)])

    # Enlarge cube size by `pad` in each dimension, so raytraced cube
    # sits exactly within the gridded points.
    phys_dim = phys_dim * (n + pad).astype(float) / n.astype(float)
    c1 = c_center - (c_center - c1) * (n[0] + pad) / float(n[0])
    c2 = c_center + (c2 - c_center) * (n[0] + pad) / float(n[0])
    n = n + pad
    # now multiply by scaling for a finer sub-grid
    n = refinement * n

    print "converting from obs. to physical coord refinement=%s, pad=%s" % \
                      (refinement, pad)

    print "(%d, %d, %d)->(%f to %f) x %f x %f (%d, %d, %d) (h^-1 cMpc)^3" % \
                      (numz, numx, numy, c1, c2, \
                       phys_dim[1], phys_dim[2], \
                       n[0], n[1], n[2])

    # this is wasteful in memory, but numpy can be pickled
    phys_map_npy = np.zeros(n)
    phys_map = algebra.make_vect(phys_map_npy, axis_names=('freq', 'ra', 'dec'))
    #mask = np.ones_like(phys_map)
    mask = np.ones_like(phys_map_npy)

    # TODO: should this be more sophisticated? N-1 or N?
    info = {}
    info['axes'] = ('freq', 'ra', 'dec')
    info['type'] = 'vect'

    #info = {'freq_delta': abs(phys_dim[0])/float(n[0]),
    #        'freq_centre': abs(c2+c1)/2.,
    info['freq_delta'] = abs(c2 - c1) / float(n[0] - 1)
    info['freq_centre'] = c1 + info['freq_delta'] * float(n[0] // 2)

    info['ra_delta'] = abs(phys_dim[1]) / float(n[1] - 1)
    #info['ra_centre'] = info['ra_delta'] * float(n[1] // 2)
    info['ra_centre'] = 0.

    info['dec_delta'] = abs(phys_dim[2]) / float(n[2] - 1)
    #info['dec_centre'] = info['dec_delta'] * float(n[2] // 2)
    info['dec_centre'] = 0.

    phys_map.info = info
    print info

    # same as np.linspace(c1, c2, n[0], endpoint=True)
    radius_axis = phys_map.get_axis("freq")
    x_axis = phys_map.get_axis("ra")
    y_axis = phys_map.get_axis("dec")

    # Construct an array of the redshifts on each slice of the cube.
    comoving_inv = cosmo.inverse_approx(cosmology.comoving_distance, z1, z2)
    za = comoving_inv(radius_axis)  # redshifts on the constant-D spacing
    nua = units.nu21 / (1. + za)

    gridy, gridx = np.meshgrid(y_axis, x_axis)
    interpol_grid = np.zeros((3, n[1], n[2]))

    for i in range(n[0]):
        # nua[0] = nu_upper, nua[1] = nu_lower
        #print nua[i], freq_axis[0], freq_axis[-1], (nua[i] - freq_axis[0]) / \
        #                                (freq_axis[-1] - freq_axis[0]) * numz

        interpol_grid[0, :, :] = (nua[i] - freq_axis[0]) / \
                               (freq_axis[-1] - freq_axis[0]) * numz
        proper_z = cosmology.proper_distance(za[i])

        angscale = proper_z * units.degree
        interpol_grid[1, :, :] = gridx / angscale / thetax * numx + numx / 2
        interpol_grid[2, :, :] = gridy / angscale / thetay * numy + numy / 2

        phys_map_npy[i, :, :] = sp.ndimage.map_coordinates(input_array,
                                                           interpol_grid,
                                                           order=order)

        interpol_grid[1, :, :] = np.logical_or(interpol_grid[1, :, :] > numx,
                                             interpol_grid[1, :, :] < 0)
        interpol_grid[2, :, :] = np.logical_or(interpol_grid[1, :, :] > numy,
                                             interpol_grid[1, :, :] < 0)
        mask = np.logical_not(np.logical_or(interpol_grid[1, :, :],
                                            interpol_grid[2, :, :]))
        phys_map_npy *= mask

    return phys_map_npy, info


if __name__ == '__main__':
    datapath_db = data_paths.DataPath()
    filename = datapath_db.fetch('simideal_15hr_beam', intend_read=True,
                                 pick='1')
    print filename
    cube = algebra.make_vect(algebra.load(filename))
    pcube = physical_grid(cube, refinement=2)
    algebra.save("physical_cube_beam.npy", pcube)

    filename = datapath_db.fetch('simideal_15hr', intend_read=True,
                                 pick='1')
    print filename
    cube = algebra.make_vect(algebra.load(filename))
    pcube = physical_grid(cube, refinement=2)
    algebra.save("physical_cube.npy", pcube)
