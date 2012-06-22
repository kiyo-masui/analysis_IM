from utils import fftutil
import numpy as np
import scipy as sp
from map import physical_gridding as pg
from utils import batch_handler as bh
from quadratic_products import pwrspec_estimator as pe
from utils import binning
from core import algebra
import math
import multiprocessing
import shelve
import copy
# 1: movie of bin in k-space
# 2: recovery of diagonal mixing
# 3: movie of bin convolved by the window
# 4: summary of mixing matrix
# TODO is it ever correct to make this unitless?
# TODO: make 1D and 2D versions of this (shells vs. disks)
    # for each k roll the weight and sum, mult by  1/N_k
    # test this in the case that the window = 1 at the origin only
    # this assumes periodicity in the region!
    # can analytically find mixing from uniform cube to full space
    # sum(w1*s2) is the assumption that W = delta(k) P(w1*w2, k=0)
    # dignostic: plot of 3D w(k), k-discs before/after mixing

# given a radius array in 2D and bins, return 
# flatten the perp, parallel (or mag), and digitize
# dictionary out: for each bin (key) associate numpy vect of 3D indices
# roll with 0 padding

def sum_window(argt):
    """A given bin in 2D k-space (labelled by bin_index_2d) is the sum over a
    "washer" in 3D k-space, a band in k_parallel and an annulus in k_x, k_y.
    Let all of the 3D bins in k_space be indexed by bin_3d. The window function
    is centered at k_3d=0, and let these indices defining the center of the 3d
    volume be given in center_3d.
    TODO: replace this with NlogN convolution
    TODO: implement 0-padded roll instead of np.roll, algebra.roll_zeropad()
    """
    (filename, bin_index_2d, k_2d, bin_3d, center_3d) = argt
    # load the cross-power of the weighting functions
    xspec = algebra.make_vect(algebra.load(filename))
    windowsum = algebra.zeros_like(xspec)

    num_3dbins_in_2dbin = bin_3d.shape[0]

    print "%d: summing over %d bins" % (bin_index_2d, num_3dbins_in_2dbin)

    for bin_3dind in range(num_3dbins_in_2dbin):
        # TODO: is this sign right, does it matter?
        off = bin_3d[bin_3dind] - center_3d
        #print off
        windowsum += np.roll(np.roll(np.roll(xspec,
                                             off[0], axis=0),
                                             off[1], axis=1),
                                             off[2], axis=2)

    k_perp_arr = binning.radius_array(xspec, zero_axes=[0])
    k_parallel_arr = binning.radius_array(xspec, zero_axes=[1, 2])
    kx_2d = copy.deepcopy(k_2d)
    ky_2d = copy.deepcopy(k_2d)
    counts_histo_2d, binavg_2d = binning.bin_an_array_2d(windowsum,
                                                         k_perp_arr,
                                                         k_parallel_arr,
                                                         kx_2d, ky_2d)

    return (bin_index_2d, counts_histo_2d, binavg_2d)


def bin_indices_2d(k_perp_arr, k_parallel_arr, k_perp_bins, k_parallel_bins,
                   debug=False):
    r"""This partitions up the 3D P(k) space into disks of k in a given k_perp
    and k_parallel bin.
    k_perp_arr is the 3D array which gives the value of k_perp at each point
    k_parallel_arr is the 3D array which gives the value of k_parallel
    k_perp_bins is the vector of k_perp bins
    k_parallel_bins is the vector of k_parallel bins
    """
    # bin0 defines the left edge, and the final value determines the rightmost
    # edge, so the number of actual kbins is n-1
    # e.g. if kbins = 1<->2<->3, then there are two bins: [1:2), and [2,3]
    n_k_perp = len(k_perp_bins) - 1
    n_k_parallel = len(k_parallel_bins) - 1
    arr_shp = k_perp_arr.shape

    # assign an index to each 2d k bin and a Cartesian product giving the set
    # of possible indices to bins in 2d k-space
    kflat = algebra.cartesian((range(n_k_perp), range(n_k_parallel)))

    # find which points in 3d k-space (flattened onto one index) fall in the
    # various k perp and parallel bins
    k_perp_indices = np.digitize(k_perp_arr.flatten(), k_perp_bins)
    k_parallel_indices = np.digitize(k_parallel_arr.flatten(), k_parallel_bins)

    # note that np.digitize([0.5,1.5,2.5,3.5,4.5], [1,2,3]) = [0 1 2 3 3]
    # but we want e.g. 1.5 to go in the 0th bin
    k_perp_indices -= 1
    k_parallel_indices -= 1

    ret_indices = {}
    # loop through each index to the 2D k perp and parallel vectors
    # and return the list of indices to the points in 3D k-space that fall into
    # that bin.
    for i_perp in range(n_k_perp):
        print "doing perp index %d" % i_perp
        #i_perp_shell = np.where(k_perp_indices == i_perp)[0]
        i_perp_shell = k_perp_indices == i_perp
        for i_parallel in range(n_k_parallel):
            #i_parallel_shell = np.where(k_parallel_indices == i_parallel)[0]
            #index_in_shell = [ind for ind in i_perp_shell \
            #                       if ind in i_parallel_shell]
            i_parallel_shell = k_parallel_indices == i_parallel
            index_in_shell = np.where(np.logical_and(i_parallel_shell,
                                                     i_perp_shell))[0]

            # find the index in the flattened 2D k-space for this bin
            flatk_ind = np.where(np.logical_and(kflat[:,0] == i_perp,
                                        kflat[:,1] == i_parallel))[0][0]
            # writing to a dictionary, so make it a string key
            flatk_ind = repr(flatk_ind)
            # numpy > 1.6
            #flatk_ind = np.ravel_multi_index([[i_perp], [i_parallel]],
            #                                 mode='raise')[0]

            # for each flattened index in 3d k-space, find its unraveled 3d
            # index; if there are no 3d k-space bins in this 2d bin, None
            num_in_shell = len(index_in_shell)
            if num_in_shell > 0:
                # numpy > 1.6
                #ret_indices[flatk_ind] = np.unravel_index(
                #                                index_in_shell, arr_shp)
                ind_in_bin = np.zeros((num_in_shell, k_perp_arr.ndim),
                                      dtype=int)

                for (count, ind3d) in zip(range(num_in_shell), index_in_shell):
                    ind_in_bin[count, :] = \
                                np.unravel_index(index_in_shell[count], arr_shp)

                ret_indices[flatk_ind] = ind_in_bin
            else:
                ret_indices[flatk_ind] = None

            if debug:
                print arr_shp, i_perp, i_parallel, index_in_shell
                print ret_indices[flatk_ind]

    return (kflat, ret_indices)


# TODO memoize this
def calculate_mixing(weight_file1, weight_file2, bins, xspec_fileout,
                     mixing_fileout,
                     unitless=False, refinement=2, pad=5, order=1,
                     window='blackman', zero_pad=False, unity_test=False):
    print "loading the weights and converting to physical coordinates"
    weight1_obs = algebra.make_vect(algebra.load(weight_file1))
    weight1 = bh.repackage_kiyo(pg.physical_grid(
                                weight1_obs,
                                refinement=refinement,
                                pad=pad, order=order))

    weight2_obs = algebra.make_vect(algebra.load(weight_file2))
    weight2 = bh.repackage_kiyo(pg.physical_grid(
                                weight2_obs,
                                refinement=refinement,
                                pad=pad, order=order))

    if window:
        window_function = fftutil.window_nd(weight1.shape, name=window)
        weight1 *= window_function
        weight2 *= window_function

    print "calculating the cross-power of the spatial weighting functions"
    arr1 = algebra.ones_like(weight1)
    arr2 = algebra.ones_like(weight2)
    xspec = pe.cross_power_est(weight1, weight2, arr1, arr2)

    # for each point in the cube, find |k|, k_perp, k_parallel
    # TODO: speed this up by using one direct numpy call (not limiting)
    k_mag_arr = binning.radius_array(xspec)
    k_perp_arr = binning.radius_array(xspec, zero_axes=[0])
    k_parallel_arr = binning.radius_array(xspec, zero_axes=[1, 2])

    if unitless:
        xspec = pe.make_unitless(xspec, radius_arr=k_mag_arr)

    print "partitioning the 3D kspace up into the 2D k bins"
    (kflat, ret_indices) = bin_indices_2d(k_perp_arr, k_parallel_arr,
                                          bins, bins)

    # NOTE: assuming lowest k bin has only one point in 3D k-space
    # could make this floor of dimensions divided by 2 also
    zerobin = ret_indices[repr(0)]
    if zerobin.shape == (1,3):
        center_3d = zerobin[0,:]
    else:
        print "you have no bin at exactly = 0"
        return

    # perform a test where the window function is a delta function at the
    # origin so that the mixing matrix is unity
    if unity_test:
        xspec = algebra.zeros_like(xspec)
        xspec[center_3d] = 1.

    # now save the window cross-power for downstream pooled users
    algebra.save(xspec_fileout, xspec)

    runlist = []
    for bin_index in range(kflat.shape[0]):
        bin_3d = ret_indices[repr(bin_index)]
        if bin_3d is not None:
            runlist.append((xspec_fileout, bin_index, bins, bin_3d, center_3d))

    pool = multiprocessing.Pool(processes=(multiprocessing.cpu_count() - 4))
    results = pool.map(sum_window, runlist)
    #gnuplot_single_slice(runlist[0])  # for troubleshooting

    # now save the results for post-processing
    # TODO: move this into this function -- one stop shopping
    outshelve = shelve.open(mixing_fileout, "n")
    outshelve["results"] = results
    outshelve["bins"] = bins
    outshelve.close()


bin_spec = [0.00765314, 2.49977141, 35]
bins = np.logspace(math.log10(bin_spec[0]),
                       math.log10(bin_spec[1]),
                       num=bin_spec[2], endpoint=True)

weightdir = "/mnt/raid-project/gmrt/eswitzer/GBT/maps/june15/"
weight_file1 = weightdir + "secA_15hr_41-90_noise_weight_I_762.npy"
weight_file2 = weightdir + "secB_15hr_41-90_noise_weight_I_762.npy"
xspec_fileout = "./xspec_unity.npy"
mixing_fileout = "mixing_summary_unity.shelve"
calculate_mixing(weight_file1, weight_file2, bins, xspec_fileout,
                mixing_fileout, unity_test=True)
