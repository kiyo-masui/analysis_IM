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
import random
import sys
import numpy.ma as ma
from plotting import plot_slice
from numpy import linalg as LA
from utils import data_paths as dp
from kiyopy import parse_ini
# TODO: implement freq cut list: but this should already zero out the input weights
# TODO: make 1D version of this (shells vs. disks)
# TODO: is counts normalization OK; identity matrix recovered
# TODO: confirm that P(k=0) should be less than sum(w1*s2) for any mixing
# TODO: mixing of uniform weighting = identity? embedded in infinite region?
# TODO: analytic: mixing of uniform cube, asymptotically large?
# TODO: movie of bin in 3D k-space, bin convolved by the window
# TODO: is it ever correct to make this unitless?

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
                     window='blackman', zero_pad=False, identity_test=False):
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

    # no window applied here (applied above)
    xspec = pe.cross_power_est(weight1, weight2, arr1, arr2,
                               window=None, nonorm=True)

    # for each point in the cube, find |k|, k_perp, k_parallel
    # TODO: speed this up by using one direct numpy call (not limiting)
    k_mag_arr = binning.radius_array(xspec)
    k_perp_arr = binning.radius_array(xspec, zero_axes=[0])
    k_parallel_arr = binning.radius_array(xspec, zero_axes=[1, 2])

    if unitless:
        xspec = pe.make_unitless(xspec, radius_arr=k_mag_arr)

    # NOTE: assuming lowest k bin has only one point in 3D k-space
    # could make this floor of dimensions divided by 2 also
    center_3d = np.transpose(np.transpose(np.where(k_mag_arr == 0.))[0])

    # In the estimator, we devide by 1/sum(w1 * w2) to get most of the effect
    # of the weighing. The mixing matrix here can be thought of as a correction
    # that that diagonal-only estimate.
    leakage_ratio = xspec[center_3d[0], center_3d[1], center_3d[2]] / \
                    np.sum(weight1 * weight2)
    print "power leakage ratio: %10.5g" % leakage_ratio

    xspec /= np.sum(weight1 * weight2)

    print "partitioning the 3D kspace up into the 2D k bins"
    (kflat, ret_indices) = bin_indices_2d(k_perp_arr, k_parallel_arr,
                                          bins, bins)

    # perform a test where the window function is a delta function at the
    # origin so that the mixing matrix is identity
    if identity_test:
        xspec = algebra.zeros_like(xspec)
        xspec[center_3d[0], center_3d[1], center_3d[2]] = 1.

    # now save the window cross-power for downstream pooled users
    algebra.save(xspec_fileout, xspec)

    runlist = []
    for bin_index in range(kflat.shape[0]):
        bin_3d = ret_indices[repr(bin_index)]
        if bin_3d is not None:
            runlist.append((xspec_fileout, bin_index, bins, bin_3d, center_3d))

    pool = multiprocessing.Pool(processes=(multiprocessing.cpu_count() - 4))
    # the longest runs get pushed to the end; randomize for better job packing
    random.shuffle(runlist)
    results = pool.map(sum_window, runlist)
    #gnuplot_single_slice(runlist[0])  # for troubleshooting

    # now save the results for post-processing
    # TODO: move this into this function -- one stop shopping
    outshelve = shelve.open(mixing_fileout, "n")
    outshelve["results"] = results
    outshelve["xspec"] = xspec
    outshelve["bins"] = bins
    outshelve["weight_pair"] = [weight_file1, weight_file2]
    outshelve.close()


def load_mixing_matrix(filename):
    r"""load a mixing matrix produced by 'calculate_mixing'
    The goal is to unmixing the 2D k-bins, but many of the k-bins have zero
    counts, e.g. no bins in 3d k-space for that 2d bin. The inversion is
    therefore performed on a subset of the indices.
    """
    mixing_shelve = shelve.open(filename, "r")

    # first find the number of 2D k bins that have more than 0 modes
    counts_reference = mixing_shelve['results'][0][1].flatten()
    counts_nonzero = (counts_reference != 0)
    n_nonzero_counts = np.sum(counts_nonzero)
    where_nonzero_counts = np.where(counts_nonzero)[0]

    # alternately: (these should be the same)
    bins3d_nonzero = []
    for results_entry in mixing_shelve['results']:
        kindex = results_entry[0]
        bins3d_nonzero.append(kindex)

    # these are produced in a batch which is not necessarily ordered
    bins3d_nonzero.sort()
    bins3d_trans = {}
    nonzero_index = 0
    for full_index in bins3d_nonzero:
        bins3d_trans[repr(full_index)] = nonzero_index
        nonzero_index += 1

    n_nonzero_bins3d = len(bins3d_nonzero)

    print bins3d_nonzero
    print where_nonzero_counts
    assert np.array_equal(bins3d_nonzero, where_nonzero_counts), \
           "load_mixing_matrix: unmatching bins"

    mixing_array = np.zeros((n_nonzero_counts, n_nonzero_counts))
    counts_array = np.zeros((n_nonzero_counts, n_nonzero_counts))

    for result_entry in mixing_shelve['results']:
        kindex = result_entry[0]
        if kindex in bins3d_nonzero:
            counts_vec = result_entry[1].flatten()
            mixing_vec = result_entry[2].flatten()
            counts_vec_nonzero = counts_vec[bins3d_nonzero]
            mixing_vec_nonzero = mixing_vec[bins3d_nonzero]

            nonzero_index = bins3d_trans[repr(kindex)]
            counts_array[nonzero_index, :] = counts_vec_nonzero
            mixing_array[nonzero_index, :] = mixing_vec_nonzero

    mixing_shelve.close()

    #masked_mixing = ma.masked_invalid(mixing_array)
    zeroed_mixing = copy.deepcopy(mixing_array)
    zeroed_mixing[np.isnan(zeroed_mixing)] = 0.
    zeroed_mixing[np.isinf(zeroed_mixing)] = 0.
    #print zeroed_mixing

    inv_mixing = LA.pinv(zeroed_mixing)

    plot_slice.simpleplot_2D("mixing.png", zeroed_mixing,
                         range(n_nonzero_bins3d), range(n_nonzero_bins3d),
                         ["X", "Y"], 1., "2D mixing", "mixing", logscale=False)

    print np.diag(inv_mixing)
    inv_mixing[inv_mixing < 0.] = -np.sqrt(-inv_mixing[inv_mixing < 0.])
    inv_mixing[inv_mixing > 0.] = np.sqrt(inv_mixing[inv_mixing > 0.])
    plot_slice.simpleplot_2D("inv_mixing.png", inv_mixing,
                         range(n_nonzero_bins3d), range(n_nonzero_bins3d),
                         ["X", "Y"], 1., "2D mixing", "mixing", logscale=False)

    plot_slice.simpleplot_2D("counts.png", counts_array,
                         range(n_nonzero_bins3d), range(n_nonzero_bins3d),
                         ["X", "Y"], 1., "2D counts", "counts", logscale=False)


calc_mixing_init = {
        "map_key": "",
        "combined_map_key": "test",
        "wigglez_sel_key": "test",
        "perpair_base": "test_file",
        "outfile": "mixing_matrix.shelve",
        "xspec_file": "./xspec.npy", 
        "refinement": 2,
        "pad": 5,
        "order": 1,
        "window": "blackman",
        "summary_only": False,
        "bins": [0.00765314, 2.49977141, 35]
               }
calc_mixing_prefix = 'cm_'

class CalcMixingMatrix(object):
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        self.datapath_db = dp.DataPath()

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          calc_mixing_init,
                                          prefix=calc_mixing_prefix)

        bin_spec = self.params["bins"]
        self.bins = np.logspace(math.log10(bin_spec[0]),
                           math.log10(bin_spec[1]),
                           num=bin_spec[2], endpoint=True)

    def execute(self, processes):
        self.generate_runlist()
        if not self.params['summary_only']:
            if self.params['wigglez_sel_key']:
                self.execute_wigglez_calc()

            self.execute_calc()

        self.execute_summary()

    def generate_runlist(self):
        map_key = self.params['map_key']
        map_cases = self.datapath_db.fileset_cases(map_key,
                                                   "pair;type;treatment")

        # first compute the mixing matrix for the crossed (AxB) weightings
        unique_pairs = dp.GBTauto_cross_pairs(map_cases['pair'],
                                              map_cases['pair'],
                                              cross_sym="_with_")

        # assume the weights are the same for all cleaning treatments
        # TODO: this may change in the future
        treatment = "0modes"
        self.all_pairs = {}
        for item in unique_pairs:
            dbkeydict = {}
            mapset0 = (map_key, item[0], treatment)
            mapset1 = (map_key, item[1], treatment)
            dbkeydict['noiseinv1_key'] = "%s:%s;noise_inv;%s" % mapset0
            dbkeydict['noiseinv2_key'] = "%s:%s;noise_inv;%s" % mapset1
            files = dp.convert_dbkeydict_to_filedict(dbkeydict,
                                                     datapath_db=self.datapath_db)

            self.all_pairs[item[0]] = (files['noiseinv1_key'],
                                   files['noiseinv2_key'])

        # For the autopower (in noise assessment), we use the same cleaned maps
        # and the weights are the same for various pairs, e.g.
        # A_with_B is the same as A_with_C, etc. because the mode cleaning does
        # not impact the weighting functions
        A_file = (self.datapath_db.fetch(
                '%s:A_with_B;noise_inv;0modes' % map_key, silent=True))
        self.all_pairs["A_with_A"] = (A_file, A_file)

        B_file = (self.datapath_db.fetch(
                '%s:B_with_A;noise_inv;0modes' % map_key, silent=True))
        self.all_pairs["B_with_B"] = (B_file, B_file)

        C_file = (self.datapath_db.fetch(
                '%s:C_with_A;noise_inv;0modes' % map_key, silent=True))
        self.all_pairs["C_with_C"] = (C_file, C_file)

        D_file = (self.datapath_db.fetch(
                '%s:D_with_A;noise_inv;0modes' % map_key, silent=True))
        self.all_pairs["D_with_D"] = (D_file, D_file)

        self.mixing_fileout = {}
        for pair in self.all_pairs:
            self.mixing_fileout[pair] = "%s_%s.shelve" % \
                                   (self.params['perpair_base'], pair)


    def execute_calc(self):
        for pair in self.all_pairs:
            weight_file1, weight_file2 = self.all_pairs[pair]
            print pair, weight_file1, weight_file2, self.mixing_fileout[pair]

            calculate_mixing(weight_file1, weight_file2, self.bins,
                             self.params['xspec_file'],
                             self.mixing_fileout[pair],
                             unitless=False,
                             refinement=self.params['refinement'],
                             pad=self.params['pad'],
                             order=self.params['order'],
                             window=self.params['window'],
                             zero_pad=False, identity_test=False)

    def execute_wigglez_calc(self):
        r"""TODO: finish this once we need mixing matrices for xspec"""
        map_files = self.datapath_db.fetch(self.params['combined_map_key'])
        combined_weightfile = map_files["weight;0modes"]
        wigglez_files = self.datapath_db.fetch(self.params['wigglez_sel_key'])

    def execute_summary(self):
        #load_mixing_matrix(self.cross_pairs)
        print "ok"

#bins = np.logspace(math.log10(0.00765314),
#                           math.log10(2.49977141),
#                           num=35, endpoint=True)

#weightdir = "/mnt/raid-project/gmrt/eswitzer/GBT/maps/june15/"
#weight_file1 = weightdir + "secA_15hr_41-90_noise_weight_I_762.npy"
#weight_file2 = weightdir + "secB_15hr_41-90_noise_weight_I_762.npy"
#xspec_fileout = "./xspec_identity.npy"
#mixing_fileout = "mixing_summary_identity.shelve"
#calculate_mixing(weight_file1, weight_file2, bins, xspec_fileout,
#                 mixing_fileout, identity_test=False, refinement=1)

#load_mixing_matrix("mixing_summary_identity.shelve")

