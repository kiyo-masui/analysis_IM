r"""
Cross-calibrate a large set of maps
"""
import numpy as np
from utils import data_paths as dp
from correlate import corr_estimation as ce
from core import algebra
from correlate import map_pair
import shelve

input_root = "/mnt/raid-project/gmrt/nbanavar/gbt_out/maps/sessions/"
output_root = "/mnt/raid-project/gmrt/eswitzer/GBT/maps/self_calibrated/"

live_sessions = range(41, 74)
live_sessions.remove(54)
live_sessions.remove(58)
live_sessions.remove(60)
live_sessions.remove(63)
live_sessions.remove(67)
live_sessions.remove(70)
live_sessions.remove(72)


maplist = ["sess_%d_calib_15hr_41-73_clean_map_I.npy" % item for item in live_sessions]
maplist_fullpath = [ input_root + item for item in maplist]

dirtymaplist = ["sess_%d_calib_15hr_41-73_dirty_map_I.npy" % item for item in live_sessions]
dirtymaplist_fullpath = [ input_root + item for item in dirtymaplist]
recal_dirtymaplist_fullpath = [ output_root + item for item in dirtymaplist]

noiseinvlist = ["sess_%d_calib_15hr_41-73_noise_inv_I.npy" % item for item in live_sessions]
noiseinvlist_fullpath = [ input_root + item for item in noiseinvlist]

weightlist = ["sess_%d_calib_15hr_41-73_weight_I.npy" % item for item in live_sessions]
weightlist_fullpath = [ output_root + item for item in weightlist]

calfactorlist = ["sess_%d_calib_15hr_41-73_cal_factor_I.txt" % item for item in live_sessions]
calfactorlist_fullpath = [ output_root + item for item in calfactorlist]

standard_freq_list = (0, 1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 14, 17, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 172, 173, 174, 176, 178, 180, 181, 184, 185, 186, 188, 190, 191, 199, 200, 202, 203, 205, 206, 207, 210, 211, 214, 215, 216, 217, 220, 221, 222, 223, 224, 225, 226, 227, 228, 230, 231, 232, 234, 235, 236, 238, 239, 240, 241, 242, 243, 245, 246, 247, 248, 249, 250, 251, 252, 253)

def noise_inv_to_weight(noiseinvlist_in, weightlist_in):
    for (noiseinv_item, weight_item) in zip(noiseinvlist_in, weightlist_in):
        print noiseinv_item, weight_item
        noise_inv = algebra.make_mat(algebra.open_memmap(noiseinv_item, mode='r'))
        noise_inv_diag = noise_inv.mat_diag()
        algebra.save(weight_item, noise_inv_diag)


def map_pair_cal(uncal_maplist, uncal_weightlist, calfactor_outlist,
                 dirtymap_inlist, dirtymap_outlist,
                 convolve=True, factorizable_noise=True,
                 sub_weighted_mean=True, freq_list=range(256)):

    map1file = uncal_maplist.pop(0)
    weight1file = uncal_weightlist.pop(0)
    calfactor_outlist.pop(0)
    dirtymap_out0 = dirtymap_outlist.pop(0)
    dirtymap_in0 = dirtymap_inlist.pop(0)

    # do nothing to the reference map
    ref_dirtymap = algebra.make_vect(algebra.load(dirtymap_in0))
    algebra.save(dirtymap_out0, ref_dirtymap)

    # load maps into pairs
    svdout = shelve.open("correlation_pairs.shelve")
    for map2file, weight2file, calfactor_outfile, \
        dirty_infile, dirty_outfile in zip(uncal_maplist, \
            uncal_weightlist, calfactor_outlist,
            dirtymap_inlist, dirtymap_outlist):

        print map1file, weight1file, map2file, weight2file

        pair = map_pair.MapPair(map1file, map2file,
                                weight1file, weight2file,
                                freq_list, avoid_db=True)

        if factorizable_noise:
            pair.make_noise_factorizable()

        if sub_weighted_mean:
            pair.subtract_weighted_mean()

        if convolve:
            pair.degrade_resolution()

        (corr, counts) = pair.correlate()
        svd_info = ce.get_freq_svd_modes(corr, len(freq_list))
        svdout[map2file] = svd_info

        # write out the left right and cal factors
        leftmode = svd_info[1][0]
        rightmode = svd_info[2][0]
        calfactor = leftmode/rightmode

        facout = open(calfactor_outfile, "w")
        for outvals in zip(leftmode, rightmode, calfactor):
            facout.write("%10.15g %10.15g %10.15g\n" % outvals)

        facout.close()

        newmap = algebra.make_vect(algebra.load(dirty_infile))
        newmap[freq_list, :, :] *= calfactor[:,np.newaxis,np.newaxis]
        algebra.save(dirty_outfile, newmap)
        print dirty_outfile

    svdout.close()

if __name__ == '__main__':
    #noise_inv_to_weight(noiseinvlist_fullpath, weightlist_fullpath)
    map_pair_cal(maplist_fullpath, weightlist_fullpath, calfactorlist_fullpath,
                 dirtymaplist_fullpath, recal_dirtymaplist_fullpath,
                 freq_list=standard_freq_list)
    #map_pair_cal(maplist_fullpath, weightlist_fullpath)
