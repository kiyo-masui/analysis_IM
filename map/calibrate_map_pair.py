r"""
Cross-calibrate a large set of maps
"""
import numpy as np
from utils import data_paths as dp
from correlate import corr_estimation as ce
from core import algebra
from correlate import map_pair
import shelve
import copy


def noise_inv_to_weight(noiseinvlist_in, weightlist_in):
    print reference_noise_inv, reference_weight
    noise_inv = algebra.make_mat(algebra.open_memmap(reference_noise_inv, mode='r'))
    noise_inv_diag = noise_inv.mat_diag()
    algebra.save(reference_weight, noise_inv_diag)
    print "done with reference weights"

    for (noiseinv_item, weight_item) in zip(noiseinvlist_in, weightlist_in):
        print noiseinv_item, weight_item
        noise_inv = algebra.make_mat(algebra.open_memmap(noiseinv_item, mode='r'))
        noise_inv_diag = noise_inv.mat_diag()
        algebra.save(weight_item, noise_inv_diag)


def remove_mean(map_in, weight):
    meansmap = np.sum(np.sum(weight * map_in, -1), -1)
    meansmap /= np.sum(np.sum(weight, -1), -1)
    meansmap.shape += (1, 1)
    map_return = copy.deepcopy(map_in)
    map_return -= meansmap

    return map_return


# TODO: implement bias correction
def template_fit(session_map, template_map, weight,
                 bias_num=None, bias_denom=None):
    freq_axis = template_map.get_axis('freq')
    n_freq = len(freq_axis)

    cal_factor = np.zeros((n_freq))

    for index in range(len(freq_axis)):
        cal_factor[index] = np.sum(session_map[index, :, :] *
                                   weight[index, :, :] *
                                   template_map[index, :, :])

        cal_factor[index] /= np.sum(template_map[index, :, :] *
                                    weight[index, :, :] *
                                    template_map[index, :, :])

    cal_factor[np.isnan(cal_factor)] = 0.

    return cal_factor


def map_pair_cal(uncal_maplist, uncal_weightlist, calfactor_outlist,
                 dirtymap_inlist, dirtymap_outlist,
                 reference_mapfile, reference_weightfile,
                 sub_weighted_mean=True, freq_list=range(256)):

    reference_map = algebra.make_vect(algebra.load(reference_mapfile))
    reference_weight = algebra.make_vect(algebra.load(reference_weightfile))

    reference_map = remove_mean(reference_map, reference_weight)

    # load maps into pairs
    for mapfile, weightfile, calfactor_outfile, \
        dirty_infile, dirty_outfile in zip(uncal_maplist, \
            uncal_weightlist, calfactor_outlist,
            dirtymap_inlist, dirtymap_outlist):

        print mapfile, weightfile

        session_map = algebra.make_vect(algebra.load(mapfile))
        session_weight = algebra.make_vect(algebra.load(weightfile))

        session_map = remove_mean(session_map, session_weight)

        calfactor = template_fit(session_map, reference_map, session_weight)

        newmap = algebra.make_vect(algebra.load(dirty_infile))
        newmap[freq_list, :, :] /= calfactor[:, np.newaxis, np.newaxis]
        algebra.save(dirty_outfile, newmap)
        print dirty_outfile

        # optional test by applying the factor to the maps
        #session_map[freq_list, :, :] /= calfactor[:, np.newaxis, np.newaxis]
        #calfactor = template_fit(session_map, reference_map, session_weight)

        facout = open(calfactor_outfile, "w")
        for outvals in calfactor:
            facout.write("%10.15g\n" % outvals)

        facout.close()


if __name__ == '__main__':
    input_root = "/mnt/raid-project/gmrt/kiyo/gbt_out/maps/apr16.2012/"
    output_root = "/mnt/raid-project/gmrt/eswitzer/GBT/maps/self_calibrated_2/"

    live_sessions = range(41, 81+1)
    live_sessions.remove(54)
    live_sessions.remove(58)
    live_sessions.remove(60)
    live_sessions.remove(63)
    live_sessions.remove(67)
    live_sessions.remove(70)
    live_sessions.remove(72)
    live_sessions.remove(78)
    live_sessions.remove(79)

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

    reference_clean = input_root + "sess_most_calib_15hr_41-73_clean_map_I.npy"
    reference_dirty = input_root + "sess_most_calib_15hr_41-73_dirty_map_I.npy"
    reference_noise_diag = input_root + "sess_most_calib_15hr_41-73_noise_diag_I.npy"
    reference_noise_inv = input_root + "sess_most_calib_15hr_41-73_noise_inv_I.npy"
    reference_weight = output_root + "sess_most_calib_15hr_41-73_weight_I.npy"

    standard_freq_list = (0, 1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 14, 17, 23, 24, 25,
                      26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 39, 40,
                      41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
                      55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68,
                      69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 81, 82, 83,
                      84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97,
                      98, 99, 100, 101, 102, 109, 110, 111, 112, 113, 114, 115,
                      116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126,
                      127, 128, 129, 135, 136, 137, 138, 139, 140, 141, 142,
                      143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153,
                      154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164,
                      165, 166, 167, 168, 169, 170, 172, 173, 174, 176, 178,
                      180, 181, 184, 185, 186, 188, 190, 191, 199, 200, 202,
                      203, 205, 206, 207, 210, 211, 214, 215, 216, 217, 220,
                      221, 222, 223, 224, 225, 226, 227, 228, 230, 231, 232,
                      234, 235, 236, 238, 239, 240, 241, 242, 243, 245, 246,
                      247, 248, 249, 250, 251, 252, 253)

    #noise_inv_to_weight(noiseinvlist_fullpath, weightlist_fullpath)
    map_pair_cal(maplist_fullpath, weightlist_fullpath, calfactorlist_fullpath,
                 dirtymaplist_fullpath, recal_dirtymaplist_fullpath,
                 reference_clean, reference_weight, freq_list=range(256))
