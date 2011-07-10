"""make weighted average of data cubes"""
import numpy as np
from core import algebra
from correlate import wigglez_xcorr as wxc

fourway_split = {
    'root_data': "/mnt/raid-project/gmrt/calinliv/wiggleZ/corr/test/",
    'maplist': ["sec_A_15hr_41-73_cleaned_clean_map_I_with_B",
               "sec_A_15hr_41-73_cleaned_clean_map_I_with_C",
               "sec_A_15hr_41-73_cleaned_clean_map_I_with_D",
               "sec_B_15hr_41-73_cleaned_clean_map_I_with_A",
               "sec_B_15hr_41-73_cleaned_clean_map_I_with_C",
               "sec_B_15hr_41-73_cleaned_clean_map_I_with_D",
               "sec_C_15hr_41-73_cleaned_clean_map_I_with_A",
               "sec_C_15hr_41-73_cleaned_clean_map_I_with_B",
               "sec_C_15hr_41-73_cleaned_clean_map_I_with_D",
               "sec_D_15hr_41-73_cleaned_clean_map_I_with_A",
               "sec_D_15hr_41-73_cleaned_clean_map_I_with_B",
               "sec_D_15hr_41-73_cleaned_clean_map_I_with_C"],
    'covlist': ["sec_A_15hr_41-73_cleaned_noise_inv_I_with_B",
               "sec_A_15hr_41-73_cleaned_noise_inv_I_with_C",
               "sec_A_15hr_41-73_cleaned_noise_inv_I_with_D",
               "sec_B_15hr_41-73_cleaned_noise_inv_I_with_A",
               "sec_B_15hr_41-73_cleaned_noise_inv_I_with_C",
               "sec_B_15hr_41-73_cleaned_noise_inv_I_with_D",
               "sec_C_15hr_41-73_cleaned_noise_inv_I_with_A",
               "sec_C_15hr_41-73_cleaned_noise_inv_I_with_B",
               "sec_C_15hr_41-73_cleaned_noise_inv_I_with_D",
               "sec_D_15hr_41-73_cleaned_noise_inv_I_with_A",
               "sec_D_15hr_41-73_cleaned_noise_inv_I_with_B",
               "sec_D_15hr_41-73_cleaned_noise_inv_I_with_C"]
}

# TODO: Sec B N^-1 is erroneous so we use Sec A N^-1
# TODO: Sec A N^-1 is inverted
old_twoway = {
    'root_data': "/mnt/raid-project/gmrt/kiyo/wiggleZ/corr/",
    'maplist': ["sec_A_15hr_41-69_cleaned_clean_map_I",
              "sec_B_15hr_41-69_cleaned_clean_map_I"],
    'covlist': ["sec_A_15hr_41-69_cleaned_noise_inv_I",
              "sec_A_15hr_41-69_cleaned_noise_inv_I"],
    'multiply_cov': [-1., -1.]
}

def combine_maps(param_dict):
    root_data = param_dict["root_data"]
    maplist = param_dict["maplist"]
    covlist = param_dict["covlist"]

    try: 
        multiply_cov = param_dict["multiply_cov"]
        print "using user-specified covariance multipliers"+repr(multiply_cov)
    except:
        multiply_cov = [1.]*len(covlist)
    mul_cov_list = zip(covlist, multiply_cov) 
    print mul_cov_list

    maps = []
    for tagname in maplist:
        filename = root_data+tagname+".npy"
        print "loading "+filename
        maps.append(algebra.make_vect(algebra.load(filename)))

    weights = []
    for cov_entry in mul_cov_list:
        (tagname, multiplier) = cov_entry 
        filename = root_data+tagname+".npy"
        print "loading " + filename
        print "multiplier " + repr(multiplier)

        # zero out any messy stuff
        raw_weight = algebra.make_vect(algebra.load(filename))
        raw_weight *= multiplier
        raw_weight[raw_weight < 1.e-20] = 0.
        nan_array = np.isnan(raw_weight)
        raw_weight[nan_array] = 0.
        inf_array = np.isinf(raw_weight)
        raw_weight[inf_array] = 0.
        weights.append(raw_weight)

    num_maps = len(maps)
    prodmap = []
    for mapind in range(0, num_maps):
        prodmap.append(maps[mapind]*weights[mapind])

    for mapind in range(1, num_maps):
        prodmap[0] += prodmap[mapind]
        weights[0] += weights[mapind]

    wxc.compressed_array_summary(weights[0], "weight map")
    wxc.compressed_array_summary(prodmap[0], "product map")

    newmap = prodmap[0]/weights[0]

    newweights = weights[0]
    newweights[newweights < 1.e-20] = 0.
    # if the new map is nan or inf, set it and the wieghts to zero
    nan_array = np.isnan(newmap)
    newmap[nan_array] = 0.
    newweights[nan_array] = 0.
    inf_array = np.isinf(newmap)
    newmap[inf_array] = 0.
    newweights[inf_array] = 0.
    wxc.compressed_array_summary(newmap, "new map")
    wxc.compressed_array_summary(newweights, "final weight map")

    return (newmap, newweights, prodmap[0])

if __name__ == '__main__':
    (map, weights, prodmap) = combine_maps(fourway_split) 
    algebra.save("combined_41-73_cleaned_clean_test.npy", map)
    algebra.save("combined_41-73_cleaned_noise_inv_test.npy", weights)
    algebra.save("combined_41-73_cleaned_product_test.npy", prodmap)
