from utils import data_paths
from core import algebra
import numpy as np
import sys

def unique_list(listin):
    seen = set()
    seen_add = seen.add
    uniq = [ x for x in listin if x not in seen and not seen_add(x)]
    return sorted(uniq)

def combine_maps_driver(inputmap_dict, inputweight_dict, output_dict,
                        fullcov=False, datapath_db=None):
    r"""Combine a list of weights, maps specified by their database keys
    """
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    signal_list = []
    weight_list = []
    for mapkey in inputmap_dict:
        signalfile = inputmap_dict[mapkey]
        weightfile = inputweight_dict[mapkey]
        print "loading pair: %s %s" % (signalfile, weightfile)
        signal_list.append(algebra.make_vect(algebra.load(signalfile)))

        if fullcov:
            raw_weight = algebra.make_mat(
                            algebra.open_memmap(weightfile))
            raw_weight = raw_weight.mat_diag()
        else:
            raw_weight = algebra.make_vect(algebra.load(weightfile))

        # zero out any messy stuff
        raw_weight[raw_weight < 1.e-20] = 0.
        raw_weight[np.isnan(raw_weight)] = 0.
        raw_weight[np.isinf(raw_weight)] = 0.
        weight_list.append(raw_weight)

    prodmap = []
    for mapind in range(0, len(signal_list)):
        prodmap.append(signal_list[mapind] * weight_list[mapind])

    print "CHECK THESE: %d %d %d" % (len(signal_list), len(weight_list),
                                     len(prodmap))

    cumulative_product = algebra.zeros_like(prodmap[0])
    cumulative_weight = algebra.zeros_like(prodmap[0])
    for mapind in range(0, len(signal_list)):
        cumulative_product += prodmap[mapind]
        cumulative_weight += weight_list[mapind]

    algebra.compressed_array_summary(cumulative_weight, "weight map")
    algebra.compressed_array_summary(cumulative_product, "product map")

    newmap = cumulative_product / cumulative_weight

    cumulative_weight[cumulative_weight < 1.e-20] = 0.
    cumulative_product[cumulative_weight < 1.e-20] = 0.

    # if the new map is nan or inf, set it and the wieghts to zero
    nan_array = np.isnan(newmap)
    newmap[nan_array] = 0.
    cumulative_product[nan_array] = 0.
    cumulative_weight[nan_array] = 0.
    inf_array = np.isinf(newmap)
    newmap[inf_array] = 0.
    cumulative_product[inf_array] = 0.
    cumulative_weight[inf_array] = 0.
    algebra.compressed_array_summary(newmap, "new map")
    algebra.compressed_array_summary(cumulative_product, "final map * weight")
    algebra.compressed_array_summary(cumulative_weight, "final weight map")

    print output_dict
    algebra.save(output_dict['map'], newmap)
    algebra.save(output_dict['product'], cumulative_product)
    algebra.save(output_dict['weight'], cumulative_weight)
    algebra.save(output_dict['ones'], algebra.ones_like(newmap))

def combine_maps(source_key, combined_key,
                 signal='map', weight='noise_inv', divider_token=";",
                 fullcov=False):
    r"""
    `source_key` is the file db key for the maps to combine
    `combined_key` is the file db key for the combined maps
    `signal` is the tag in the file db entry for the signal maps
    `weight` is the tag in the file db entry for the N^-1 weights
    `fullcov` uses a memory map for large N^-1 and pulls the diagonal
    `divider_token` is the token that divides the map section name
            from the data type e.g. "A_with_B;noise_inv"
    """
    datapath_db = data_paths.DataPath()
    input_fdb = datapath_db.fetch(source_key, intend_read=True,
                                   silent=True)
    output_fdb = datapath_db.fetch(combined_key, intend_write=True,
                                   silent=True)
    (input_fdict, output_fdict) = (input_fdb[1], output_fdb[1])

    # determine the number of foreground cleaning options (modes removed)
    input_map_treatments = []
    input_map_splits = []
    for filekey in input_fdb[0]:
        # the key for a cleaned map is e.g. "C_with_B;map;45modes"
        keyparse = filekey.split(divider_token)
        if len(keyparse) == 3:
            input_map_treatments.append(keyparse[-1])
            input_map_splits.append(keyparse[0])

    input_map_treatments = unique_list(input_map_treatments)
    input_map_splits = unique_list(input_map_splits)

    output_map_treatments = []
    output_map_products = []
    for filekey in output_fdb[0]:
        # the key for a cleaned map is e.g. "map;45modes"
        keyparse = filekey.split(divider_token)
        if len(keyparse) == 2:
            output_map_treatments.append(keyparse[-1])
            output_map_products.append(keyparse[0])

    output_map_treatments = unique_list(output_map_treatments)
    output_map_products = unique_list(output_map_products)

    if input_map_treatments != output_map_treatments:
        print "the source map does not match the requested combined map output"
        sys.exit()

    print input_map_splits, input_map_treatments, output_map_products

    for treatment in input_map_treatments:
        inputmap_dict = {}
        inputweight_dict = {}
        for split in input_map_splits:
            mapkey = "%s;%s;%s" % (split, signal, treatment)
            weightkey = "%s;%s;%s" % (split, weight, treatment)
            inputmap_dict[split] = input_fdict[mapkey]
            inputweight_dict[split] = input_fdict[weightkey]

        output_dict = {}
        for product in output_map_products:
            mapkey = "%s;%s" % (product, treatment)
            output_dict[product] = output_fdict[mapkey]

        #print "-"*80
        #print inputmap_dict
        #print inputweight_dict
        #print output_dict

        combine_maps_driver(inputmap_dict, inputweight_dict, output_dict,
                            fullcov=fullcov, datapath_db=datapath_db)

if __name__ == '__main__':
    combine_maps("GBT_15hr_map_fdgcal_cleaned",
                 "GBT_15hr_map_fdgcal_cleaned_combined")

    combine_maps("GBT_15hr_map_fdgcal_cleaned_sims",
                 "GBT_15hr_map_fdgcal_cleaned_sims_combined")

    combine_maps("GBT_15hr_map_fdgcal_cleaned_noconv",
                 "GBT_15hr_map_fdgcal_cleaned_noconv_combined")

    combine_maps("GBT_15hr_map_fdgcal_cleaned_sims_noconv",
                 "GBT_15hr_map_fdgcal_cleaned_sims_noconv_combined")

    sys.exit()
    combine_maps("GBT_22hr_map_cleaned",
                 "GBT_22hr_map_cleaned_combined")

    combine_maps("GBT_1hr_map_cleaned",
                 "GBT_1hr_map_cleaned_combined")

    # ADD SIMS HERE

    sys.exit()
    combine_maps("GBT_15hr_optimalmap_glued_cleaned",
                 "GBT_15hr_optimalmap_glued_cleaned_combined")

    combine_maps("GBT_15hr_optimalmap_glued_cleaned_sims",
                 "GBT_15hr_optimalmap_glued_cleaned_sims_combined")

    combine_maps("GBT_15hr_optimalmap_glued_cleaned_noconv",
                 "GBT_15hr_optimalmap_glued_cleaned_noconv_combined")

    combine_maps("GBT_15hr_optimalmap_glued_cleaned_sims_noconv",
                 "GBT_15hr_optimalmap_glued_cleaned_sims_noconv_combined")

    sys.exit()
    combine_maps("GBT_15hr_map_cleaned",
                 "GBT_15hr_map_cleaned_combined")

    combine_maps("GBT_15hr_map_cleaned_sims",
                 "GBT_15hr_map_cleaned_sims_combined")

    combine_maps("GBT_15hr_map_cleaned_noconv",
                 "GBT_15hr_map_cleaned_noconv_combined")

    combine_maps("GBT_15hr_map_cleaned_sims_noconv",
                 "GBT_15hr_map_cleaned_sims_noconv_combined")

    combine_maps("GBT_15hr_map_oldcal_cleaned",
                 "GBT_15hr_map_oldcal_cleaned_combined")

    combine_maps("GBT_15hr_map_oldcal_cleaned_sims",
                 "GBT_15hr_map_oldcal_cleaned_sims_combined")

    combine_maps("GBT_15hr_map_oldcal_cleaned_noconv",
                 "GBT_15hr_map_oldcal_cleaned_noconv_combined")

    combine_maps("GBT_15hr_map_oldcal_cleaned_sims_noconv",
                 "GBT_15hr_map_oldcal_cleaned_sims_noconv_combined")

    combine_maps("GBT_15hr_map_fluxcal_cleaned",
                 "GBT_15hr_map_fluxcal_cleaned_combined")

    combine_maps("GBT_15hr_map_fluxcal_cleaned_sims",
                 "GBT_15hr_map_fluxcal_cleaned_sims_combined")

    combine_maps("GBT_15hr_map_fluxcal_cleaned_noconv",
                 "GBT_15hr_map_fluxcal_cleaned_noconv_combined")

    combine_maps("GBT_15hr_map_fluxcal_cleaned_sims_noconv",
                 "GBT_15hr_map_fluxcal_cleaned_sims_noconv_combined")

