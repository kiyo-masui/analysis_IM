from utils import data_paths
from core import algebra
import numpy as np
import sys


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
                 signal='map', weight='noise_inv', divider=";",
                 fullcov=False):
    r"""
    `source_key` is the file db key for the maps to combine
    `combined_key` is the file db key for the combined maps
    `signal` is the tag in the file db entry for the signal maps
    `weight` is the tag in the file db entry for the N^-1 weights
    `fullcov` uses a memory map for large N^-1 and pulls the diagonal
    `divider` is the token that divides the map section name
            from the data type e.g. "A_with_B;noise_inv"
    """
    datapath_db = data_paths.DataPath()
    input_fdb = datapath_db.fetch(source_key, intend_read=True,
                                   silent=True)
    output_fdb = datapath_db.fetch(combined_key, intend_write=True,
                                   silent=True)
    (input_fdict, output_fdict) = (input_fdb[1], output_fdb[1])

    # determine the number of foreground cleaning options (modes removed)
    input_map_cases = datapath_db.fileset_cases(source_key,
                      "pair;type;treatment", divider=divider)
    output_map_cases = datapath_db.fileset_cases(combined_key,
                      "type;treatment", divider=divider)

    if input_map_cases['treatment'] != output_map_cases['treatment']:
        print "the source map does not match the requested combined map output"
        sys.exit()

    print input_map_cases['pair'], input_map_cases['treatment'], \
          output_map_cases['type']

    for treatment in input_map_cases['treatment']:
        inputmap_dict = {}
        inputweight_dict = {}
        for split in input_map_cases['pair']:
            mapkey = "%s;%s;%s" % (split, signal, treatment)
            weightkey = "%s;%s;%s" % (split, weight, treatment)
            inputmap_dict[split] = input_fdict[mapkey]
            inputweight_dict[split] = input_fdict[weightkey]

        output_dict = {}
        for product in output_map_cases['type']:
            mapkey = "%s;%s" % (product, treatment)
            output_dict[product] = output_fdict[mapkey]

        #print "-"*80
        #print inputmap_dict
        #print inputweight_dict
        #print output_dict

        combine_maps_driver(inputmap_dict, inputweight_dict, output_dict,
                            fullcov=fullcov, datapath_db=datapath_db)

def wrap_combine(basemap, skipsims=False, skipmaps=False):
    if not skipmaps:
        combine_maps("%s_cleaned" % basemap,
                     "%s_cleaned_combined" % basemap)

        combine_maps("%s_cleaned_noconv" % basemap,
                     "%s_cleaned_noconv_combined" % basemap)

    if not skipsims:
        combine_maps("%s_cleaned_sims" % basemap,
                     "%s_cleaned_sims_combined" % basemap)

        combine_maps("%s_cleaned_sims_noconv" % basemap,
                     "%s_cleaned_sims_noconv_combined" % basemap)


if __name__ == '__main__':
    wrap_combine("GBT_15hr_optimalmap_mapv2fdgcal", skipmaps=True)
    wrap_combine("GBT_15hr_optimalmap_mapv2oldcal", skipmaps=True)
    #wrap_combine("GBT_15hr_optimalmap_fluxpolcal")
    #wrap_combine("GBT_22hr_map_fluxpolcal")
    #wrap_combine("GBT_15hr_map_fluxpolcal")
    #wrap_combine("GBT_1hr_map_fluxpolcal")
    #wrap_combine("GBT_15hr_map_fdgcal")
    #wrap_combine("GBT_15hr_map_oldcal")
    #wrap_combine("GBT_15hr_map_fluxcal")
