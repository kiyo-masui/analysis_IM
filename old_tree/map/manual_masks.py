from utils import data_paths
from core import algebra
import numpy as np

def add_manual_mask(source_key, cut_freq_list=None,
                    signal_name='map', noise_inv_name='noise_inv',
                    weight_name='weight', divider_token=";"):
    r"""
    `source_key` is the file db key for the maps to combine
    `signal_name` is the tag in the file db entry for the signal maps
    `noise_inv_name` is the tag in the file db entry for the N^-1 weights
    `weight_name` is the tag in the file db entry for the weights to write out
    `divider_token` is the token that divides the map section name
            from the data type e.g. "A_with_B;noise_inv"
    """
    datapath_db = data_paths.DataPath()
    source_fdb = datapath_db.fetch(source_key, silent=True)
    source_fdict = source_fdb[1]

    # accumulate all the files to combine
    noise_inv_keys = {}
    weight_keys = {}
    signal_keys = {}
    for filekey in source_fdb[0]:
        if divider_token in filekey:
            data_type = filekey.split(divider_token)[1]
            map_section = filekey.split(divider_token)[0]

            if data_type == signal_name:
                signal_keys[map_section] = source_fdict[filekey]

            if data_type == noise_inv_name:
                noise_inv_keys[map_section] = source_fdict[filekey]

            if data_type == weight_name:
                weight_keys[map_section] = source_fdict[filekey]

    for mapkey in signal_keys:
        signal_file = signal_keys[mapkey]
        noise_inv_file = noise_inv_keys[mapkey]
        weight_file = weight_keys[mapkey]
        print "loading pair: %s %s -> %s" % \
                (signal_file, noise_inv_file, weight_file)
        signal_map = algebra.make_vect(algebra.load(signal_file))
        weightmap = algebra.make_vect(algebra.load(noise_inv_file))

        # set the new weights to zero where the N^-1 is small
        # or the signal map is inf or nan
        weightmap[np.isnan(weightmap)] = 0.
        weightmap[np.isinf(weightmap)] = 0.
        weightmap[np.isnan(signal_map)] = 0.
        weightmap[np.isinf(signal_map)] = 0.
        weightmap[weightmap < 1.e-20] = 0.

        if cut_freq_list is not None:
            for cutindex in cut_freq_list:
                weightmap[cutindex, :, :] = 0.

        # could also determine the filename here, outside of the database
        #outputdir = datapath_db.fetch_parent(source_key, return_path=True)
        #weight_out = "%s/%s" % (outputdir, source_key)
        algebra.compressed_array_summary(weightmap, "new weight map")
        algebra.save(weight_file, weightmap)


if __name__ == '__main__':
    #fieldlist = ['15hr', '22hr']
    #fieldlist = ['22hr']
    fieldlist = ['15hr']

    cutlist = [6, 7, 8, 15, 16, 18, 19, 20, 21, 22, 37, 103, 104, 105, 106,
               107, 108, 130, 131, 132, 133, 134, 237, 244, 254, 255]

    augmented = [177, 194, 195, 196, 197, 198, 201, 204, 209, 213, 229]
    for entry in augmented:
        cutlist.append(entry)

    for field in fieldlist:
        for mode_num in range (0, 55, 5):
            source_key = 'GBT_%s_map_cleaned_%dmode' % (field, mode_num)
            add_manual_mask(source_key, cut_freq_list=cutlist)
