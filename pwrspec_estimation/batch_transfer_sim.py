from optparse import OptionParser
from kiyopy import parse_ini
from correlate import pair_set
from map import combine_maps
from utils import data_paths

# write function to cache noise_inv diag operation (memoize?)
# initialize pair_set

# consider writing the whole thing as a shell script which takes 1 argument:
# the number of the sim to add

# the output root

# how do we parallelize this? least effort is to be set up for tpb160. can
# either write a shell script and run a bunch in parallel, or spawn
# multiprocessing; same difference. in the end, need a function which takes
# basically one value and an ini file. this function prepares a param for the
# extened version of pairset. it then calls the map combiner and the
# cross-power estimator.

# lower prioritry
# add optparse to pairset (run at the command prompt)

params_init = {
               'output_root': "./data_test",
               'map1': 'GBT_15hr_map',
               'map2': 'GBT_15hr_map',
               'noise_inv1': 'GBT_15hr_map',
               'noise_inv2': 'GBT_15hr_map',
               'freq_list': (),
               'lags': (0.1, 0.2),
               'convolve': True,
               'factorizable_noise': True,
               'sub_weighted_mean': True,
               'modes': [10, 15],
               'no_weights': False,
               'SVD_root': None,
               'regenerate_noise_inv': False,
               'subtract_inputmap_from_sim': False,
               'subtract_sim_from_inputmap': False
               }
prefix = 'fs_'

def run_sim_batch(params, optdict):
    print optdict
    datapath_db = data_paths.DataPath()
    params["simfile"] = datapath_db.fetch(optdict['sim_dbkey'],
                                          intend_read=True,
                                          pick=optdict['sim_number'])

    print "using simulation: %s" % params["simfile"]

    if optdict['sub_map']:
        print "subtracting the cleaned map from the cleaned map+sim"
        params["subtract_inputmap_from_sim"] = True

    if optdict['sub_sim']:
        print "subtracting the cleaned sim from the cleaned map+sim"
        params["subtract_sim_from_inputmap"] = True

    #pair_set.PairSet(params_dict=params).execute()

    sim_id = "sim%s" % optdict['sim_number']
    combine_maps.combine_maps("batch_cleaned_cache", "GBT_15hr_oldcal_corrfg",
                              batchsim=sim_id)


def main():
    r"""main command-line interface"""

    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")

    parser.add_option("-d", "--sim-dbkey",
                      action="store",
                      dest="sim_dbkey",
                      default="sim_15hr_oldmap_str_beam",
                      help="db key of simulation to inject")

    parser.add_option("-n", "--sim-num",
                      action="store",
                      dest="sim_number",
                      default="0",
                      help="Simulation number to inject")

    parser.add_option("-m", "--subtract-map",
                      action="store_true",
                      dest="sub_map",
                      default=False,
                      help="Subtracted the cleaned map from the sim?")

    parser.add_option("-s", "--subtract-sim",
                      action="store_true",
                      dest="sub_sim",
                      default=False,
                      help="Subtracted the cleaned sim from the map?")

    (options, args) = parser.parse_args()
    #print options, args
    optdict = vars(options)

    if len(args) != 1:
        parser.error("wrong number of arguments")

    params = parse_ini.parse(args[0], params_init,
                             prefix=prefix)

    run_sim_batch(params, optdict)


if __name__ == "__main__":
    main()
