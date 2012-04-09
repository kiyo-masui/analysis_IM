from optparse import OptionParser
from kiyopy import parse_ini
from correlate import pair_set

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
               'output_root': "local_test",
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
               'subtract_inputmap_from_sim': False
               }
prefix = 'fs_'

def main():
    r"""main command-line interface"""

    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")

    parser.add_option("-n", "--sim-num",
                      action="store_const",
                      dest="sim_number",
                      default="0",
                      help="Simulation number to inject")

    parser.add_option("-s", "--subtract-signal",
                      action="store_true",
                      dest="sub_signal",
                      default=False,
                      help="Subtracted the cleaned signal from the sim?")

    (options, args) = parser.parse_args()
    optdict = vars(options)

    if len(args) != 1:
        parser.error("wrong number of arguments")

    print options
    print args

    params = parse_ini.parse(args[0], params_init,
                             prefix=prefix)

    params["simfile"] = '/mnt/raid-project/gmrt/eswitzer/GBT/simulations/15hr_oldmap_str/sim_beam_000.npy'
    print params

    pair_set.PairSet(params_dict=params).execute()


if __name__ == "__main__":
    main()
