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

