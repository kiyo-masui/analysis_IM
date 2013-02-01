r"""
Cross-calibrate a large set of maps
"""
from utils import data_paths as dp
from correlate import corr_estimation as ce
from core import algebra
from correlate import map_pair

# function to convert list of big covs to diagonal only

def noise_inv_to_weight(uncal_weightlist):
        noise_inv = algebra.make_mat(
                                    algebra.open_memmap(filename, mode='r'))
        self.noisefiledict[filename] = noise_inv.mat_diag()
        algebra.save(filename_diag, self.noisefiledict[filename])

def map_pair_cal(uncal_maplist, uncal_weightlist, in_path, out_path,
                 convolve=False, factorizable_noise=True,
                 sub_weighted_mean=True):

    # load maps into pairs
    for mapname, noisename in zip(uncal_maplist, uncal_weightlist):
        pair = map_pair.MapPair(map1, map2,
                                noise_inv1, noise_inv2,
                                self.freq_list)

        pair.set_names(pdict['tag1'], pdict['tag2'])

        pair.lags = self.lags
        pair.params = self.params
        self.pairs[pairitem] = pair

        (corr, counts) = pair.correlate(pair.lags, speedup=True)
        svd_info = ce.get_freq_svd_modes(corr, len(self.freq_list))
        leftmode = svd_info[1][:0]
        rightmode =svd_info[2][:0]

    # write out the maps and noise
        algebra.save(map1_file, pair.map1)
        algebra.save(map2_file, pair.map2)
        algebra.save(noise_inv1_file, pair.noise_inv1)
        algebra.save(noise_inv2_file, pair.noise_inv2)

if __name__ == '__main__':
    if len(sys.argv) == 2:
        PairSet(str(sys.argv[1])).execute()
        #PairSet(str(sys.argv[1])).clean_maps()
    else:
        print 'Need one argument: parameter file name.'


