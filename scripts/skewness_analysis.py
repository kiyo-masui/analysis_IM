import numpy as np
from scipy import stats
import random
import h5py
import bin_map

def make_noise_map(std, bins):
    '''(list) -> ndarray
    Precondition: list of floats
    Return 3D array shaped as (bins_x, bins_y, bins_z) with numbers drawn from
    Gaussian(mean = 0, std.deviation = dev) multiplied by delta_T
    '''
    array = np.random.normal(0, std, size = (bins, bins, bins))
    return array

def noise_skewness(n, bins):
    '''(int, int) -> list
    Return list of mean and std from n realisations of noise
    '''
    array = np.zeros(n)
    for index in range(len(array)):
        map = make_noise_map(3.5, bins)
        map = map.flatten()
        array[index] = stats.skew(map)
    list = array.mean(), array.std()
    return list

def simulation_skewness(catalog, bins, group_name, subgroup_name):
    '''(h5py) -> float
    Return skewness of a simulation map
    '''
    map = bin_map.data_binning(catalog, bins, group_name, subgroup_name)
    map = map.flatten()
    skewness = stats.skew(map)
    return skewness

def noise_add_sim(catalog, bins, group_name, subgroup_name):
    '''(h5py) -> ndarray
    Return sum of the representaion noise map and simulation map from the
    catalog
    '''
    noise_map = make_noise_map(3.5, bins)
    sim_map = bin_map.data_binning(catalog, bins, group_name, subgroup_name)
    total_map = sim_map + noise_map
    return total_map

def noise_sim_skewness(catalog, bins, group_name, subgroup_name, n):
    '''(h5py) -> ndarray
    Return sum of the representaion noise map and simulation map from the
    catalog
    '''
    sim_map = bin_map.data_binning(catalog, bins, group_name, subgroup_name)
    array = np.zeros(n)
    for index in range(len(array)):
        map = make_noise_map(3.5, bins-1)
        map = map + sim_map
        map = map.flatten()
        array[index] = stats.skew(map)
    list = array.mean(), array.std()
    return list

if __name__ == '__main__':
    directory = '/cita/h/home-2/mufma/code/analysis_IM/Tryout.hdf5'
    bins = 30
    print simulation_skewness('/cita/h/home-2/mufma/code/analysis_IM/Tryout.hdf5',
                              bins + 1, '21_brightness_probe', '21_probe')
    print noise_skewness(100, bins)
    print noise_sim_skewness('/cita/h/home-2/mufma/code/analysis_IM/Tryout.hdf5', bins + 1,
'21_brightness_probe', '21_probe', 100)
