import numpy as np
import random


def make_noise_map(delta_T, dev, bins_x, bins_y, bins_z):
    '''(list) -> ndarray
    Precondition: list of floats
    Returns 3D array shaped as (bins_x, bins_y, bins_z) with numbers drawn from
    Gaussian(mean = 0, std.deviation = dev) multiplied by delta_T
    '''
    return array = delta_T*np.random.normal(0, dev, size = (bins_x, bins_y, bins_z))

def calculate_skewness(array):
    '''(array) -> number
    Returns skewness of an array of numbers
    '''
    N = len(array)
    mean1 = np.sum(array)/N
    mean2 = np.sum(array**2)/N
    mean3 = np.sum(array**3)/N
    top = (mean3 - 3*mean1*mean2 + 2*mean1**3)/N
    down = ((mean2 - mean1**2)/N)**(1.5)
    return top/down

def calculate_mean_std(array):
    '''(array) -> list
    Returns mean and standard deviation of the array.
    '''
    N = len(array)
    mean1 = np.sum(array)/N
    mean2 = np.sum(array**2)/N
    std = (mean2 - mean1**2)/N
    return mean1,std


