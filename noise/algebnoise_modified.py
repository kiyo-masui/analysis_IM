from core import algebra
import numpy as np
from math import sqrt
import copy
import random
import matplotlib.pyplot as plt


def make_real_noise(file_directory, file_name):
    # opening files
    weight_filename = open(file_directory + file_name, "r")
    # making a matrix
    weight = algebra.make_vect(algebra.load(weight_filename)) + 1e-9
    # creating an array representing frequencies of the slices
    freq_axis = weight.get_axis("freq")
    #calculation of frequency binning
    array = np.roll(freq_axis, 1) - freq_axis
    delta_freq = array[4]
    # calculating the total time in weight map units
    total_units = np.sum(weight[10, :, :])
    tot = np.sum(weight[:, :, :])
    print total_units/tot
    # creating a matrix for noise connecting to it an info package(row,col,...)
    coef = 30000.0/np.sqrt(delta_freq*(360000.0/total_units))
    Tnoise = np.random.standard_normal(weight.shape)
    Tnoise = Tnoise * coef / np.sqrt(weight)
    Tnoise.copy_axis_info(weight)
    # saving cube to the directory
    file = open("/tmp/mufma/data/DeltaT.npy", "w")
    algebra.save(file, Tnoise)
    file.close()

def calc_typical_time(file_directory, file_name, int_time=36e4):
    '''
    Returns typical value of weight from the file_name, from the region in the
    centre of the map.
    '''
    weight_filename = open(file_directory + file_name, "r")
    weight = algebra.make_vect(algebra.load(weight_filename))
    max = weight.max()
    weight_slice = weight[10, :, :]
    mean = np.sum(weight_slice)/(len(weight_slice)*len(weight_slice[0]))
    weight_slice = weight_slice.flatten()
    slice = weight_slice[np.where(weight_slice>mean)]
    return np.sum(slice)/(len(slice))*int_time/np.sum(weight_slice)

def calc_std_deviation(file_directory, file_name, delta_T=3e4,
                       int_time=36e4):
    '''
    Returns standard deviation of the weight map from the file_name
    '''
    weight_filename = open(file_directory + file_name, "r")
    weight = algebra.make_vect(algebra.load(weight_filename)) + 1e-9
    freq_axis = weight.get_axis("freq")
    array = np.roll(freq_axis, 1) - freq_axis
    delta_freq = array[4]
    total_units = np.sum(weight[10, :, :])
    coef = delta_T/np.sqrt(delta_freq*(int_time/total_units))
    array = np.zeros(200)
    for index in range(len(array)):
        Tnoise = np.random.standard_normal(weight.shape)
        Tnoise = Tnoise * coef / np.sqrt(weight)
        Tnoise = Tnoise.flatten()
        Tnoise = Tnoise[np.where(abs(Tnoise)<20.-index*0.1)]
        array[index] = Tnoise.std()
    return array


if __name__ == '__main__':
    directory = "/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_oldcal/"
    name = "sec_A_15hr_41-90_noise_weight_I.npy"
    print calc_std_deviation(directory, name)
    print calc_typical_time(directory, name)
    #make_real_noise(directory, name)
