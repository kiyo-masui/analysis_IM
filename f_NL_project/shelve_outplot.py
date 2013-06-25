import shelve
import matplotlib.pyplot as plt
import numpy as np

def open_shelve(filename):
    """
    Open shelve file and return list of data
    """
    file = shelve.open(filename)
    data = file['sim:phys'][1][1]
    points = data['bin_center'][:]
    values = data['binavg'][:]
    file.close()
    return points, values

def plot_r_val(filename_1, filename_2, filename_3):
    """
    Plot r values and output it for analysis
    """
    cross = open_shelve(filename_1)
    delta = open_shelve(filename_2)
    kappa = open_shelve(filename_3)
    r_val = cross[1]/np.sqrt(delta[1]*kappa[1])
    plt.plot(cross[0], r_val, 'b')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
