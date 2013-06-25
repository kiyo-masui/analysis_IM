import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
# useful function is array.flatten()


def make_hist(catalog, num_bin):
    """
    Output histogram data for the plot
    """
    hubble = 0.7
    volume = 8.e9
    hist, bin_edges = np.histogram(catalog, bins=num_bin, density=False)
    mass = (np.roll(bin_edges, -1) + bin_edges)/2.
    mass = mass[:len(hist)]*hubble
    delta_M = np.roll(bin_edges, -1) - bin_edges
    delta_M = delta_M[:len(hist)]*hubble
    return mass, hist/(delta_M*volume)

def catalog_reader(catalog_name, group, subgroup):
    """
    Returns array of data from the catalog[group][subgroup]
    """
    catalog = h5py.File(catalog_name)
    data = catalog[group][subgroup][:]
    catalog.close()
    return data

def tinker_reader():
    """
    Read tinker halo mass function.
    """
    directory = "/cita/h/home-2/mufma/code/Tinker/test0_696.dndM"
    file = open(directory, "r")
    halo_mf_data = np.genfromtxt(file)
    mass = halo_mf_data[0:99,0]
    dndM = halo_mf_data[0:99,1]
    file.close()
    return mass, dndM

def append_joe_catalog(redshift):
    """
    Read catalog for a given redshift. Put all info in one array.
    """
    cumul = np.zeros(0)
    for list in range(1,66),range(67,77),range(78,101):
        for num in list:
            data = catalog_reader('/mnt/raid-project/gmrt/mufma/h5py_catalogs/'
                                  + 'simulation%d'%num
                                  +"/%.3fhalo_catalog.hdf5"%redshift,
                                  'Halo_Masses','Halo_Masses')
        cumul = np.append(cumul, data)
    print "Summation of catalogs finished"
    return cumul

def append_jd_catalog(redshift):
    """
    Read catalog for a given redshift. Put all info in one array.
    """
    cumul = np.zeros(0)
    for num in [0,1,2,3,4,5,6,7]:
        data = catalog_reader('/mnt/raid-project/gmrt/mufma/JD_h5py_catalogs/0.800halo_catalog000%d.hdf5'%num,'Halo_Masses','Halo_Masses')
        cumul = np.append(cumul, data)
    return cumul

def plot(mass_1, dndM_1, mass_2, dndM_2):
    plt.plot(mass_1, dndM_1, 'r')
    plt.plot(mass_2, dndM_2, 'b')
    plt.xscale('log')
    plt.yscale('log')

if __name__ == '__main__':
    t_mass, t_dndM = tinker_reader()
    cumul = append_jd_catalog(0.800)
    #cumul = append_joe_catalog(0.696)
    h_mass, h_dndM = make_hist(cumul, 1000)
    f = interp1d(h_mass, h_dndM)
    print "h_max", h_mass.max(), ""
    print "t_max", t_mass.max(), ""
    y_axis = f(t_mass)
    plot(t_mass, t_dndM, t_mass, y_axis)
    #plot(h_mass, h_dndM, t_mass, t_dndM)
    plt.show()
