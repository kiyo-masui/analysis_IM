from matplotlib import pyplot as plt
#from utils import binning
#import bin_map

# Regarding histograms
def make_hist(dmap, num_bin):
    """
    Plot a histogram in the figure object, using log scale.
    """
    hist, bin_edges, patches = plt.hist(dmap, bins=num_bin, normed=True,
histtype = 'step')
    plt.xscale('log')
    plt.yscale('log')

# Regarding hdf5 catalogs
def dat_mass(adress, group = "Halo_Masses", subgroup = "Halo_Masses"):
    """() -> np.array
    Return masses of objects from the catalog.
    """
    catalog = h5py.File(adress, 'r')
    mass = catalog[group][subroup][:]
    catalog.close()
    return mass

def dat_pos(adress):
    """() -> list of np.array
    Return list of positions (x,y,z) of objects from the catalog.
    """
    catalog = h5py.File(adress, 'r')
    pos_x = catalog['Positions']['x'][:]
    pos_y = catalog['Positions']['y'][:]
    pos_z = catalog['Positions']['z'][:]
    catalog.close()
    return pos_x, pos_y, pos_z

def dat_sum(adress ,redshift = 0.696 , group = "Halo_Masses",
            subgroup = "Halo_Masses"):
    """
    Return numpy array of masses from all catalogs.
    """
    cumul = np.zeros(0)
    adress = '/mnt/raid-project/gmrt/mufma/h5py_catalogs/'
    for list in range(1,66),range(67,77),range(78,101):
        for num in list:
            data = dat_mass(adress + 'simulation%d'%num +
                            '/%.3fhalo_catalog.hdf5'%redshift)
            cumul = cumul + data
    return cumul

# Regarding Tinker mass function
def tinker(directory = "/cita/h/home-2/mufma/code/Tinker/test.dndM"):
    """() -> list of np.array
    Return list of output (mass, dndm) from the Tinker code.
    """
    file = open(directory, "r")
    halo_mf_data = np.genfromtxt(file)
    mass = halo_mf_data[:,0]
    dndm = halo_mf_data[:,1]
    return mass, dndm

