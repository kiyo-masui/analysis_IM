import h5py
import numpy as np
import matplotlib.pyplot as plt
from core import algebra
from utils import binning
import bin_map
import skewness_analysis

def make_hist(catalog, num_bin, name, params=['X-axis', 'Y-axis', 'Title'],
              set_logx=False, set_logy=False, savefig=True,
              plotdir="/cita/d/www/home/mufma/plots/"):
    '''
    Make histogram of the map
    '''
    print "Starting to plot histogram"
    #weight = np.zeros(len(map))
    #for k in range(len(weight)):
    #    weight[k] = 1.e-20
    hist, bin_edges = np.histogram(catalog, bins=num_bin)
    mass = np.zeros(len(hist))
    hist = np.array(hist,dtype='Float64')
    for k in range(len(hist)):
        mass[k] = (bin_edges[k+1] + bin_edges[k])/2.
        num = bin_edges[k+1] - bin_edges[k]
        hist[k] = hist[k]/(num*1.7e6)
    plt.plot(mass,hist,linestyle='steps')
    y_min = hist.min()
    y_max = hist.max()
    plt.xlabel(params[0])
    plt.ylabel(params[1])
    plt.title(params[2])
    #plt.ylim(y_min,2*y_max)
    plt.ylim(1e-24,1e2)
    plt.xlim(1e12,1.5e16)
    if set_logx == True:
        plt.xscale('log')
    if set_logy == True:
        plt.yscale('log')
    if savefig == True:
        plt.savefig(plotdir + name, format="png")

def open_data_hdf5(catalog_name, group, subgroup):
    '''
    Returns array of data from the catalog[group][subgroup]
    '''
    catalog = h5py.File(catalog_name)
    data = catalog[group][subgroup][:]
    catalog.close()
    return data

if __name__ == '__main__':
    directory = "/cita/h/home-2/mufma/code/Tinker/test.dndM"
    file = open(directory, "r")
    halo_mf_data = np.genfromtxt(file)
    mass = np.zeros(len(halo_mf_data))
    dndM = np.zeros(len(halo_mf_data))
    plot = np.zeros(len(halo_mf_data))
    for k in range(len(halo_mf_data)):
        mass[k] = halo_mf_data[k][0]
        dndM[k] = halo_mf_data[k][1]
        plot[k] = mass[k]*dndM[k]
    summ = 0
    for k in range(87,99):
        summ = summ + dndM[k]*(mass[k+1]-mass[k])
    print summ
    #plt.plot(mass,plot)
    #plt.xlim(1e10,1e16)
    #plt.yscale('log')
    #plt.xscale('log')
    #plt.show()
    cumul = np.zeros(0)
    for list in range(1,66),range(67,77),range(78,101):
        for num in list:
            data = open_data_hdf5('/mnt/raid-project/gmrt/mufma/h5py_catalogs/'
                                  + 'simulation%d'%num +"/0.696halo_catalog.hdf5",                                      'Halo_Masses','Halo_Masses')
        cumul = np.append(cumul, data)
    print "Opened catalogs"
    
    make_hist(cumul, 100, "aver_halo_vs_tinker",
              params=['Masses','Diff. distrib.','Tinker vs Halo_aver'],
              savefig=False)
    
    plt.plot(mass, dndM, linestyle = 'steps')
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig("tink.png")
    plt.show()
    

    """
    data_HI = open_data_hdf5('/cita/h/home-2/mufma/code/analysis_IM/Tryout.hdf5'                             ,'HI_Masses', 'HI_Masses')
    binned_map = bin_map.data_binning('/cita/h/home-2/mufma/code/\
                                      analysis_IM/Tryout.hdf5',41.,
                                      'HI_Masses', 'HI_Masses')

    binned_map = binned_map.flatten()
    binned_map = binned_map[np.where(binned_map!=0)]
    make_hist([binned_map,data_HI], 500, 'Binned_HI.png',
              ['M_solar','Number of clusters','HI PDF'],
              set_logy = True, set_logx = False)
    """

# THIS IS PART FOR NOISE(1) + BIN_SIM(1) HISTOGRAMS
"""
noise = skewness_analysis.make_noise_map(3.5, 40)
bin_sim = bin_map.data_binning('/cita/h/home-2/mufma/code/analysis_IM/Tryout.hdf5',41.,'21_brightness_probe','21_probe')
sim_noise = bin_sim + noise
bin_sim = bin_sim.flatten()
noise = noise.flatten()
sim_noise = sim_noise.flatten()
make_hist([noise], 100, 'PDF_sim.png',
          ['mK', 'Number', 'PDF'],
          set_logy = True)
"""

