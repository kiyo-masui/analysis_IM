import scipy as sp
import cPickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from utils import data_paths

def plot_six_svd_amplitudes(vals, fig_filename):
    n_vals = len(vals)

    plt.xlim((0.8,300))
    plt.loglog(abs(sp.sort(-vals[0] / float(n_vals))), 'b.',
               abs(sp.sort(-vals[1] / float(n_vals))), 'g.',
               abs(sp.sort(-vals[2] / float(n_vals))), 'r.',
               abs(sp.sort(-vals[3] / float(n_vals))), 'c.',
               abs(sp.sort(-vals[4] / float(n_vals))), 'm.',
               abs(sp.sort(-vals[5] / float(n_vals))), 'y.',
               linestyle='None')
    plt.title("Mode amplitudes")
    plt.xlabel("Mode number")
    plt.ylabel("Amplitude")
    plt.savefig(fig_filename, dpi=200)
    plt.clf()
    #print 'Mean noise: ', sp.sum(vals) / n_vals
    #print 'Largest eigenvalues/vals: ',
    #print sp.sort(vals / n_vals)[-10:]

def plot_six_svd_modes(vals, fig_filename):
    n_vals = len(vals)

    plt.xlim((0,231))
    plt.plot(vals[0], 'b-',
               vals[1], 'g-',
               vals[2], 'r-',
               vals[3], 'c-',
               vals[4], 'm-',
               vals[5], 'y-',
            )
    plt.title("Modes")
    plt.xlabel("Frequency index")
    plt.ylabel("Mode amplitude")
    plt.savefig(fig_filename, dpi=200)
    plt.clf()
    #print 'Mean noise: ', sp.sum(vals) / n_vals
    #print 'Largest eigenvalues/vals: ',
    #print sp.sort(vals / n_vals)[-10:]


def plot_svd_series():
    root = "/mnt/raid-project/gmrt/eswitzer/GBT/"
    #root += "cleaned_maps/freq_slices_refactor_tests/"
    root += "cleaned_maps/15hr/"
    svdblob = cPickle.load(open(root + "svd_info.pkl", "r"))
    (amp_ind, mode0_ind, mode1_ind) = (0,1,2)

    amp_set = [ svdblob[pair_ind][amp_ind] for pair_ind in range(0,6) ]
    plot_six_svd_amplitudes(amp_set, "mode_amplitudes.png")

    moderange = range(0,15)
    filelist = [ "svd_mode_%d.png" % mode_index for mode_index in moderange ]
    for (filename, index) in zip(filelist, moderange):
        print filename
        mode_set = [ svdblob[pair_ind][mode0_ind][index] for pair_ind in range(0,6) ]
        plot_six_svd_modes(mode_set, filename)

plot_svd_series()
