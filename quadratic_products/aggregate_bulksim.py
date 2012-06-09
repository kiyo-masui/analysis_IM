import numpy as np
import numpy.ma as ma
from quadratic_products import power_spectrum as ps
import glob
from plotting import plot_slice
import shelve
from kiyopy import parse_ini

aggregatesummary_init = {
        "directory": "dir",
        "basefile": "file",
        "outfile": "file"
    }

aggregatesummary_prefix = 'as_'

class AggregateSummary(object):
    """take a directory of simulation outputs and merge them into a single
    numpy array/shelve file"""
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          aggregatesummary_init,
                                          prefix=aggregatesummary_prefix)

    def execute(self, processes):
        """produce a list of files to combine and run"""
        outfile = self.params["outfile"]
        filelist = glob.glob('%s/%s*.shelve' % \
                             (self.params['directory'],
                              self.params['basefile']))

        self.produce_summary(filelist, outfile)

        bulksubdir = "optimal_762_beam_meansub/"
        basename = "GBT_15hr_optimalmap_selfcal_762_sim_"
        outfile = bulksim_dir + basename + "_beam_meansub.shelve"

    def produce_summary(self, filelist, outfile, debug=False):
        num_sim = len(filelist)
        print "number of simulations to aggregate: %d" % num_sim

        sim_toread = ps.PowerSpectrum(filelist[0])
        # run this here just to get the 2D->1D bin centers
        sim_toread.convert_2d_to_1d()

        k1d = sim_toread.bin_center_1d
        k1d_from_2d = sim_toread.bin_center_1d_from_2d
        num_k1d = sim_toread.num_k1d
        num_k1d_from_2d = sim_toread.num_k1d_from_2d

        (kx, ky) = (sim_toread.bin_center_x, sim_toread.bin_center_y)
        (num_kx, num_ky) = (sim_toread.num_kx, sim_toread.num_ky)

        print "k_1d bins: %d, %d, kx bins: %d, ky bins: %d" % \
               (num_k1d, num_k1d_from_2d, num_kx, num_ky)
        del sim_toread

        trial_array_1d = np.zeros((num_sim, num_k1d))
        trial_array_1d_from_2d = np.zeros((num_sim, num_k1d_from_2d))
        trial_array_2d = np.zeros((num_sim, num_kx, num_ky))

        for (simfile, index) in zip(filelist, range(num_sim)):
            print simfile, index
            sim_toread = ps.PowerSpectrum(simfile)
            sim_toread.convert_2d_to_1d()

            agg1d = sim_toread.agg_stat_1d_pwrspec()
            agg1d_from_2d = sim_toread.agg_stat_1d_pwrspec(from_2d=True)
            agg2d = sim_toread.agg_stat_2d_pwrspec()

            trial_array_1d[index, :] = agg1d['0modes']['mean']
            trial_array_1d_from_2d[index, :] = agg1d_from_2d['0modes']['mean']
            trial_array_2d[index, :, :] = agg2d['0modes']['mean']
            if debug:
                nanarray = np.isnan(trial_array_2d[index, :, :])
                print "is nan", np.sum(nanarray)

            del sim_toread

        outshelve = shelve.open(outfile, "n", protocol=-1)
        outshelve["k1d"] = k1d
        outshelve["k1d_from_2d"] = k1d_from_2d
        outshelve["kx"] = kx
        outshelve["ky"] = ky
        outshelve["pk_1d"] = trial_array_1d
        outshelve["pk_1d_from_2d"] = trial_array_1d_from_2d
        outshelve["pk_2d"] = trial_array_2d
        outshelve.close()


def aggregate_1d_statistics(summary_shelvefile, from2d = False):
    r"""Find the mean, std, cov, corr etc of 1D p(k)'s"""
    summary = shelve.open(summary_shelvefile, "r")

    if from2d:
        k_1d = summary["k1d_from_2d"]
        trial_array_1d = summary["pk_1d_from_2d"]
    else:
        k_1d = summary["k1d"]
        trial_array_1d = summary["pk_1d"]

    summary.close()

    mtrial_array_1d = ma.array(trial_array_1d, mask=np.isnan(trial_array_1d))

    mean1d = np.ma.mean(mtrial_array_1d, axis=0)
    std1d = np.ma.std(mtrial_array_1d, axis=0, ddof=1)
    corr1d = np.ma.corrcoef(np.ma.transpose(mtrial_array_1d), ddof=1)
    cov1d = np.ma.cov(np.ma.transpose(mtrial_array_1d), ddof=1)

    logk_1d = np.log10(k_1d)
    plot_slice.simpleplot_2D("sim_cov_1d.png", corr1d, logk_1d, logk_1d,
                         ["logk", "logk"], 1., "1D covariance", "corr")

    for (k, pk, err) in zip(k_1d, mean1d, std1d):
        print k, pk, err


def aggregate_2d_statistics(summary_shelvefile):
    summary = shelve.open(summary_shelvefile, "r")
    k_x = summary["k_x"]
    k_y = summary["k_y"]
    trial_array_2d = summary["pk_2d"]
    summary.close()

    mtrial_array_2d = ma.array(trial_array_2d, mask=np.isnan(trial_array_2d))

    mean2d = np.ma.mean(mtrial_array_2d, axis=0)
    #mean2d[np.isnan(mean2d)] = 0.
    mean2d[np.isnan(mean2d)] = 1.e-16
    np.set_printoptions(threshold='nan')
    print mean2d
    logkx = np.log10(kx)
    logky = np.log10(ky)
    # TODO: why not logky?
    plot_slice.simpleplot_2D("sim_2d.png", mean2d*10**12., logkx, logky,
                             ["logkx", "logky"], 1., "2D power", "logP(k)",
                             logscale=True)

    # can use C or F to do column or row-major
    num_flat = num_kx * num_ky
    flat_axis = range(num_flat)
    mtrial_array_2d_flat = np.ma.reshape(mtrial_array_2d, (num_sim, num_flat))
    print mtrial_array_2d_flat.shape
    #corr2d = np.corrcoef(np.ma.transpose(mtrial_array_2d_flat), ddof=1)
    corr2d = np.ma.cov(np.ma.transpose(mtrial_array_2d_flat), ddof=1)
    print corr2d.shape
    plot_slice.simpleplot_2D("sim_2d_corr.png", corr2d, flat_axis, flat_axis,
                             ["k", "k"], 1., "2D power corr", "corr")

if __name__ == "__main__":
    bulksim_dir = "/mnt/raid-project/gmrt/eswitzer/GBT/simulations/bulk/"

    summary_shelvefile = bulksim_dir + \
                         "GBT_15hr_optimalmap_selfcal_762_nobeam.shelve"

    # GBT_15hr_optimalmap_selfcal_762_beam_meansub.shelve

    aggregate_1d_statistics(summary_shelvefile, from2d = False)
    aggregate_1d_statistics(summary_shelvefile, from2d = True)
