import scipy as sp
import matplotlib.pyplot as plt
import numpy.ma as ma
from numpy import random

# TODO: extend this to a correlation object returned by correlate (that can
# plot itself
# TODO: extend this to allow different freq axes per correlation pair
class CorrelationPlot(object) :
    def __init__(self, lags, corr, freq, freq_axis):
        """
        freq is the set of frequency indices considered
        freq_axis is e.g. Map1.get_axis('freq')
        """
        self.corr = corr
        self.lags = sp.array(lags)
        self.freq_inds = sp.array(freq, dtype=int)
        self.freq1 = freq_axis[self.freq_inds]
        self.freq2 = freq_axis[self.freq_inds]
        self.norms = corr[:, :, 0].diagonal()

    def plot(self, finds1, finds2):
        markers = ['o', 's', '^', 'h', 'p', 'v']
        colours = ['b', 'g', 'r', 'c', 'm', 'y']
        nmarkers = len(markers)
        ncolours = len(colours)

        for find1 in finds1:
            # Plot the correlations.
            plt.figure()
            for ii, find2 in enumerate(finds2):
                plt.plot(self.lags, self.corr[find1, find2, :],
                         linestyle='None',
                         marker=markers[ii%nmarkers],
                         markerfacecolor=colours[ii%ncolours])
            plt.xlabel('lag (degrees)')
            plt.ylabel('normalized correlation with ' +
                       str(int(round(self.freq1[find1] / 1.e6))) + ' MHz bin')
            #Plot the normalizations to give an absolute scale.
            plt.figure()
            for ii, find2 in enumerate(finds2):
                plt.plot(self.freq2[find2] / 1.e6, self.norms[find1, find2],
                         linestyle='None',
                         marker=markers[ii%nmarkers],
                         markerfacecolor=colours[ii%ncolours])
            plt.xlabel('Frequency (MHz)')
            plt.ylabel('Normalization of correlation with ' +
                       str(int(round(self.freq1[find1] / 1.e6))) +
                       ' MHz bin (K^2)')


    def plot_freq(self, norms=False, lag_ind=0):
        # Set up binning.
        nf = len(self.freq_inds)
        n_bins = 12
        factor = 1.5
        start = 2.1e6
        freq_diffs = sp.empty(n_bins)
        freq_diffs[0] = start
        for ii in range(1, n_bins):
            freq_diffs[ii] = factor * freq_diffs[ii - 1]
        n_diffs = len(freq_diffs)
        # Allowcate memory.
        corrf = sp.zeros(n_diffs)
        countsf = sp.zeros(n_diffs, dtype=int)
        for ii in range(nf):
            for jj in range(nf):
                if norms:
                    thiscorr = self.corr[ii, jj, lag_ind] * self.norms[ii, jj]
                else:
                    thiscorr = self.corr[ii, jj, lag_ind]
                df = abs(self.freq1[ii] - self.freq2[jj])
                for kk in range(1, n_diffs - 1):
                    if (abs(freq_diffs[kk] - df) <= abs(freq_diffs[kk - 1] - df)
                        and abs(freq_diffs[kk] - df) < abs(freq_diffs[kk + 1] - df)):
                        d_ind = kk
                if abs(freq_diffs[0] - df) < abs(freq_diffs[1] - df):
                    d_ind = 0
                if abs(freq_diffs[-1] - df) <= abs(freq_diffs[-2] - df):
                    d_ind = n_diffs - 1
                corrf[d_ind] += thiscorr
                countsf[d_ind] += 1
        corrf /= countsf
        plt.semilogx(freq_diffs / 1e6, corrf * 1e6, linestyle='None',
                 marker='o')
        plt.xlabel('Frequency Lag (MHz)')
        plt.ylabel('Correlation (mK$^2$)')


    def plot_contour(self, filename, norms=False, lag_inds=(0),
                     cross_power=False, title=None, coloraxis=[]):
        lag_inds = list(lag_inds)
        # Set up binning.
        nf = len(self.freq_inds)
        n_bins = 20
        factor = 1.5
        #start = 2.1e6
        freq_diffs = sp.empty(n_bins)
        freq_diffs[0] = 0.0001
        freq_diffs[1] = 2.5 * 200.0 / 256
        freq_diffs[2] = 4.5 * 200.0 / 256
        freq_diffs[3] = 6.5 * 200.0 / 256
        for ii in range(4, n_bins):
            freq_diffs[ii] = factor * freq_diffs[ii - 1]
        freq_diffs *= 1e6
        n_diffs = len(freq_diffs)
        # Allowcate memory.
        corrf = sp.zeros((n_diffs, len(lag_inds)))
        countsf = sp.zeros(n_diffs, dtype=int)
        for ii in range(nf):
            for jj in range(nf):
                if norms:
                    thiscorr = (self.corr[ii, jj, lag_inds] *
                                self.norms[ii, jj, sp.newaxis])
                else:
                    thiscorr = self.corr[ii, jj, lag_inds]
                df = abs(self.freq1[ii] - self.freq2[jj])
                for kk in range(1, n_diffs - 1):
                    if (abs(freq_diffs[kk] - df) <= abs(freq_diffs[kk - 1] - df)
                        and abs(freq_diffs[kk] - df) < abs(freq_diffs[kk + 1] - df)):
                        d_ind = kk
                if abs(freq_diffs[0] - df) < abs(freq_diffs[1] - df):
                    d_ind = 0
                if abs(freq_diffs[-1] - df) <= abs(freq_diffs[-2] - df):
                    d_ind = n_diffs - 1
                corrf[d_ind, :] += thiscorr
                countsf[d_ind] += 1
        corrf /= countsf[:, sp.newaxis]
        if cross_power:
            pdat = corrf * 1e3
        else:
            pdat = sp.sign(corrf) * sp.sqrt(abs(corrf)) * 1e3
        a = plt.figure()
        #a.set_figwidth(a.get_figwidth() / 3.0)
        if len(coloraxis) > 0:
            f = plt.contourf(self.lags[lag_inds], (freq_diffs) / 1e6, pdat, coloraxis)
        else:
            f = plt.contourf(self.lags[lag_inds], (freq_diffs) / 1e6, pdat)

        f.ax.set_xscale('log')
        f.ax.set_yscale('log')

        #im = plt.pcolormesh(self.lags[lag_inds], (freq_diffs) / 1e6, pdat,
        #                    shading='gouraud')
        #im.axes.set_xscale('log')
        #im.axes.set_yscale('log')

        plt.axis('scaled')
        plt.xlim((0.05, 0.9))
        plt.ylim((0.8, 100))
        plt.xlabel("angular lag, $\sigma$ (degrees, 34$\cdotp$Mpc/h)")
        plt.ylabel("frequency lag, $\pi$ (MHz, 4.5$\cdotp$Mpc/h)")
        plt.title(title)
        #c = plt.colorbar(f, ticks=coloraxis)
        c = plt.colorbar(f)

        c.ax.set_ylabel("correlation (mK)")
        plt.savefig(filename)

    def plot_collapsed(self, filename, norms=False, lag_inds=(0), save_old=False,
                       plot_old=False, cross_power=False, title=None):
        lag_inds = list(lag_inds)
        # Set up binning.
        nf = len(self.freq_inds)
        freq_diffs = sp.arange(0.1e6, 100e6, 200.0 / 256 * 1e6)
        n_diffs = len(freq_diffs)
        # Allowcate memory.
        corrf = sp.zeros((n_diffs, len(lag_inds)))
        countsf = sp.zeros(n_diffs, dtype=int)
        for ii in range(nf):
            for jj in range(nf):
                if norms:
                    thiscorr = (self.corr[ii, jj, lag_inds] *
                                self.norms[ii, jj, sp.newaxis])
                else:
                    thiscorr = self.corr[ii, jj, lag_inds]
                df = abs(self.freq1[ii] -self.freq2[jj])
                for kk in range(1, n_diffs - 1):
                    if (abs(freq_diffs[kk] - df) <= abs(freq_diffs[kk - 1] - df)
                        and abs(freq_diffs[kk] - df) < abs(freq_diffs[kk + 1] - df)):
                        d_ind = kk
                if abs(freq_diffs[0] - df) < abs(freq_diffs[1] - df):
                    d_ind = 0
                if abs(freq_diffs[-1] - df) <= abs(freq_diffs[-2] - df):
                    d_ind = n_diffs - 1
                corrf[d_ind, :] += thiscorr
                countsf[d_ind] += 1
        corrf /= countsf[:, sp.newaxis]
        # Now collapse to 1 axis:
        # Cosmology dependant conersion to MPc.
        a_fact = 34.0
        f_fact = 4.5
        nbins = 10
        lags = sp.empty(nbins)
        lags[0] = 2.0
        lags[1] = 4.0
        for ii in range(1, nbins):
            lags[ii] = 1.5 * lags[ii - 1]
        R = self.lags[lag_inds]
        R = (a_fact * R[:, sp.newaxis])**2
        R = R + (f_fact * freq_diffs[sp.newaxis, :] / 1.0e6)**2
        R = sp.sqrt(R)
        corr = sp.zeros(nbins)
        counts = sp.zeros(nbins, dtype=int)
        lag_weights = sp.arange(1, len(lag_inds) + 1)**2
        for ii in range(len(lag_inds)):
            for jj in range(n_diffs):
                dR = R[ii, jj]
                for kk in range(1, nbins - 1):
                    if (abs(lags[kk] - dR) <= abs(lags[kk - 1] - dR)
                        and abs(lags[kk] - dR) < abs(lags[kk + 1] - dR)):
                        ind = kk
                if abs(lags[0] - dR) < abs(lags[1] - dR):
                    ind = 0
                if abs(lags[-1] - dR) <= abs(lags[-2] - dR):
                    ind = nbins - 1
                counts[ind] += lag_weights[ii]
                corr[ind] += lag_weights[ii] * corrf[jj, ii]
        n_boot = 1000
        corrb = ma.zeros((nbins, n_boot))
        countsb = ma.zeros((nbins, n_boot), dtype=int)
        for qq in range(n_boot):
            for mm in range(n_diffs * len(lag_inds)):
                ii = random.random_integers(0, len(lag_inds) - 1)
                jj = random.random_integers(0, n_diffs - 1)
                dR = R[ii, jj]
                for kk in range(1, nbins - 1):
                    if (abs(lags[kk] - dR) <= abs(lags[kk - 1] - dR)
                        and abs(lags[kk] - dR) < abs(lags[kk + 1] - dR)):
                        ind = kk
                if abs(lags[0] - dR) < abs(lags[1] - dR):
                    ind = 0
                if abs(lags[-1] - dR) <= abs(lags[-2] - dR):
                    ind = nbins - 1
                countsb[ind, qq] += lag_weights[ii]
                corrb[ind, qq] += lag_weights[ii] * corrf[jj, ii]
        corr = corr / counts
        corrb = corrb / countsb
        if cross_power:
            pdat = corr * 1e3
            pdatb = corrb * 1e3
        else:
            pdat = sp.sign(corr) * sp.sqrt(abs(corr)) * 1e3
            pdatb = sp.sign(corrb) * sp.sqrt(abs(corrb)) * 1e3

        pdatb_mean = sp.mean(pdatb, -1)
        pdatb_sig = sp.std(pdatb, -1)
        a = plt.figure()
        ax = plt.gca()
        ax.set_yscale("log")
        ax.set_xscale("log")
        errors = sp.empty((2, nbins))
        errors[:, :] = pdatb_sig
        if save_old:
            self.old_pdat = pdat
            self.old_errors = errors
        elin = 2
        msize = 6
        inds = pdat - errors[0, :] > 0
        if sp.any(inds):
            f = plt.errorbar(lags[inds], pdat[inds],
                             errors[:,inds], linestyle='None', marker='o',
                             color='b', elinewidth=elin, markersize=msize)
        inds = pdat + errors[1, :] < 0
        if sp.any(inds):
            f = plt.errorbar(lags[inds], -pdat[inds], errors[:, inds],
                             linestyle='None', marker='o', color='r',
                             elinewidth=elin, markersize=msize)
        inds = sp.logical_and(pdat - errors[0, :] <= 0, pdat > 0)
        if sp.any(inds):
            vals = pdat[inds] + 2 * errors[1, inds]
            es = sp.zeros((2, len(vals)))
            es[0, :] = 0.25 * abs(vals)
            f = plt.errorbar(lags[inds], vals,
                             es, linestyle='None', marker='None',
                             color='b', lolims=True, elinewidth=elin,
                             markersize=msize)
        inds = sp.logical_and(pdat + errors[1, :] >= 0, pdat < 0)
        if sp.any(inds):
            vals = pdat[inds] - 2 * errors[0, inds]
            es = sp.zeros((2, len(vals)))
            es[0, :] = 0.25 * abs(vals)
            f = plt.errorbar(lags[inds], -vals,
                             es, linestyle='None', marker='None',
                             color='r', lolims=True, elinewidth=elin,
                             markersize=msize)
        t_lags = sp.arange(0.1, 100, 0.1)
        r0 = 5.5
        rb = 7.0
        t = (sp.sqrt(((rb + t_lags) / r0)**(-1.8)))
        t = t * 0.15 / t[0]
        f = plt.plot(t_lags, t, marker='None', color='k', linestyle='-')
        if plot_old:
            elin = 0.4
            msize = 6
            mfc = 'w'
            pdat = self.old_pdat
            errors = self.old_errors
            inds = pdat - errors[0, :] > 0
            if sp.any(inds):
                f = plt.errorbar(lags[inds], pdat[inds],
                                 errors[:, inds], linestyle='None',
                                 marker='o', color='b', elinewidth=elin, mfc=mfc,
                                 markersize=msize)
            inds = pdat + errors[1, :] < 0
            if sp.any(inds):
                f = plt.errorbar(lags[inds], -pdat[inds],
                                 errors[:, inds], linestyle='None',
                                 marker='o', color='r', elinewidth=elin, mfc=mfc,
                                 markersize=msize)
            inds = sp.logical_and(pdat - errors[0, :] <= 0, pdat > 0)
            if sp.any(inds):
                vals = pdat[inds] + 2 * errors[1, inds]
                es = sp.zeros((2, len(vals)))
                es[0, :] = 0.25 * abs(vals)
                f = plt.errorbar(lags[inds], vals,
                                 es, linestyle='None', marker='None',
                                 color='b', lolims=True, elinewidth=elin,
                                 markersize=msize, mfc=mfc)
            inds = sp.logical_and(pdat + errors[1, :] >= 0, pdat < 0)
            if sp.any(inds):
                vals = pdat[inds] - 2 * errors[0, inds]
                es = sp.zeros((2, len(vals)))
                es[0, :] = 0.25 * abs(vals)
                f = plt.errorbar(lags[inds], -vals,
                                 es, linestyle='None', marker='None',
                                 color='r', lolims=True, elinewidth=elin,
                                 markersize=msize, mfc=mfc)
        if cross_power:
            plt.axis([1.5, 100, 0.0001, 10.])
        else:
            plt.axis([1.5, 100, 0.01, 500.0])
        plt.xlabel('lag (Mpc/h)')
        plt.ylabel('correlation (mK)')
        plt.title(title)
        plt.savefig(filename)
