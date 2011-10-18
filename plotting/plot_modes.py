import matplotlib.pyplot as plt
import scipy as sp

def plot_svd(vals):
    r"""Deprecated.

    Plots the svd values and prints out some statistics."""

    n_vals = len(vals)
    plt.semilogy(abs(sp.sort(-vals / n_vals)), marker='o',
                 linestyle='None')
    print 'Mean noise: ', sp.sum(vals) / n_vals
    print 'Largest eigenvalues/vals: ',
    print sp.sort(vals / n_vals)[-10:]



