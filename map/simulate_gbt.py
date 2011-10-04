"""Make simulated signal realizations in the GBT survey volume"""
import numpy as np
import scipy as sp
from simulations import corr21cm
import core.algebra as algebra
import map.beam as beam
from core import constants as cc


# TODO: interpolate const redshift from Corr21cm.realisation onto const freq
# TODO: confirm ra=x, dec=y (thetax = 5, thetay = 3 in 15hr)
def realize_simulation(freq_axis, ra_axis, dec_axis, verbose=True):
    """do basic handling to call Richard's simulation code"""
    thetax = max(ra_axis) - min(ra_axis)
    thetay = max(dec_axis) - min(dec_axis)
    (num_ra, num_dec) = (len(ra_axis), len(dec_axis)

    freq_near = max(freq_axis)
    freq_far = min(freq_axis)
    num_freq = len(freq_axis)
    z_near = cc.freq_21cm_MHz / freq_near
    z_far = cc.freq_21cm_MHz / freq_far

    cr = corr21cm.Corr21cm()

    if verbose:
        print "Sim: %dx%d field (%fx%f degrees) from z=%f to z=%f (%d bins)" %
              (num_ra, num_dec, thetax, thetay, z_near, z_far, num_freq)

    # note that the temperature there is in mK, and one wants K locally
    return cr.realisation(z_near, z_far, thetax, thetay,
                          num_freq, num_ra, num_dec)[::-1, ...] * 0.001


def make_simulation_set(template_file, outfile_raw, outfile_beam,
                        outfile_beam_plus_data):
    """Produce simulated GBT data volumes of three types:
    (from dimensions of a given template file)
    1. the raw simulation
    2. the simulation convolved by the instrumental beam
    3. the simulation convolved by the instrumental beam plus the template
    """
    exMap = algebra.make_vect(algebra.load(template_file))
    # TODO: get axes

    rf = realize_simulation(freq_axis, ra_axis, dec_axis, verbose=verbose)
    a = algebra.make_vect(rf, axis_names=('freq', 'ra', 'dec'))
    a = a[:, :nx, :ny]
    # TODO: copy axis
    a.info = info

    beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031])
    freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
    freq_data *= 1.0e6
    #b = beam.GaussianBeam(width = [0.25, 0.25*freq_far/freq_near], freq = [freq_far, freq_near])
    b = beam.GaussianBeam(beam_data, freq_data)
    ab = b.apply(a)

    algebra.save(outfile_raw, a)
    algebra.save(outfile_beam, ab)

    ab += exMap

    algebra.save(outfile_beam_plus_data, ab)



