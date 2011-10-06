"""Make simulated signal realizations in the GBT survey volume"""
import scipy as sp
from simulations import corr21cm
#from simulations import foregroundsck
#from simulations import pointsource
import core.algebra as algebra
import map.beam as beam
from core import constants as cc


# TODO: TEMPORARY duplicate of func in optical_catalog; move more generic
def template_map_axes(filename):
    """Open a numpy array map and extract its axis/etc. information
    """
    print "using the volume template file: " + filename
    template_map = algebra.make_vect(algebra.load(filename))
    freq_axis = template_map.get_axis('freq')
    ra_axis = template_map.get_axis('ra')
    dec_axis = template_map.get_axis('dec')
    return (freq_axis, ra_axis, dec_axis, template_map.shape, template_map)


# TODO: confirm ra=x, dec=y (thetax = 5, thetay = 3 in 15hr)
def realize_simulation(freq_axis, ra_axis, dec_axis, verbose=True):
    """do basic handling to call Richard's simulation code
    Notes on foreground calls for the future (also copy param member func.):
    syn = foregroundsck.Synchrotron()
    synfield = syn.getfield()
    ps = pointsource.DiMatteo()
    psf = ps.getfield()
    """
    simobj = corr21cm.Corr21cm()
    simobj.x_width = max(ra_axis) - min(ra_axis)
    simobj.y_width = max(dec_axis) - min(dec_axis)
    (simobj.x_num, simobj.y_num) = (len(ra_axis), len(dec_axis))

    simobj.nu_lower = min(freq_axis)
    simobj.nu_upper = max(freq_axis)
    simobj.nu_num = len(freq_axis)

    if verbose:
        print "Sim: %dx%d field (%fx%f deg) from nu=%f to nu=%f (%d bins)" % \
              (simobj.x_num, simobj.y_num, simobj.x_width, simobj.y_width,
               simobj.nu_lower, simobj.nu_upper, simobj.nu_num)

    # note that the temperature there is in mK, and one wants K locally
    return simobj.getfield() * 0.001


def realize_simulation_old(freq_axis, ra_axis, dec_axis, verbose=True):
    """do basic handling to call Richard's simulation code"""
    thetax = max(ra_axis) - min(ra_axis)
    thetay = max(dec_axis) - min(dec_axis)
    (num_ra, num_dec) = (len(ra_axis), len(dec_axis))

    freq_near = max(freq_axis)
    freq_far = min(freq_axis)
    num_freq = len(freq_axis)
    z_near = cc.freq_21cm_MHz / freq_near
    z_far = cc.freq_21cm_MHz / freq_far

    simobj = corr21cm.Corr21cm()

    if verbose:
        print "Sim: %dx%d field (%fx%f deg) from z=%f to z=%f (%d bins)" % \
              (num_ra, num_dec, thetax, thetay, z_near, z_far, num_freq)

    # note that the temperature there is in mK, and one wants K locally
    return simobj.realisation(z_near, z_far, thetax, thetay,
                          num_freq, num_ra, num_dec)[::-1, ...] * 0.001


def make_simulation_set(template_file, outfile_raw, outfile_beam,
                        outfile_beam_plus_data, verbose=True):
    """Produce simulated GBT data volumes of three types:
    (from dimensions of a given template file)
    1. the raw simulation
    2. the simulation convolved by the instrumental beam
    3. the simulation convolved by the instrumental beam plus the template
    """
    (freq_axis, ra_axis, dec_axis, gbt_map_shape, gbt_map) = \
        template_map_axes(template_file)

    gbtsim = realize_simulation(freq_axis, ra_axis, dec_axis, verbose=verbose)
    sim_map = algebra.make_vect(gbtsim, axis_names=('freq', 'ra', 'dec'))

    # TODO: need this?
    (num_ra, num_dec) = (len(ra_axis), len(dec_axis))
    sim_map = sim_map[:, :num_ra, :num_dec]

    sim_map.copy_axis_info(gbt_map)

    beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031])
    freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
    freq_data *= 1.0e6
    beamobj = beam.GaussianBeam(beam_data, freq_data)
    sim_map_withbeam = beamobj.apply(sim_map)

    algebra.save(outfile_raw, sim_map)
    algebra.save(outfile_beam, sim_map_withbeam)

    sim_map_withbeam += gbt_map

    algebra.save(outfile_beam_plus_data, sim_map_withbeam)


def generate_sim_15hr_gbtregion():
    """generate simulations
    """

    calinliv_dir = '/mnt/raid-project/gmrt/calinliv/wiggleZ/'
    root_template = calinliv_dir + "corr/test/"
    template_mapname = root_template + \
                    'sec_A_15hr_41-73_cleaned_clean_map_I_with_B.npy'

    simdir = '/mnt/raid-project/gmrt/eswitzer/wiggleZ/simulations/15hr/'
    outfile_raw = simdir + 'sim_test.npy'
    outfile_beam = simdir + 'sim_beam_test.npy'
    outfile_beam_plus_data = simdir + 'sim_beam_plus_data_test.npy'

    make_simulation_set(template_mapname, outfile_raw, outfile_beam,
                        outfile_beam_plus_data)


if __name__ == '__main__':
    generate_sim_15hr_gbtregion()
