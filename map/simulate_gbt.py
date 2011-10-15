"""Make simulated signal realizations in the GBT survey volume"""
import scipy as sp
from simulations import corr21cm
#from simulations import foregroundsck
#from simulations import pointsource
import core.algebra as algebra
import map.beam as beam
from core import constants as cc
import multiprocessing


# TODO: TEMPORARY duplicate of func in optical_catalog; move more generic
def template_map_axes(filename):
    """Open a numpy array map and extract its axis/etc. information
    """
    template_map = algebra.make_vect(algebra.load(filename))
    freq_axis = template_map.get_axis('freq')
    ra_axis = template_map.get_axis('ra')
    dec_axis = template_map.get_axis('dec')
    return (freq_axis, ra_axis, dec_axis, template_map.shape, template_map)


# TODO: confirm ra=x, dec=y (thetax = 5, thetay = 3 in 15hr)
def realize_simulation(freq_axis, ra_axis, dec_axis, verbose=True,
                       streaming=True):
    """do basic handling to call Richard's simulation code
    here we use 300 h km/s from WiggleZ for streaming dispersion
    Notes on foreground calls for the future (also copy param member func.):
    syn = foregroundsck.Synchrotron()
    synfield = syn.getfield()
    ps = pointsource.DiMatteo()
    psf = ps.getfield()
    """
    if streaming:
        simobj = corr21cm.Corr21cm(sigma_v=300.*0.72)
    else:
        simobj = corr21cm.Corr21cm()

    simobj.x_width = max(ra_axis) - min(ra_axis)
    simobj.y_width = max(dec_axis) - min(dec_axis)
    (simobj.x_num, simobj.y_num) = (len(ra_axis), len(dec_axis))

    simobj.nu_lower = min(freq_axis)/1.e6
    simobj.nu_upper = max(freq_axis)/1.e6
    simobj.nu_num = len(freq_axis)

    if verbose:
        print "Sim: %dx%d field (%fx%f deg) from nu=%f to nu=%f (%d bins)" % \
              (simobj.x_num, simobj.y_num, simobj.x_width, simobj.y_width,
               simobj.nu_lower, simobj.nu_upper, simobj.nu_num)

    # note that the temperature there is in mK, and one wants K locally
    return simobj.getfield() * 0.001


def make_simulation_set(template_file, outfile_raw, outfile_beam,
                        outfile_beam_plus_data, verbose=True, streaming=True):
    """Produce simulated GBT data volumes of three types:
    (from dimensions of a given template file)
    1. the raw simulation
    2. the simulation convolved by the instrumental beam
    3. the simulation convolved by the instrumental beam plus the template
    """
    (freq_axis, ra_axis, dec_axis, gbt_map_shape, gbt_map) = \
        template_map_axes(template_file)

    gbtsim = realize_simulation(freq_axis, ra_axis, dec_axis,
                                streaming=streaming, verbose=verbose)
    sim_map = algebra.make_vect(gbtsim, axis_names=('freq', 'ra', 'dec'))
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


def wrap_sim(runitem):
    (template_mapname, outfile_raw, outfile_beam, outfile_beam_plus_data, streaming) = runitem
    print "using template: " + template_mapname
    print "using raw output cube file: " + outfile_raw
    print "using raw output cube file conv. by beam: " + outfile_beam
    print "using raw output cube file conv. by beam plus data: " + outfile_beam_plus_data

    make_simulation_set(template_mapname, outfile_raw, outfile_beam,
                        outfile_beam_plus_data, streaming=streaming)


# TODO: have this use the path data base
def generate_sim_15hr_gbtregion(numsim=100, streaming=True):
    """generate simulations
    """

    calinliv_dir = '/mnt/raid-project/gmrt/calinliv/wiggleZ/'
    root_template = calinliv_dir + "corr/test/"
    template_mapname = root_template + \
                    'sec_A_15hr_41-73_cleaned_clean_map_I_with_B.npy'

    simdir = '/mnt/raid-project/gmrt/eswitzer/wiggleZ/simulations/15hr/'
    if streaming:
        simdir += "simvel_"
    else:
        simdir += "sim_"

    templatelist = [ template_mapname for index in range(numsim) ]

    rawlist = [ simdir + ('%03d' % index) + '.npy' \
                for index in range(numsim) ]

    beamlist = [ simdir + ('beam_%03d' % index) + '.npy' \
                for index in range(numsim) ]

    bpdlist = [ simdir + ('beam_plus_data_%03d' % index) + '.npy' \
                for index in range(numsim) ]

    streamlist = [ streaming for index in range(numsim) ]

    runlist = zip(templatelist, rawlist, beamlist, bpdlist, streamlist)

    for runitem in runlist:
        wrap_sim(runitem)

    # the version of python which supports the simulation code does not support
    # multiprocessing...
    #pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-3)
    #pool.map(wrap_sim, runlist)


if __name__ == '__main__':
    generate_sim_15hr_gbtregion(streaming=False)
