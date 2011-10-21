"""Make simulated signal realizations in the GBT survey volume"""
import scipy as sp
from simulations import corr21cm
#from simulations import foregroundsck
#from simulations import pointsource
import core.algebra as algebra
import map.beam as beam
from core import constants as cc
import multiprocessing
from utils import data_paths
import sys


# TODO: confirm ra=x, dec=y (thetax = 5, thetay = 3 in 15hr)
def realize_simulation(template_map, streaming=True):
    """do basic handling to call Richard's simulation code
    here we use 300 h km/s from WiggleZ for streaming dispersion
    Notes on foreground calls for the future (also copy param member func.):
    syn = foregroundsck.Synchrotron()
    synfield = syn.getfield()
    ps = pointsource.DiMatteo()
    psf = ps.getfield()
    """
    if streaming:
        simobj = corr21cm.Corr21cm.like_kiyo_map(template_map, sigma_v=300.*0.72)
    else:
        simobj = corr21cm.Corr21cm.like_kiyo_map(template_map)

    return simobj.get_kiyo_field()


def make_simulation_set(template_file, outfile_raw, outfile_beam,
                        outfile_beam_plus_data, verbose=True, streaming=True):
    """Produce simulated GBT data volumes of three types:
    (from dimensions of a given template file)
    1. the raw simulation
    2. the simulation convolved by the instrumental beam
    3. the simulation convolved by the instrumental beam plus the template
    """
    #(freq_axis, ra_axis, dec_axis, gbt_map_shape, gbt_map) = \
    #    template_map_axes(template_file)
    #gbtsim = realize_simulation(freq_axis, ra_axis, dec_axis,
    #                            streaming=streaming, verbose=verbose)
    template_map = algebra.make_vect(algebra.load(template_file))
    gbtsim = realize_simulation(template_map, streaming=streaming)

    sim_map = algebra.make_vect(gbtsim, axis_names=('freq', 'ra', 'dec'))
    sim_map.copy_axis_info(template_map)
    algebra.save(outfile_raw, sim_map)

    beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031])
    freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
    freq_data *= 1.0e6

    beamobj = beam.GaussianBeam(beam_data, freq_data)
    sim_map_withbeam = beamobj.apply(sim_map)
    algebra.save(outfile_beam, sim_map_withbeam)

    sim_map_withbeam += template_map
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
def generate_sim(template_key, output_key, output_beam_key,
                 output_beam_plus_data_key, streaming=True, parallel=True):
    """generate simulations
    here, assuming the sec A of the set is the template map
    """

    datapath_db = data_paths.DataPath()
    template_mapname = datapath_db.fetch(template_key, pick='A_clean_map',
                                         purpose="template for sim. output",
                                         intend_read=True)

    rawlist = datapath_db.fetch(output_key, intend_write=True,
                                purpose="output sim+", silent=True)

    beamlist = datapath_db.fetch(output_beam_key, intend_write=True,
                                purpose="output sim+beam", silent=True)

    bpdlist = datapath_db.fetch(output_beam_plus_data_key, intend_write=True,
                                purpose="output sim+beam+data", silent=True)

    runlist = [(template_mapname, rawlist[1][index], beamlist[1][index],
                bpdlist[1][index], streaming) for index in rawlist[0]]

    #sys.exit()
    if parallel:
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-3)
        pool.map(wrap_sim, runlist)
    else:
        for runitem in runlist:
            wrap_sim(runitem)


if __name__ == '__main__':

    # alternative template: 'GBT_15hr_Liviu_15mode'

    # 15 hr
    generate_sim('GBT_15hr_map', 'sim_15hr',
                 'sim_15hr_beam', 'sim_15hr_beam_plus_data',
                  streaming=False, parallel=True)

    # 15 hr with velocity streaming
    generate_sim('GBT_15hr_map', 'simvel_15hr',
                 'simvel_15hr_beam', 'simvel_15hr_beam_plus_data',
                  streaming=True, parallel=True)

    # 22 hr
    generate_sim('GBT_22hr_map', 'sim_22hr',
                 'sim_22hr_beam', 'sim_22hr_beam_plus_data',
                  streaming=False, parallel=True)

    # 22 hr with velocity streaming
    generate_sim('GBT_22hr_map', 'simvel_22hr',
                 'simvel_22hr_beam', 'simvel_22hr_beam_plus_data',
                  streaming=True, parallel=True)

    # 1 hr
    generate_sim('GBT_1hr_map', 'sim_1hr',
                 'sim_1hr_beam', 'sim_1hr_beam_plus_data',
                  streaming=False, parallel=True)

    # 1 hr with velocity streaming
    generate_sim('GBT_1hr_map', 'simvel_1hr',
                 'simvel_1hr_beam', 'simvel_1hr_beam_plus_data',
                  streaming=True, parallel=True)
