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
from numpy import random
import struct


# TODO: confirm ra=x, dec=y (thetax = 5, thetay = 3 in 15hr)
def realize_simulation(template_map, scenario=None, seed=None, refinement=2):
    """do basic handling to call Richard's simulation code
    here we use 300 h km/s from WiggleZ for streaming dispersion
    Notes on foreground calls for the future (also copy param member func.):
    syn = foregroundsck.Synchrotron()
    synfield = syn.getfield()
    ps = pointsource.DiMatteo()
    psf = ps.getfield()
    """
    if seed is not None:
        random.seed(seed)

    if scenario is None:
        print "running standard simulation"
        simobj = corr21cm.Corr21cm.like_kiyo_map(template_map)
        maps = simobj.get_kiyo_field_physical(refinement=refinement)

    else:
        if scenario == "streaming":
            print "running streaming simulation"
            simobj = corr21cm.Corr21cm.like_kiyo_map(template_map,
                                                     sigma_v=300.*0.72)
            maps = simobj.get_kiyo_field_physical(refinement=refinement)

        if scenario == "ideal":
            print "running ideal simulation"
            simobj = corr21cm.Corr21cm.like_kiyo_map(template_map)
            maps = simobj.get_kiyo_field_physical(refinement=refinement,
                                density_only=True,
                                no_mean=True,
                                no_evolution=True)

    return maps


def make_simulation_set(template_file, outfile_physical,
                        outfile_raw, outfile_beam, outfile_beam_plus_data,
                        verbose=True, scenario=None):
    """Produce simulated GBT data volumes of three types:
    (from dimensions of a given template file)
    0. the simulation in a physical volume
    1. the raw simulation in z
    2. the simulation convolved by the instrumental beam
    3. the simulation convolved by the instrumental beam plus the template
    """
    template_map = algebra.make_vect(algebra.load(template_file))

    # The usual seed is not fine enough for parallel jobs
    randsource = open("/dev/random", "rb")
    seed = struct.unpack("I", randsource.read(4))[0]
    #seed = abs(long(outfile_physical.__hash__()))
    (gbtsim, gbtphys, physdim) = realize_simulation(template_map,
                                                    scenario=scenario,
                                                    seed=seed)

    phys_map = algebra.make_vect(gbtphys, axis_names=('freq', 'ra', 'dec'))
    pshp = phys_map.shape

    # define the axes of the physical map; several alternatives are commented
    info = {}
    info['axes'] = ('freq', 'ra', 'dec')
    info['type'] = 'vect'
    info['freq_delta'] = abs(physdim[0] - physdim[1]) / float(pshp[0] - 1)
    info['freq_centre'] = physdim[0] + info['freq_delta'] * float(pshp[0] // 2)
    #        'freq_centre': abs(physdim[0] + physdim[1]) / 2.,

    info['ra_delta'] = abs(physdim[2]) / float(pshp[1] - 1)
    #info['ra_centre'] = info['ra_delta'] * float(pshp[1] // 2)
    #        'ra_centre': abs(physdim[2]) / 2.,
    info['ra_centre'] = 0.

    info['dec_delta'] = abs(physdim[3]) / float(pshp[2] - 1)
    #info['dec_centre'] = info['dec_delta'] * float(pshp[2] // 2)
    #        'dec_centre': abs(physdim[3]) / 2.,
    info['dec_centre'] = 0.

    phys_map.info = info

    sim_map = algebra.make_vect(gbtsim, axis_names=('freq', 'ra', 'dec'))
    sim_map.copy_axis_info(template_map)

    algebra.save(outfile_raw, sim_map)
    algebra.save(outfile_physical, phys_map)

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
    (template_mapname, outfile_physical, outfile_raw, outfile_beam, \
                                outfile_beam_plus_data, scenario) = runitem

    print "using template: " + template_mapname
    print "using physical-space cube file: " + outfile_physical
    print "using raw output cube file: " + outfile_raw
    print "using raw output cube file conv. by beam: " + outfile_beam
    print "using raw output cube file conv. by beam plus data: " + outfile_beam_plus_data

    make_simulation_set(template_mapname, outfile_physical, outfile_raw, outfile_beam,
                        outfile_beam_plus_data, scenario=scenario)


def test_scheme(template_file, sim_filename1, sim_filename2):
    r"""look at some differences between maps"""
    template_map = algebra.make_vect(algebra.load(template_file))
    gbtsim1 = realize_simulation(template_map, scenario='streaming',
                                seed=5489, refinement=1.)
    gbtsim2 = realize_simulation(template_map,
                                    seed=5489, refinement=1.)

    sim_map1 = algebra.make_vect(gbtsim1, axis_names=('freq', 'ra', 'dec'))
    sim_map2 = algebra.make_vect(gbtsim2, axis_names=('freq', 'ra', 'dec'))
    sim_map1.copy_axis_info(template_map)
    sim_map2.copy_axis_info(template_map)
    algebra.save(sim_filename1, sim_map1)
    algebra.save(sim_filename2, sim_map2)

def generate_sim(template_key, output_physical_key, output_key,
                 output_beam_key, output_beam_plus_data_key,
                 scenario=None, parallel=True):
    """generate simulations
    here, assuming the sec A of the set is the template map
    """

    datapath_db = data_paths.DataPath()
    template_mapname = datapath_db.fetch(template_key, pick='A;clean_map',
                                         purpose="template for sim. output",
                                         intend_read=True)

    physlist = datapath_db.fetch(output_physical_key, intend_write=True,
                                purpose="output sim+", silent=True)

    rawlist = datapath_db.fetch(output_key, intend_write=True,
                                purpose="output sim+", silent=True)

    beamlist = datapath_db.fetch(output_beam_key, intend_write=True,
                                purpose="output sim+beam", silent=True)

    bpdlist = datapath_db.fetch(output_beam_plus_data_key, intend_write=True,
                                purpose="output sim+beam+data", silent=True)

    runlist = [(template_mapname, physlist[1][index], rawlist[1][index],
                beamlist[1][index], bpdlist[1][index], scenario)
                for index in rawlist[0]]

    #sys.exit()
    if parallel:
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-4)
        pool.map(wrap_sim, runlist)
    else:
        for runitem in runlist:
            wrap_sim(runitem)


def generate_simset(fieldname, scenario=None):
    template_key = 'GBT_%s_map' % fieldname

    if scenario is None:
        tag= "sim"
    else:
        if scenario == "streaming":
            tag = "simvel"

        if scenario == "ideal":
            tag = "simideal"

    output_physical_key = '%s_%s_physical' % (tag, fieldname)
    output_key = '%s_%s' % (tag, fieldname)
    output_beam_key = '%s_%s_beam' % (tag, fieldname)
    output_beam_plus_data_key = '%s_%s_beam_plus_data' % (tag, fieldname)

    generate_sim(template_key, output_physical_key, output_key,
                 output_beam_key, output_beam_plus_data_key,
                 scenario=scenario, parallel=True)


def generate_full_simset(fieldlist):
    for fieldname in fieldlist:
        generate_simset(fieldname, scenario="ideal")
        generate_simset(fieldname)
        generate_simset(fieldname, scenario="streaming")


def run_scheme_test():
    template_file = "/mnt/raid-project/gmrt/tcv/maps/sec_A_15hr_41-90_clean_map_I.npy"
    sim_filename1 = "sim_streaming1.npy"
    sim_filename2 = "sim_streaming2.npy"
    test_scheme(template_file, sim_filename1, sim_filename2)


if __name__ == '__main__':
    generate_full_simset(['15hr', '22hr', '1hr'])
    #generate_full_simset(['15hr', '22hr'])
    #generate_full_simset(['15hr'])
    #run_scheme_test()
