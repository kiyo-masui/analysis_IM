import scipy as sp
import numpy as np
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
from kiyopy import parse_ini
from utils import file_tools
from utils import units
from map import simulate_gbt as sg

def wrap_sim(runitem):
    (template_mapname, outfile_physical, outfile_raw, outfile_beam, \
                                outfile_beam_plus_data, scenario) = runitem

    print "using template: " + template_mapname
    print "using physical-space cube file: " + outfile_physical
    print "using raw output cube file: " + outfile_raw
    print "using raw output cube file conv. by beam: " + outfile_beam
    print "using raw output cube file conv. by beam plus data: " + outfile_beam_plus_data

    sg.make_simulation_set(template_mapname, outfile_physical, outfile_raw, outfile_beam,
                        outfile_beam_plus_data, scenario=scenario)


def generate_sim(params, parallel=True, silent=True, datapath_db=None):
    """generate simulations
    here, assuming the sec A of the set is the template map
    """

    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    template_mapname = datapath_db.fetch(params['template_key'],
                                         intend_read=True,
                                         purpose="template_map",
                                         silent=silent)

    physlist = datapath_db.fetch(params['sim_physical_key'],
                                 intend_write=True,
                                 purpose="output sim (physical)",
                                 silent=silent)

    rawlist = datapath_db.fetch(params['sim_key'],
                                intend_write=True,
                                purpose="output sim",
                                silent=silent)

    beamlist = datapath_db.fetch(params['sim_beam_key'],
                                 intend_write=True,
                                 purpose="output sim+beam",
                                 silent=silent)

    bpdlist = datapath_db.fetch(params['sim_beam_plus_data_key'],
                                intend_write=True,
                                purpose="output sim+beam+data",
                                silent=silent)

    runlist = [(template_mapname, physlist[1][index], rawlist[1][index],
                beamlist[1][index], bpdlist[1][index],
                params['pwrspec_scenario'])
                for index in rawlist[0]]

    print runlist
    #sys.exit()
    if parallel:
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-4)
        pool.map(wrap_sim, runlist)
    else:
        for runitem in runlist:
            wrap_sim(runitem)


def generate_aux_simset(params, silent=False, datapath_db=None):

    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    weightfile = datapath_db.fetch(params['weight_key'],
                                         intend_read=True,
                                         purpose="weight map",
                                         silent=silent)

    input_rawsimset = datapath_db.fetch(params['sim_key'],
                                        intend_read=True, silent=silent)

    output_deltasimset = datapath_db.fetch(params['sim_delta_key'],
                                           intend_write=True, silent=silent)

    input_beamsimset = datapath_db.fetch(params['sim_beam_key'],
                                         intend_read=True, silent=silent)

    output_meansubsimset = datapath_db.fetch(params['sim_beam_meansub_key'],
                                         intend_write=True, silent=silent)

    output_convsimset = datapath_db.fetch(params['sim_beam_conv_key'],
                                         intend_write=True, silent=silent)

    output_meansubconvsimset = datapath_db.fetch(
                                         params['sim_beam_meansubconv_key'],
                                         intend_write=True, silent=silent)

    # TODO: actually implement the foreground simulations
    output_fgsimset = datapath_db.fetch(params['sim_beam_plus_fg_key'],
                                         intend_write=True, silent=silent)

    for index in input_rawsimset[0]:
        sg.generate_delta_sim(input_rawsimset[1][index],
                           output_deltasimset[1][index])

        sg.generate_proc_sim(input_beamsimset[1][index], weightfile,
                          output_meansubsimset[1][index],
                          meansub=True, degrade=False)

        sg.generate_proc_sim(input_beamsimset[1][index], weightfile,
                          output_convsimset[1][index],
                          meansub=False, degrade=True)

        sg.generate_proc_sim(input_beamsimset[1][index], weightfile,
                          output_meansubconvsimset[1][index],
                          meansub=True, degrade=True)


# cases: [15hr, 22hr, 1hr], [ideal, nostr, str]
params_init = {
               'output_root': '',
               'sim_physical_key': '',
               'sim_key': '',
               'sim_beam_key': '',
               'sim_beam_plus_fg_key': '',
               'sim_beam_plus_data_key': '',
               'sim_delta_key': '',
               'sim_beam_meansub_key': '',
               'sim_beam_conv_key': '',
               'sim_beam_meansubconv_key': '',
               'template_key': '',
               'weight_key': '',
               'pwrspec_scenario': '',
               'omega_HI': ''
               }
prefix = 'sg_'

if __name__ == '__main__':
    params = parse_ini.parse(str(sys.argv[1]), params_init,
                             prefix=prefix)
    print params

    datapath_db = data_paths.DataPath()
    output_root = datapath_db.fetch(params['output_root'])

    file_tools.mkparents(output_root)
    parse_ini.write_params(params, output_root + 'params.ini',
                           prefix=prefix)

    generate_sim(params, parallel=True, datapath_db=datapath_db)
    generate_aux_simset(params, datapath_db=datapath_db)

