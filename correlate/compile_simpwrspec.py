import numpy as np
import math
from core import algebra
from utils import data_paths
from simulations import corr21cm
import multiprocessing
from map import physical_gridding
from utils import binning
from correlate import map_pair as mp
from correlate import pwrspec_estimation as pe
from correlate import transfer_function as tf
import sys
from utils import batch_handler
import copy


def gather_batch_sim_run(simleft_key, simright_key,
                         weightleft_key, weightright_key,
                         inifile=None, datapath_db=None):

    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "reading from: " + batchquad_path

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=False, verbose=True)

    pwr_1d = []
    pwr_1d_from_2d = []
    pwr_2d = []
    mock_cases = datapath_db.fileset_cases(simleft_key, "realization")

    for index in mock_cases['realization']:
        input = {}
        input['map1_key'] = "%s:%s" % (simleft_key, index)
        input['map2_key'] = "%s:%s" % (simright_key, index)
        input['noiseinv1_key'] = "%s" % weightleft_key
        input['noiseinv2_key'] = "%s" % weightright_key
        files = convert_keydict_to_filedict(input, db=datapath_db)

        pwr2d_run, pwr1d_run = caller.execute(files['map1_key'],
                                              files['map2_key'],
                                              files['noiseinv1_key'],
                                              files['noiseinv2_key'],
                                              inifile=inifile)

        pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwr2d_run,
                                                  transfer=transfer))

        pwr_2d.append(pwr2d_run)
        pwr_1d.append(pwr1d_run)

    pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d, tag,
                         outdir="./plot_data_v2")

    return (pwr_1d, pwr_1d_from_2d, pwr_2d)


def gather_physical_sim_run(sim_key, inifile=None, datapath_db=None):

    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "reading from: " + batchquad_path

    funcname = "correlate.batch_quadratic.call_phys_space_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=False, verbose=True)

    mock_cases = datapath_db.fileset_cases(sim_key, "realization")

    pwr_1d = []
    pwr_1d_from_2d = []
    pwr_2d = []
    for index in mock_cases['realization']:
        map1_file = datapath_db.fetch("%s:%s" % (sim_key, index))
        map2_file = datapath_db.fetch("%s:%s" % (sim_key, index))

        pwr2d_run, pwr1d_run = caller.execute(map1_file, map2_file,
                                              inifile=inifile)

        pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwr2d_run,
                                                  transfer=transfer))

        pwr_2d.append(pwr2d_run)
        pwr_1d.append(pwr1d_run)

    pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d, tag,
                         outdir="./plot_data_v2")

    # now report the theory input curve over the same bins
    bin_left = pwr_1d[0]['bin_left']
    bin_center = pwr_1d[0]['bin_center']
    bin_right = pwr_1d[0]['bin_right']

    zfilename = datapath_db.fetch(redshift_filekey, intend_read=True,
                                 pick='1')

    fileout = outdir + "/" + "theory_power_spectrum.dat"

    zspace_cube = algebra.make_vect(algebra.load(zfilename))
    simobj = corr21cm.Corr21cm.like_kiyo_map(zspace_cube)
    pwrspec_input = simobj.get_pwrspec(bin_center)
    if unitless:
        pwrspec_input *= bin_center ** 3. / 2. / math.pi / math.pi

    outfile = open(fileout, "w")
    for specdata in zip(bin_left, bin_center, bin_right, pwrspec_input):
        outfile.write(("%10.15g " * 4 + "\n") % specdata)

    outfile.close()

    return (pwr_1d, pwr_1d_from_2d, pwr_2d)


def plot_all_simpwr(fileroot, basics_only=False):
    r"""currently fileroot is sim, simvel or simideal"""

    gather_physical_sim_run("%s_physical" % fileroot,
                            "%s_physical" % fileroot)

    gather_batch_sim_run(fileroot,
                         fileroot)

    if not basics_only:
        gather_batch_sim_run("%s_beam" % fileroot,
                             "%s_beam" % fileroot)

        gather_batch_sim_run("%s_beam" % fileroot,
                             "%s_beam_meanconv" % fileroot,
                             degrade_resolution=True, subtract_mean=True)

        gather_batch_sim_run("%s_beam" % fileroot,
                             "%s_beam_mean" % fileroot,
                             degrade_resolution=False, subtract_mean=True)

        gather_batch_sim_run("%s_beam" % fileroot,
                             "%s_beam_conv" % fileroot,
                             degrade_resolution=True, subtract_mean=False)


if __name__ == '__main__':
    plot_all_simpwr("sim_15hr")
    plot_all_simpwr("simvel_15hr")
    plot_all_simpwr("simideal_15hr")
    sys.exit()


