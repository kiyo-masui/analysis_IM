r"""Several scripts to call plot_cube_movie and make movies of GBT data,
simulations, etc.
"""
from plotting import plot_cube as pc
from utils import data_paths
# TODO: reimplement transverse plotting

def unique_list(listin):
    seen = set()
    seen_add = seen.add
    uniq = [ x for x in listin if x not in seen and not seen_add(x)]
    return sorted(uniq)

def plot_gbt_mapset(outputdir="/cita/d/www/home/eswitzer/movies/",
                    make_map=True, make_dirty_map=False,
                    make_noise_inv=False, make_noise_diag=True):

    pc.plot_gbt_maps('GBT_15hr_map_oldcalmoon', outputdir=outputdir,
                     make_map=make_map, make_dirty_map=make_dirty_map,
                     make_noise_inv=make_noise_inv,
                     make_noise_diag=make_noise_diag)

    #pc.plot_gbt_maps('GBT_15hr_map_oldcal', outputdir=outputdir,
    #                 make_map=make_map, make_dirty_map=make_dirty_map,
    #                 make_noise_inv=make_noise_inv,
    #                 make_noise_diag=make_noise_diag)

    #pc.plot_gbt_maps('GBT_15hr_optimalmap_fluxpolcal', outputdir=outputdir,
    #                 make_map=make_map, make_dirty_map=make_dirty_map,
    #                 make_noise_inv=make_noise_inv,
    #                 make_noise_diag=make_noise_diag)

    #pc.plot_gbt_maps('GBT_15hr_optimalmap_mapv2oldcal', outputdir=outputdir,
    #                 make_map=make_map, make_dirty_map=make_dirty_map,
    #                 make_noise_inv=make_noise_inv,
    #                 make_noise_diag=make_noise_diag)

    #pc.plot_gbt_maps('GBT_15hr_optimalmap_mapv2fdgcal', outputdir=outputdir,
    #                 make_map=make_map, make_dirty_map=make_dirty_map,
    #                 make_noise_inv=make_noise_inv,
    #                 make_noise_diag=make_noise_diag)

    #pc.plot_gbt_maps('GBT_15hr_optimalmap_mapv2fdgcalmoderm', outputdir=outputdir,
    #                 make_map=make_map, make_dirty_map=make_dirty_map,
    #                 make_noise_inv=make_noise_inv,
    #                 make_noise_diag=make_noise_diag)


    #pc.plot_gbt_maps('GBT_15hr_map_fdgcal', outputdir=outputdir, transverse=False)
    #pc.plot_gbt_maps('GBT_15hr_oldcal', outputdir=outputdir, transverse=False)
    #pc.plot_gbt_maps('GBT_15hr_map', outputdir=outputdir, transverse=False)
    #pc.plot_gbt_maps('GBT_22hr_map', outputdir=outputdir, transverse=False)
    #pc.plot_gbt_maps('GBT_1hr_map', outputdir=outputdir, transverse=False)
    #pc.plot_gbt_maps('GBT_15hr_map_proposal', transverse=False,
    #                 outputdir=outputdir, skip_noise=True)

def plot_gbt_newmapset(outputdir="/cita/d/www/home/eswitzer/movies/"):
    pc.plot_gbt_maps('GBT_15hr_newmap737', outputdir=outputdir, transverse=False)
    pc.plot_gbt_maps('GBT_15hr_newmap799', outputdir=outputdir, transverse=False)
    pc.plot_gbt_maps('GBT_15hr_newmap862', outputdir=outputdir, transverse=False)

def plot_cleaned_gbt15hr_mapset(outputdir="/cita/d/www/home/eswitzer/movies/"):
    pc.plot_cleaned_maps('GBT_15hr_map_cleaned_noconv_25mode',
                         outputdir=outputdir, transverse=False, convolve=True)
    #pc.plot_cleaned_maps('GBT_15hr_map_cleaned_0mode', outputdir=outputdir, transverse=False)
    #pc.plot_cleaned_maps('GBT_15hr_map_cleaned_5mode', outputdir=outputdir, transverse=False)
    #pc.plot_cleaned_maps('GBT_15hr_map_cleaned_10mode', outputdir=outputdir, transverse=False)
    #pc.plot_cleaned_maps('GBT_15hr_map_cleaned_15mode', outputdir=outputdir, transverse=False)
    #pc.plot_cleaned_maps('GBT_15hr_map_cleaned_20mode', outputdir=outputdir, transverse=False)
    #pc.plot_cleaned_maps('GBT_15hr_map_cleaned_25mode', outputdir=outputdir, transverse=False)
    #pc.plot_cleaned_maps('GBT_15hr_map_cleaned_30mode', outputdir=outputdir, transverse=False)
    #pc.plot_cleaned_maps('GBT_15hr_map_cleaned_35mode', outputdir=outputdir, transverse=False)
    #pc.plot_cleaned_maps('GBT_15hr_map_cleaned_40mode', outputdir=outputdir, transverse=False)
    #pc.plot_cleaned_maps('GBT_15hr_map_cleaned_45mode', outputdir=outputdir, transverse=False)
    #pc.plot_cleaned_maps('GBT_15hr_map_cleaned_50mode', outputdir=outputdir, transverse=False)

def plot_gbt_simset(fieldname, outputdir="/cita/d/www/home/eswitzer/movies/"):
    datapath_db = data_paths.DataPath()

    #keyname = "%s" % fieldname
    #filename = datapath_db.fetch(keyname, pick='0')
    #pc.make_cube_movie(filename, "Temperature (mK)", pc.cube_frame_dir,
    #                    sigmarange=3., outputdir=outputdir, multiplier=1000.,
    #                    transverse=False, filetag_suffix="_"+fieldname)

    #keyname = "%s_beam" % fieldname
    #filename = datapath_db.fetch(keyname, pick='0')
    #pc.make_cube_movie(filename, "Temperature (mK)", pc.cube_frame_dir,
    #                    sigmarange=3., outputdir=outputdir, multiplier=1000.,
    #                    transverse=False, filetag_suffix="_"+fieldname)

    #keyname = "%s_beam_conv" % fieldname
    #filename = datapath_db.fetch(keyname, pick='0')
    #pc.make_cube_movie(filename, "Temperature (mK)", pc.cube_frame_dir,
    #                    sigmarange=3., outputdir=outputdir, multiplier=1000.,
    #                    transverse=False, filetag_suffix="_"+fieldname)

    keyname = "%s_beam_meansub" % fieldname
    filename = datapath_db.fetch(keyname, pick='0')
    pc.make_cube_movie(filename, "Temperature (mK)", pc.cube_frame_dir,
                        sigmarange=3., outputdir=outputdir, multiplier=1000.,
                        transverse=False, filetag_suffix="_"+fieldname)

    keyname = "%s_beam_meansubconv" % fieldname
    filename = datapath_db.fetch(keyname, pick='0')
    pc.make_cube_movie(filename, "Temperature (mK)", pc.cube_frame_dir,
                        sigmarange=3., outputdir=outputdir, multiplier=1000.,
                        transverse=False, filetag_suffix="_"+fieldname)

    #keyname = "%s_delta" % fieldname
    #filename = datapath_db.fetch(keyname, pick='0')
    #pc.make_cube_movie(filename, "overdensity", pc.cube_frame_dir,
    #                    sigmarange=3., outputdir=outputdir, multiplier=1.,
    #                    transverse=False, filetag_suffix="_"+fieldname)

    #keyname = "%s_physical" % fieldname
    #filename = datapath_db.fetch(keyname, pick='0')
    #pc.make_cube_movie(filename, "Temperature (mK)", pc.cube_frame_dir,
    #                    sigmarange=3., outputdir=outputdir, multiplier=1000.,
    #                    transverse=False, filetag_suffix="_"+fieldname,
    #                    physical=True)

def plot_gbt_diff_tests(outputdir="/cita/d/www/home/eswitzer/movies/",
                        transverse=True):
    tcv_15root = "/mnt/raid-project/gmrt/tcv/"
    tcv_15root += "modetest/73_ABCD_all_15_modes_real_maponly/"
    tcv_15map = tcv_15root + "sec_A_15hr_41-90_cleaned_clean_map_I_with_B.npy"
    tcv_15noise = tcv_15root + "sec_A_15hr_41-90_cleaned_noise_inv_I_with_B.npy"
    #ers_15root = "/mnt/raid-project/gmrt/eswitzer/GBT/"
    #ers_15root += "cleaned_maps/freq_slices_refactor_tests_15modes/"
    #ers_15map = ers_15root + "sec_A_15hr_41-90_cleaned_clean_map_I_with_B.npy"
    #ers_15noise = ers_15root + "sec_A_15hr_41-90_cleaned_noise_inv_I_with_B.npy"
    ers_15root = "./data_test/"
    ers_15map = ers_15root + "sec_A_cleaned_clean_map_I_with_B_15modes.npy"
    ers_15noise = ers_15root + "sec_A_cleaned_noise_inv_I_with_B_15modes.npy"

    pc.plot_difference(tcv_15map, ers_15map, "Temperature (mK)", sigmarange=6.,
                    fractional=False, diff_filename="./map_difference.npy",
                    outputdir=outputdir, transverse=False)

    pc.plot_difference(tcv_15noise, ers_15noise, "log inv. covariance", sigmarange=-1.,
                    multiplier=1., logscale=True, fractional=True,
                    diff_filename="./noise_inv_fractional_difference.npy",
                    outputdir=outputdir, transverse=False)

    if transverse:
        pc.plot_difference(tcv_15map, ers_15map, "Temperature (mK)", sigmarange=6.,
                        fractional=False, diff_filename="./map_difference.npy",
                        outputdir=outputdir, transverse=True)

        pc.plot_difference(tcv_15noise, ers_15noise, "log inv. covariance", sigmarange=-1.,
                        multiplier=1., logscale=True, fractional=True,
                        diff_filename="./noise_inv_fractional_difference.npy",
                        outputdir=outputdir, transverse=True)


# GBT_15hr_map_combined_cleaned_nomeanconv_0mode_map
def plot_gbt_comb_modeset(map_key, outputdir="/cita/d/www/home/eswitzer/movies/",
                          convolve=False, divider_token=";"):
    datapath_db = data_paths.DataPath()
    input_fdb = datapath_db.fetch(map_key, intend_read=True, silent=True)

    map_treatments = []
    for filekey in input_fdb[0]:
        keyparse = filekey.split(divider_token)
        if len(keyparse) == 2:
            treatment = keyparse[-1]
            if "modes" in treatment:
                map_treatments.append(treatment)

    map_treatments = unique_list(map_treatments)

    for treatment in map_treatments:
        source_key = "db:%s:map;%s" % (map_key, treatment)
        print source_key, treatment
        pc.make_cube_movie(source_key, "Temperature (mK)", pc.cube_frame_dir,
                        sigmarange=3., outputdir=outputdir, multiplier=1000.,
                        transverse=False, convolve=convolve)

        source_key = "db:%s:product;%s" % (map_key, treatment)
        pc.make_cube_movie(source_key, "Cleaned map times weights", pc.cube_frame_dir,
                        sigmarange=-1, outputdir=outputdir, multiplier=1000.,
                        transverse=False, convolve=convolve)

        source_key = "db:%s:weight;%s" % (map_key, treatment)
        pc.make_cube_movie(source_key, "inverse variance weight", pc.cube_frame_dir,
                        sigmarange=-1, outputdir=outputdir, multiplier=1.,
                        transverse=False)


def plot_sim_scheme(outputdir="/cita/d/www/home/eswitzer/movies/"):
    sim1 = "sim_streaming1.npy"
    sim2 = "sim_streaming2.npy"
    pc.plot_difference(sim1, sim2, "Temperature (mK)", sigmarange=6.,
                    fractional=False, diff_filename="./sim_difference.npy",
                    outputdir=outputdir, transverse=False)


def plot_manual(fieldname, outputdir="/cita/d/www/home/eswitzer/movies/"):
    datapath_db = data_paths.DataPath()
    file2 = './physical_cube.npy'

    keyname = "simideal_%s_physical" % fieldname
    filename = datapath_db.fetch(keyname, pick='1')
    pc.make_cube_movie(filename, "Temperature (mK)", pc.cube_frame_dir,
                        sigmarange=3., outputdir=outputdir, multiplier=1000.,
                        transverse=False, filetag_suffix="_"+fieldname,
                        physical=True)
    pc.make_cube_movie(file2, "Temperature (mK)", pc.cube_frame_dir,
                        sigmarange=3., outputdir=outputdir, multiplier=1000.,
                        transverse=False, filetag_suffix="_"+fieldname,
                        physical=True)


def plot_individual(filename, outputdir="/cita/d/www/home/eswitzer/movies/"):
    pc.make_cube_movie(filename, "Temperature (mK)", pc.cube_frame_dir,
                        sigmarange=3., outputdir=outputdir, multiplier=1000.,
                        transverse=False, filetag_suffix="_simloss",
                        physical=True)


def plot_window(outputdir="/cita/d/www/home/eswitzer/movies/"):
    file1 = './observed_window.npy'
    file2 = './physical_window.npy'
    fieldname = '15hr'

    pc.make_cube_movie(file1, "Window", pc.cube_frame_dir,
                        sigmarange=-1, outputdir=outputdir, multiplier=1.,
                        transverse=False, filetag_suffix="_"+fieldname,
                        physical=True)
    pc.make_cube_movie(file2, "Window", pc.cube_frame_dir,
                        sigmarange=-1, outputdir=outputdir, multiplier=1.,
                        transverse=False, filetag_suffix="_"+fieldname,
                        physical=True)


def plot_wigglez(fieldname, outputdir="/cita/d/www/home/eswitzer/movies/",
                 complete=False):
    datapath_db = data_paths.DataPath()
    if complete:
        ctag = "complete_"
    else:
        ctag = ""

    db_key = "WiggleZ_%s_%sdelta_binned_data" % (fieldname, ctag)
    filename = datapath_db.fetch(db_key)
    pc.make_cube_movie(filename, "counts", pc.cube_frame_dir,
                        sigmarange=-1, outputdir=outputdir, multiplier=1.,
                        transverse=False, filetag_suffix="_"+fieldname)

    db_key = "WiggleZ_%s_%sbinned_data" % (fieldname, ctag)
    filename = datapath_db.fetch(db_key)
    pc.make_cube_movie(filename, "counts", pc.cube_frame_dir,
                        sigmarange=-1, outputdir=outputdir, multiplier=1.,
                        transverse=False, filetag_suffix="_"+fieldname)

    db_key = "WiggleZ_%s_%sselection" % (fieldname, ctag)
    filename = datapath_db.fetch(db_key)
    pc.make_cube_movie(filename, "selection", pc.cube_frame_dir,
                        sigmarange=-1, outputdir=outputdir, multiplier=1.,
                        transverse=False, filetag_suffix="_"+fieldname)

    db_key = "WiggleZ_%s_%sseparable_selection" % (fieldname, ctag)
    filename = datapath_db.fetch(db_key)
    pc.make_cube_movie(filename, "selection", pc.cube_frame_dir,
                        sigmarange=-1, outputdir=outputdir, multiplier=1.,
                        transverse=False, filetag_suffix="_"+fieldname)

    db_key = "WiggleZ_%s_%smontecarlo" % (fieldname, ctag)
    filename = datapath_db.fetch(db_key)
    pc.make_cube_movie(filename, "selection", pc.cube_frame_dir,
                        sigmarange=-1, outputdir=outputdir, multiplier=1.,
                        transverse=False, filetag_suffix="_"+fieldname)


def plot_mode_amplitudes(mapkey, outputdir="/cita/d/www/home/eswitzer/movies/"):
    datapath_db = data_paths.DataPath()
    map_cases = datapath_db.fileset_cases(mapkey, "pair;product;treatment")
    for pair in map_cases['pair']:
        source_key = "db:%s:%s;modes;100modes" % (mapkey, pair)
        pc.make_cube_movie(source_key, "Mode amp", pc.cube_frame_dir,
                           sigmarange=-1, outputdir=outputdir, multiplier=1.,
                           transverse=False, convolve=False, logscale=True)


if __name__ == "__main__":

    #plot_mode_amplitudes('GBT_15hr_map_fdgcal_cleaned')
    #plot_mode_amplitudes('GBT_15hr_map_fdgcal_cleaned_noconv')
    plot_gbt_mapset()
    #plot_cleaned_gbt15hr_mapset()
    #plot_gbt_simset('sim_15hr')
    #plot_gbt_simset('simvel_15hr')
    #plot_gbt_simset('simideal_15hr')
    #plot_gbt_simset('22hr')
    #plot_gbt_newmapset()

    #plot_gbt_comb_modeset('GBT_1hr_map_oldcalkiyo_cleaned_combined', convolve=True,
    #                      outputdir="/mnt/raid-project/gmrt/eswitzer/GBT/movies/")

    #plot_gbt_comb_modeset('GBT_15hr_optimalmap_cleaned_combined', convolve=True,
    #                      outputdir="/mnt/raid-project/gmrt/eswitzer/GBT/movies/")

    #plot_gbt_comb_modeset('GBT_15hr_optimalmap_selfcal_762_cleaned_combined', convolve=True)
    #plot_gbt_comb_modeset('GBT_15hr_optimalmap_selfcal_762_cleaned_combined', convolve=False)
    #plot_gbt_comb_modeset('GBT_15hr_map_mapcal2_cleaned_combined', convolve=True)
    #plot_gbt_comb_modeset('GBT_15hr_map_mapcal2_cleaned_combined', convolve=False)
    #plot_gbt_comb_modeset('GBT_15hr_map_mapcal_cleaned_combined', convolve=True)
    #plot_gbt_comb_modeset('GBT_15hr_map_mapcal_cleaned_combined', convolve=False)

    #plot_gbt_comb_modeset('GBT_15hr_map_fluxcal_cleaned_noconv_combined', convolve=False)
    #plot_gbt_comb_modeset('GBT_15hr_map_fluxcal_cleaned_noconv_combined', convolve=True)
    #plot_gbt_comb_modeset('GBT_15hr_map_fluxcal_cleaned_combined', convolve=False)
    #plot_gbt_comb_modeset('GBT_15hr_map_fluxcal_cleaned_sims_noconv_combined', convolve=False)
    #plot_gbt_comb_modeset('GBT_15hr_map_fluxcal_cleaned_sims_noconv_combined', convolve=True)
    #plot_gbt_comb_modeset('GBT_15hr_map_fluxcal_cleaned_sims_combined', convolve=False)

    #plot_gbt_comb_modeset('GBT_15hr_optimalmap_mapv2fdgcal_cleaned_noconv_combined', convolve=False)
    #plot_gbt_comb_modeset('GBT_15hr_optimalmap_mapv2fdgcal_cleaned_noconv_combined', convolve=True)
    #plot_gbt_comb_modeset('GBT_15hr_optimalmap_mapv2fdgcal_cleaned_combined', convolve=False)
    #plot_gbt_comb_modeset('GBT_15hr_optimalmap_mapv2oldcal_cleaned_noconv_combined', convolve=False)
    #plot_gbt_comb_modeset('GBT_15hr_optimalmap_mapv2oldcal_cleaned_noconv_combined', convolve=True)
    #plot_gbt_comb_modeset('GBT_15hr_optimalmap_mapv2oldcal_cleaned_combined', convolve=False)

    #plot_gbt_comb_modeset('GBT_15hr_map_fdgcal_cleaned_noconv_combined', convolve=False)
    #plot_gbt_comb_modeset('GBT_15hr_map_fdgcal_cleaned_noconv_combined', convolve=True)
    #plot_gbt_comb_modeset('GBT_15hr_map_fdgcal_cleaned_combined', convolve=False)
    #plot_gbt_comb_modeset('GBT_15hr_map_fdgcal_cleaned_sims_noconv_combined', convolve=False)
    #plot_gbt_comb_modeset('GBT_15hr_map_fdgcal_cleaned_sims_noconv_combined', convolve=True)
    #plot_gbt_comb_modeset('GBT_15hr_map_fdgcal_cleaned_sims_combined', convolve=False)

    #plot_wigglez('15hr', complete=False)
    #plot_wigglez('22hr', complete=False)
    #plot_wigglez('1hr', complete=True)
    #plot_wigglez('22hr')
    #plot_wigglez('1hr')

    # random odds and ends
    #plot_individual("/mnt/raid-project/gmrt/eswitzer/GBT/cleaned_maps/15hr_sim/sec_A_cleaned_clean_map_I_with_B_0modes.npy")
    #plot_window()
    #plot_gbt_diff_tests()
    #plot_sim_scheme()
    #plot_manual('15hr')

