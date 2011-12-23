r"""Several scripts to call plot_cube_movie and make movies of GBT data,
simulations, etc.
"""
from plotting import plot_cube as pc
from utils import data_paths

def plot_gbt_mapset(outputdir="/cita/d/www/home/eswitzer/movies/"):
    #pc.plot_gbt_maps('GBT_15hr_map_cleaned_20mode', outputdir=outputdir, transverse=False)
    #pc.plot_gbt_maps('GBT_15hr_map_cleaned_0mode', outputdir=outputdir, transverse=False)
    #pc.plot_gbt_maps('GBT_22hr_map_cleaned_20mode', outputdir=outputdir, transverse=False)
    #pc.plot_gbt_maps('GBT_22hr_map_cleaned_0mode', outputdir=outputdir, transverse=False)

    #pc.plot_gbt_maps('GBT_15hr_map_proposal', transverse=True,
    #                 outputdir=outputdir, skip_noise=True)
    #pc.plot_gbt_maps('GBT_15hr_map_proposal', transverse=False,
    #                 outputdir=outputdir, skip_noise=True)
    ##pc.plot_gbt_maps('GBT_15hr_map', outputdir=outputdir, transverse=True)
    #pc.plot_gbt_maps('GBT_15hr_map', outputdir=outputdir, transverse=False)
    ##pc.plot_gbt_maps('GBT_22hr_map', outputdir=outputdir, transverse=True)
    #pc.plot_gbt_maps('GBT_22hr_map', outputdir=outputdir, transverse=False)
    ##pc.plot_gbt_maps('GBT_1hr_map', outputdir=outputdir, transverse=True)
    pc.plot_gbt_maps('GBT_1hr_map', outputdir=outputdir, transverse=False)


def plot_gbt_simset(fieldname, outputdir="/cita/d/www/home/eswitzer/movies/"):
    datapath_db = data_paths.DataPath()

    #keyname = "sim_%s" % fieldname
    keyname = "simideal_%s" % fieldname
    filename = datapath_db.fetch(keyname, pick='1')
    pc.make_cube_movie(filename, "Temperature (mK)", pc.cube_frame_dir,
                        sigmarange=3., outputdir=outputdir, multiplier=1000.,
                        transverse=False, filetag_suffix="_"+fieldname)

    #keyname = "sim_%s_beam" % fieldname
    keyname = "simideal_%s_beam" % fieldname
    filename = datapath_db.fetch(keyname, pick='1')
    pc.make_cube_movie(filename, "Temperature (mK)", pc.cube_frame_dir,
                        sigmarange=3., outputdir=outputdir, multiplier=1000.,
                        transverse=False, filetag_suffix="_"+fieldname)


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


def plot_gbt_comb_modeset(fieldname, outputdir="/cita/d/www/home/eswitzer/movies/"):
    datapath_db = data_paths.DataPath()

    for modenum in range(0, 55, 5):
        #keyname = "GBT_%s_combined_cleaned_%dmode_map" % (fieldname, modenum)
        #filename = datapath_db.fetch(keyname)
        #pc.make_cube_movie(filename, "Temperature (mK)", pc.cube_frame_dir,
        #                sigmarange=2.5, outputdir=outputdir, multiplier=1000.,
        #                transverse=False)

        keyname = "GBT_%s_combined_cleaned_%dmode_product" % \
                  (fieldname, modenum)
        filename = datapath_db.fetch(keyname)
        pc.make_cube_movie(filename, "Cleaned map times weights", pc.cube_frame_dir,
                        sigmarange=-1, outputdir=outputdir, multiplier=1000.,
                        transverse=False)

        #keyname = "GBT_%s_combined_cleaned_%dmode_weight" % \
        #          (fieldname, modenum)
        #filename = datapath_db.fetch(keyname)
        #pc.make_cube_movie(filename, "inverse variance weight", pc.cube_frame_dir,
        #                sigmarange=2.5, outputdir=outputdir, multiplier=1.,
        #                transverse=False)


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

    #db_key = "WiggleZ_%s_%sbinned_data" % (fieldname, ctag)
    #filename = datapath_db.fetch(db_key)
    #pc.make_cube_movie(filename, "counts", pc.cube_frame_dir,
    #                    sigmarange=-1, outputdir=outputdir, multiplier=1.,
    #                    transverse=False, filetag_suffix="_"+fieldname)

    #db_key = "WiggleZ_%s_%sselection" % (fieldname, ctag)
    #filename = datapath_db.fetch(db_key)
    #pc.make_cube_movie(filename, "selection", pc.cube_frame_dir,
    #                    sigmarange=-1, outputdir=outputdir, multiplier=1.,
    #                    transverse=False, filetag_suffix="_"+fieldname)

    #db_key = "WiggleZ_%s_%sseparable_selection" % (fieldname, ctag)
    #filename = datapath_db.fetch(db_key)
    #pc.make_cube_movie(filename, "selection", pc.cube_frame_dir,
    #                    sigmarange=-1, outputdir=outputdir, multiplier=1.,
    #                    transverse=False, filetag_suffix="_"+fieldname)

    db_key = "WiggleZ_%s_%smontecarlo" % (fieldname, ctag)
    filename = datapath_db.fetch(db_key)
    pc.make_cube_movie(filename, "selection", pc.cube_frame_dir,
                        sigmarange=-1, outputdir=outputdir, multiplier=1.,
                        transverse=False, filetag_suffix="_"+fieldname)


if __name__ == "__main__":
    #plot_gbt_mapset()
    #plot_gbt_comb_modeset('15hr')
    #plot_gbt_comb_modeset('22hr')
    #plot_gbt_simset('15hr')
    #plot_gbt_simset('22hr')

    #plot_gbt_diff_tests()
    #plot_sim_scheme()

    #plot_manual('15hr')
    #plot_window()
    #plot_wigglez('15hr', complete=False)
    plot_wigglez('22hr', complete=False)
    #plot_wigglez('1hr', complete=True)
    #plot_wigglez('22hr')
    #plot_wigglez('1hr')
