from plotting import plot_cube as pc

# TODO: delete tag from invocation
# TODO: cube root invocation
# TODO: fix error bars in selection function
# TODO: define output directory for movies
cube_frame_dir = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/cube_frames/"
output_dir = ""

def plot_15hr_selection():
    """make plots of the 15hr selection
    """
    root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned_wiggleZ/"
    pc.make_cube_movie(root_directory + "reg15montecarlo.npy",
                    "reg15montecarlo",
                    "MC Selection", cube_frame_dir + "reg15montecarlo",
                    sigmarange=None)
    root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned_wiggleZ/"
    pc.make_cube_movie(root_directory + "reg15separable.npy",
                    "15hr_Wigglez_separable_selection",
                    "# galaxies/pixel", cube_frame_dir + "15hr_Wigglez_separable_selection")
    root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned_wiggleZ/"
    pc.make_cube_movie(root_directory + "reg15selection.npy",
                    "15hr_Wigglez_selection",
                    "# galaxies/pixel", cube_frame_dir + "15hr_Wigglez_selection")

    pc.make_cube_movie("delta_selection.npy",
                    "delta_selection",
                    "Difference in selection", cube_frame_dir + "delta_selection")
    pc.make_cube_movie("reg15selection_est.npy",
                    "reg15selection_est",
                    "Full selection", cube_frame_dir + "reg15selection_est",
                    sigmarange=None)

def plot_15hr_selection():
    """make plots of the 22hr selection
    """
    root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned_wiggleZ/"

    pc.make_cube_movie(root_directory + "reg22montecarlo.npy",
                    "reg22montecarlo",
                    "MC Selection", cube_frame_dir + "reg22montecarlo",
                    sigmarange=None)

    pc.make_cube_movie(root_directory + "reg22selection.npy",
                    "reg22selection",
                    "Selection", cube_frame_dir + "reg22selection",
                    sigmarange=None)

    pc.make_cube_movie(root_directory + "reg22separable.npy",
                    "reg22separable",
                    "Separable Selection", cube_frame_dir + "reg22separable",
                    sigmarange=None)

# make plots of the 15hr field
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/combined_maps/"
pc.make_cube_movie(root_directory + "combined_41-73_cleaned_clean_test.npy",
                "combined_41-73_cleaned_clean_test",
                "Temperature (mK)", cube_frame_dir + "combined_41-73_cleaned_clean_test",
                sigmarange=6., multiplier = 1000.)
pc.make_cube_movie(root_directory + "combined_41-73_cleaned_noise_inv_test.npy",
                "combined_41-73_cleaned_noise_inv_test",
                "Covariance", cube_frame_dir + "combined_41-73_cleaned_noise_inv_test",
                sigmarange=-1)

# make plots of the 22hr field
#root_directory = "/mnt/raid-project/gmrt/calinliv/wiggleZ/corr/84_ABCD_all_15_modes/"
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/combined_maps/"
pc.make_cube_movie(root_directory + "combined_22hr_41-84_cleaned_clean.npy",
                "combined_22hr_41-84_cleaned_clean.npy",
                "Temperature (mK)", cube_frame_dir + "combined_22hr_41-84_cleaned_clean.npy",
                sigmarange=6, multiplier = 1000.)

pc.make_cube_movie(root_directory + "combined_22hr_41-84_cleaned_noise_inv.npy",
                "combined_22hr_41-84_cleaned_noise_inv.npy",
                "Covariance", cube_frame_dir + "combined_22hr_41-84_cleaned_noise_inv.npy",
                sigmarange=-1)

sys.exit()
root_directory = "/mnt/raid-project/gmrt/calinliv/wiggleZ/simulations/test100/"
pc.make_cube_movie(root_directory + "simulated_signal_map_1.npy",
                "simulated_signal_map_1",
                "Temperature", cube_frame_dir + "simulated_signal_map_1",
                sigmarange=[-1.,1.])
root_directory = "/mnt/raid-project/gmrt/calinliv/wiggleZ/simulations/test100/"
pc.make_cube_movie(root_directory + "simulated_signal_map_1_with_beam.npy",
                "simulated_signal_map_1_with_beam",
                "Temperature", cube_frame_dir + "simulated_signal_map_1_with_beam",
                sigmarange=[-1.,1.])
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/modetest/73_ABCD_all_0_modes_sim_maponly_NOCONV/"
pc.make_cube_movie(root_directory + "sec_A_15hr_41-73_cleaned_clean_map_I_with_B.npy",
                "sec_A_15hr_41-73_cleaned_clean_map_I_with_B_noconv",
                "Temperature", cube_frame_dir + "sec_A_15hr_41-73_cleaned_clean_map_I_with_B_noconv",
                sigmarange=[-1.,1.])
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/modetest/73_ABCD_all_0_modes_sim_maponly/"
pc.make_cube_movie(root_directory + "sec_A_15hr_41-73_cleaned_clean_map_I_with_B.npy",
                "sec_A_15hr_41-73_cleaned_clean_map_I_with_B",
                "Temperature", cube_frame_dir + "sec_A_15hr_41-73_cleaned_clean_map_I_with_B",
                sigmarange=[-1.,1.])
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/modetest_combined_maps_0_50/"
pc.make_cube_movie(root_directory + "combined_sim_41-73_cleaned_clean_0.npy",
                "combined_sim_41-73_cleaned_clean_0",
                "Temperature", cube_frame_dir + "combined_sim_41-73_cleaned_clean_0",
                sigmarange=[-1.,1.])
sys.exit()

root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/noise_model/"
pc.make_cube_movie(root_directory + "completeness_model_41-73_15modes.npy",
                "completeness_model_41-73_15modes",
                "Completeness", cube_frame_dir + "completeness_model_41-73_15modes",
                sigmarange=5., ignore=-1.)
sys.exit()
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/noise_model/"
pc.make_cube_movie(root_directory + "noise_model_41-73_15modes.npy",
                "noise_model_41-73_15modes",
                "Temperature", cube_frame_dir + "noise_model_41-73_15modes",
                sigmarange=5., ignore=-1.)
sys.exit()

root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/combined_maps/"
pc.make_cube_movie(root_directory + "combined_41-73_clean_test.npy",
                "combined_41-73_clean_test",
                "Temperature", cube_frame_dir + "combined_41-73_clean_test",
                sigmarange=8.)
pc.make_cube_movie(root_directory + "combined_41-73_noise_inv_test.npy",
                "combined_41-73_noise_inv_test",
                "Covariance", cube_frame_dir + "combined_41-73_noise_inv_test",
                sigmarange=-1)
sys.exit()

root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/combined_maps/"
pc.make_cube_movie(root_directory + "combined_41-73_cleaned_clean_test.npy",
                "combined_41-73_cleaned_clean_test",
                "Temperature", cube_frame_dir + "combined_41-73_cleaned_clean_test",
                sigmarange=4.)
pc.make_cube_movie(root_directory + "combined_41-73_cleaned_noise_inv_test.npy",
                "combined_41-73_cleaned_noise_inv_test",
                "Covariance", cube_frame_dir + "combined_41-73_cleaned_noise_inv_test",
                sigmarange=-1)
sys.exit()

root_directory = "/mnt/raid-project/gmrt/kiyo/wiggleZ/corr/"
pc.make_cube_movie(root_directory + \
                "sec_B_15hr_41-69_cleaned_clean_map_I_with_A.npy",
                "sec_B_15hr_41-69_cleaned_clean_map_I_with_A",
                "Temperature",
                cube_frame_dir + "sec_B_15hr_41-69_cleaned_clean_map_I_with_A",
                sigmarange=[-0.001, 0.001])
pc.make_cube_movie(root_directory + \
                "sec_A_15hr_41-69_cleaned_clean_map_I_with_B.npy",
                "sec_A_15hr_41-69_cleaned_clean_map_I_with_B",
                "Temperature",
                cube_frame_dir + "sec_A_15hr_41-69_cleaned_clean_map_I_with_B",
                sigmarange=[-0.001, 0.001])
sys.exit()

root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/combined_maps/"
pc.make_cube_movie(root_directory + "combined_41-69_cleaned_product.npy",
                "combined_41-69_cleaned_product",
                "Temperature", cube_frame_dir + "combined_41-69_cleaned_product",
                sigmarange=20.)
pc.make_cube_movie(root_directory + "combined_41-69_cleaned_clean.npy",
                "combined_41-69_cleaned_clean",
                "Temperature", cube_frame_dir + "combined_41-69_cleaned_clean",
                sigmarange=[-0.01, 0.01])
pc.make_cube_movie(root_directory + "combined_41-69_cleaned_noise_inv.npy",
                "combined_41-69_cleaned_noise_inv",
                "Covariance", cube_frame_dir + "combined_41-69_cleaned_noise_inv",
                sigmarange=-1)
sys.exit()

root_directory = "/mnt/raid-project/gmrt/kiyo/wiggleZ/corr/"
pc.make_cube_movie(root_directory + "sec_A_15hr_41-69_cleaned_clean_map_I.npy",
                "sec_A_15hr_41-69_cleaned_clean_map_I",
                "Temperature", cube_frame_dir + "sec_A_15hr_41-69_cleaned_clean_map_I",
                sigmarange=[-0.01, 0.01])
pc.make_cube_movie(root_directory + "sec_A_15hr_41-69_cleaned_noise_inv_I.npy",
                "sec_A_15hr_41-69_cleaned_noise_inv_I",
                "Covariance", cube_frame_dir + "sec_A_15hr_41-69_cleaned_noise_inv_I",
                sigmarange=-1)
pc.make_cube_movie(root_directory + "sec_B_15hr_41-69_cleaned_clean_map_I.npy",
                "sec_B_15hr_41-69_cleaned_clean_map_I",
                "Temperature", cube_frame_dir + "sec_B_15hr_41-69_cleaned_clean_map_I",
                sigmarange=[-0.01, 0.01])
pc.make_cube_movie(root_directory + "sec_B_15hr_41-69_cleaned_noise_inv_I.npy",
                "sec_B_15hr_41-69_cleaned_noise_inv_I",
                "Covariance", cube_frame_dir + "sec_B_15hr_41-69_cleaned_noise_inv_I",
                sigmarange=-1)
sys.exit()

root_directory = "/mnt/raid-project/gmrt/calinliv/wiggleZ/corr/test/"
maplist = ["sec_A_15hr_41-73_cleaned_clean_map_I_with_B",
        "sec_A_15hr_41-73_cleaned_clean_map_I_with_C",
        "sec_A_15hr_41-73_cleaned_clean_map_I_with_D",
        "sec_B_15hr_41-73_cleaned_clean_map_I_with_A",
        "sec_B_15hr_41-73_cleaned_clean_map_I_with_C",
        "sec_B_15hr_41-73_cleaned_clean_map_I_with_D",
        "sec_C_15hr_41-73_cleaned_clean_map_I_with_A",
        "sec_C_15hr_41-73_cleaned_clean_map_I_with_B",
        "sec_C_15hr_41-73_cleaned_clean_map_I_with_D",
        "sec_D_15hr_41-73_cleaned_clean_map_I_with_A",
        "sec_D_15hr_41-73_cleaned_clean_map_I_with_B",
        "sec_D_15hr_41-73_cleaned_clean_map_I_with_C"]
covlist = ["sec_A_15hr_41-73_cleaned_noise_inv_I_with_B",
        "sec_A_15hr_41-73_cleaned_noise_inv_I_with_C",
        "sec_A_15hr_41-73_cleaned_noise_inv_I_with_D",
        "sec_B_15hr_41-73_cleaned_noise_inv_I_with_A",
        "sec_B_15hr_41-73_cleaned_noise_inv_I_with_C",
        "sec_B_15hr_41-73_cleaned_noise_inv_I_with_D",
        "sec_C_15hr_41-73_cleaned_noise_inv_I_with_A",
        "sec_C_15hr_41-73_cleaned_noise_inv_I_with_B",
        "sec_C_15hr_41-73_cleaned_noise_inv_I_with_D",
        "sec_D_15hr_41-73_cleaned_noise_inv_I_with_A",
        "sec_D_15hr_41-73_cleaned_noise_inv_I_with_B",
        "sec_D_15hr_41-73_cleaned_noise_inv_I_with_C"]

for tagname in maplist:
    pc.make_cube_movie(root_directory + tagname + ".npy",
                tagname, "Temperature", cube_frame_dir + tagname)
for tagname in covlist:
    pc.make_cube_movie(root_directory + tagname + ".npy",
                tagname, "Covariance", cube_frame_dir + tagname)

