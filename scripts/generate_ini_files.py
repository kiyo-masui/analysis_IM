r"""
This writes out ini files for the map cleaning in an automated way
"""
from kiyopy import parse_ini
from utils import file_tools
# TODO what to do with this one?
# GBT_15hr_map_fdgcal_cleanedplussim_cleaned.ini

def write_map_cleanerini(mapname, cutlist, nfreq, factorizable=True, meansub=True,
                         modes = range(0, 105, 5), simfile=None, prefix="fs_",
                         inidir="./input/ers/map_cleaning_autogen/"):
    file_tools.mkparents(inidir)

    params = {}
    params["modes"] = modes
    params["map1"] = mapname
    params["map2"] = mapname
    params["noise_inv1"] = mapname
    params["noise_inv2"] = mapname
    params["no_weights"] = False
    params["sub_weighted_mean"] = meansub
    params["factorizable_noise"] = factorizable
    params["freq_list"] = tuple([ind for ind in range(nfreq) \
                                 if ind not in cutlist])

    tag = "_sims" if simfile else ""
    if simfile:
        params["simfile"] = simfile

    params["convolve"] = False
    params["output_root"] = "%s_cleaned%s_noconv_path_Eric" % (mapname, tag)
    filename = "%s/%s_cleaned%s_noconv.ini" % (inidir, mapname, tag)
    parse_ini.write_params(params, filename, prefix=prefix)

    params["convolve"] = True
    params["output_root"] = "%s_cleaned%s_path_Eric" % (mapname, tag)
    filename = "%s/%s_cleaned%s.ini" % (inidir, mapname, tag)
    parse_ini.write_params(params, filename, prefix=prefix)

if __name__ == '__main__':
    cutlist = [6, 7, 8, 15, 16, 18, 19, 20, 21, 22, 37, 80, 103, 104, 105, 106, \
               107, 108, 130, 131, 132, 133, 134, 171, 175, 177, 179, 182, 183, \
               187, 189, 192, 193, 194, 195, 196, 197, 198, 201, 204, 208, 209, \
               212, 213, 218, 219, 229, 233, 237, 244, 254, 255]

    simfile = '/mnt/raid-project/gmrt/eswitzer/GBT/simulations/15hr_oldmap_str/sim_beam_000.npy'
    maplist = ["GBT_15hr_map_fdgcal",
               "GBT_15hr_map_fdgcal_plussim",
               "GBT_15hr_map_fluxcal",
               "GBT_15hr_map_fluxpolcal",
               "GBT_15hr_map_oldcal",
               "GBT_15hr_map_oldcal_plussim",
               "GBT_1hr_map_fluxpolcal",
               "GBT_22hr_map_fluxpolcal"]

    for mapname in maplist:
        write_map_cleanerini(mapname, cutlist, 256, simfile=None)
        write_map_cleanerini(mapname, cutlist, 256, simfile=simfile)

    simfile = '/mnt/raid-project/gmrt/eswitzer/GBT/simulations/15hr_optimalmap_str/sim_beam_000.npy'
    maplist = ["GBT_15hr_optimalmap_fluxpolcal",
               "GBT_15hr_optimalmap_mapv2fdgcal",
               "GBT_15hr_optimalmap_mapv2oldcal"]

    for mapname in maplist:
        write_map_cleanerini(mapname, cutlist, 120, simfile=None)
        write_map_cleanerini(mapname, cutlist, 120, simfile=simfile)


