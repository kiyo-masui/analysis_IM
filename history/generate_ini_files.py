r"""
This writes out ini files for the map cleaning in an automated way
"""
from kiyopy import parse_ini
from utils import file_tools
# TODO what to do with this one?
# GBT_15hr_map_fdgcal_cleanedplussim_cleaned.ini

def write_map_cleanerini(mapname, cutlist, nfreq, factorizable=True, meansub=True,
                         regenerate=False, convolve=False,
                         subtract_inputmap_from_sim = True,
                         subtract_sim_from_inputmap = False,
                         sim_multiplier = 1., username="Eric",
                         modes = range(0, 105, 5), simfile=None, prefix="fs_",
                         inidir="./input/ers/map_cleaning_autogen/"):
    file_tools.mkparents(inidir)
    params = {}

    simtag = ""
    if simfile:
        simtag = "_plussim"
        if subtract_inputmap_from_sim:
            simtag = "_plussim_minusmap"
        if subtract_sim_from_inputmap:
            simtag = "_plussim_minussim"

    alt = ""
    if sim_multiplier != 1.:
        multstring = "%5.3g" % sim_multiplier
        alt = "_simx" + multstring.replace(".", "p").strip()

    key = '%s_cleaned%s%s' % (mapname, simtag, alt)
    params["output_root"] = '%s_path_%s' % (key, username)

    params["SVD_root"] = None
    params["modes"] = modes
    params["map1"] = mapname
    params["map2"] = mapname
    params["noise_inv1"] = mapname
    params["noise_inv2"] = mapname
    params["no_weights"] = False
    params["sub_weighted_mean"] = meansub
    params["factorizable_noise"] = factorizable
    params["convolve"] = convolve
    params["regenerate_noise_inv"] = regenerate
    params["freq_list"] = tuple([ind for ind in range(nfreq) \
                                 if ind not in cutlist])

    if simfile:
        params["simfile"] = simfile
        params["sim_multiplier"] = sim_multiplier
        params["subtract_inputmap_from_sim"] = subtract_inputmap_from_sim
        params["subtract_sim_from_inputmap"] = subtract_sim_from_inputmap

    filename = "%s/%s.ini" % (inidir, params["output_root"])
    print filename
    parse_ini.write_params(params, filename, prefix=prefix)

def write_map_cleanerini_old(mapname, cutlist, nfreq, factorizable=True, meansub=True,
                         regenerate=False, noconv=False,
                         subtract_inputmap_from_sim = True,
                         subtract_sim_from_inputmap = False,
                         modes = range(0, 105, 5), simfile=None, prefix="fs_",
                         inidir="./input/ers/map_cleaning_autogen_old/"):
    file_tools.mkparents(inidir)

    params = {}
    params["SVD_root"] = None
    params["modes"] = modes
    params["map1"] = mapname
    params["map2"] = mapname
    params["noise_inv1"] = mapname
    params["noise_inv2"] = mapname
    params["no_weights"] = False
    params["sub_weighted_mean"] = meansub
    params["factorizable_noise"] = factorizable
    params["regenerate_noise_inv"] = regenerate
    params["freq_list"] = tuple([ind for ind in range(nfreq) \
                                 if ind not in cutlist])

    tag = "_sims" if simfile else ""
    if simfile:
        params["simfile"] = simfile
        params["sim_multiplier"] = 1.
        params["subtract_inputmap_from_sim"] = subtract_inputmap_from_sim
        params["subtract_sim_from_inputmap"] = subtract_sim_from_inputmap

    params["convolve"] = False
    # TODO: move this to direct path rather than db
    params["output_root"] = "%s_cleaned%s_noconv_path_Eric" % (mapname, tag)
    filename = "%s/%s_cleaned%s_noconv.ini" % (inidir, mapname, tag)
    parse_ini.write_params(params, filename, prefix=prefix)

    params["convolve"] = True
    params["output_root"] = "%s_cleaned%s_path_Eric" % (mapname, tag)
    filename = "%s/%s_cleaned%s.ini" % (inidir, mapname, tag)
    parse_ini.write_params(params, filename, prefix=prefix)


def gen_old_inis():
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
               "GBT_15hr_map_oldcal_plussim_10pct",
               "GBT_1hr_map_fluxpolcal",
               "GBT_22hr_map_fluxpolcal"]

    for mapname in maplist:
        write_map_cleanerini_old(mapname, cutlist, 256, simfile=None)
        write_map_cleanerini_old(mapname, cutlist, 256, simfile=simfile)

    #simfile = '/mnt/raid-project/gmrt/eswitzer/GBT/simulations/15hr_optimalmap_str/sim_beam_000.npy'
    simfile = '/mnt/raid-project/gmrt/eswitzer/GBT/simulations/15hr_optimalmap762_str/sim_beam_000.npy'
    maplist = ["GBT_15hr_optimalmap_fluxpolcal",
               "GBT_15hr_optimalmap_mapv2fdgcal",
               "GBT_15hr_optimalmap_mapv2oldcal"]
    cutlist = []

    for mapname in maplist:
        write_map_cleanerini_old(mapname, cutlist, 40, simfile=None)
        write_map_cleanerini_old(mapname, cutlist, 40, simfile=simfile)


def gen_new_inis(sim_multiplier=1.):
    cutlist = [6, 7, 8, 15, 16, 18, 19, 20, 21, 22, 37, 80, 103, 104, 105, 106, \
               107, 108, 130, 131, 132, 133, 134, 171, 175, 177, 179, 182, 183, \
               187, 189, 192, 193, 194, 195, 196, 197, 198, 201, 204, 208, 209, \
               212, 213, 218, 219, 229, 233, 237, 244, 254, 255]

    simfile = '/mnt/raid-project/gmrt/eswitzer/GBT/simulations/15hr_oldmap_str/sim_beam_000.npy'
    maplist = ["GBT_15hr_map_oldcal", "GBT_15hr_map_mapcal",
               "GBT_15hr_map_mapcal2"]

    for mapname in maplist:
        write_map_cleanerini(mapname, cutlist, 256, simfile=None)

        write_map_cleanerini(mapname, cutlist, 256, simfile=simfile,
                             subtract_inputmap_from_sim=False,
                             subtract_sim_from_inputmap=False,
                             sim_multiplier=sim_multiplier)

        write_map_cleanerini(mapname, cutlist, 256, simfile=simfile,
                             subtract_inputmap_from_sim=True,
                             subtract_sim_from_inputmap=False,
                             sim_multiplier=sim_multiplier)

        write_map_cleanerini(mapname, cutlist, 256, simfile=simfile,
                             subtract_inputmap_from_sim=False,
                             subtract_sim_from_inputmap=True,
                             sim_multiplier=sim_multiplier)

    simfile = '/mnt/raid-project/gmrt/eswitzer/GBT/simulations/15hr_optimalmap762_str/sim_beam_000.npy'
    maplist = ["GBT_15hr_optimalmap_selfcal_762"]
    cutlist = []

    for mapname in maplist:
        write_map_cleanerini(mapname, cutlist, 40, simfile=None)

        write_map_cleanerini(mapname, cutlist, 40, simfile=simfile,
                             subtract_inputmap_from_sim=False,
                             subtract_sim_from_inputmap=False,
                             sim_multiplier=sim_multiplier)

        write_map_cleanerini(mapname, cutlist, 40, simfile=simfile,
                             subtract_inputmap_from_sim=True,
                             subtract_sim_from_inputmap=False,
                             sim_multiplier=sim_multiplier)

        write_map_cleanerini(mapname, cutlist, 40, simfile=simfile,
                             subtract_inputmap_from_sim=False,
                             subtract_sim_from_inputmap=True,
                             sim_multiplier=sim_multiplier)


if __name__ == '__main__':
    #gen_old_inis()
    gen_new_inis()
    gen_new_inis(sim_multiplier=1.1)
