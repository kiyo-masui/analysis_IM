"""
    Procedure for adding new files to the path database:
    (or email Eric Switzer eswitzer@cita.utoronto.ca)

    1. use one of the macros to add a data path/file/set
    2. update the version_tag to today's date
    3. if you want: run utils/data_paths.py, go to the web directory and type
    './build' to convert the markdown to html for documentation.

    Note that it is easier to debug this file by executing it directly
    rather than running it through the database parser.

    Environment variables for prawn/chime:
    export GBTDATA_ESWITZER='/mnt/raid-project/gmrt/eswitzer/'
    export GBTDATA_ADLEWIS='/mnt/raid-project/gmrt/alewis/'
    export GBTDATA_CALINLIV='/mnt/raid-project/gmrt/calinliv/'
    export GBTDATA_KIYO='/mnt/raid-project/gmrt/kiyo/'
    export GBTDATA_NBANAVAR='/mnt/raid-project/gmrt/nbanavar/'
    export GBTDATA_PEN='/mnt/raid-project/gmrt/pen/'
    export GBTDATA_TCHANG='/mnt/raid-project/gmrt/tchang/'
    export GBTDATA_TCV='/mnt/raid-project/gmrt/tcv/'

    This is just a script.
"""
# this is executed directly and pdb is handed into the DataPath class
# for documentation on the database format, see utils/data_paths.py
import os
import sys
import path_forms
import shelve

# local path is /cita/d/www/home/eswitzer/GBT_param/path_database.py
db_url = "http://www.cita.utoronto.ca/~eswitzer/GBT_param/path_database.py"
version_tag = "Feb 15st, 2012"
pdb = {}
groups = {}
dbcl = path_forms.PathForms(verbose=False)
field_list = ['15hr', '22hr', '1hr']

dbcl.register_list_empty_groups(['GBTmaps', 'GBTcleaned',
                                 'WiggleZ', 'Simulations'])

#-----------------------------------------------------------------------------
# user directories
#-----------------------------------------------------------------------------
userdirs = [("GBTDATA_NBANAVAR", "Nidhi B."),
            ("GBTDATA_CALINLIV", "Liviu C."),
            ("GBTDATA_TCHANG", "Tzu-Ching C."),
            ("GBTDATA_ADLEWIS", "Adam L."),
            ("GBTDATA_KIYO", "Kiyo M."),
            ("GBTDATA_PEN", "Ue-Li P."),
            ("GBTDATA_ESWITZER", "Eric S."),
            ("GBTDATA_TCV", "Tabitha V.")]

for userdiritem in userdirs:
    dbcl.register_envpath(userdiritem[0], userdiritem[1])

# a standard local file path for tests
dbcl.register_path('local_test', './data_test',
                   "user's local directory for running tests")

#-----------------------------------------------------------------------------
# paths to maps 'register_path'
# unique information: mappath directory, mappath key, mappath notes, mappath desc
#-----------------------------------------------------------------------------
pathname = dbcl.fetch_path('GBTDATA_TCV') + "maps/"
notes = '15hr, 22hr, and 1hr maps made with the new calibration'
dbcl.register_path('GBT_maps_Tabitha', pathname,
                    "Tabitha's map directory", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_ESWITZER') + "GBT/maps/"
notes = '15hr, 22hr, and 1hr maps made with the new calibration'
dbcl.register_path('GBT_maps_Eric', pathname,
                    "Eric's map directory", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_ESWITZER') + "GBT/maps/fdgcal_plussim/"
notes = 'fdgcal plus sim'
dbcl.register_path('GBT_maps_Eric_fdgcal_plussim', pathname,
                    "Eric's map directory (fdgcal plus sim)", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_ESWITZER') + "GBT/maps/signal_only/"
notes = 'signal sim only sim'
dbcl.register_path('GBT_maps_Eric_signal_only', pathname,
                    "Eric's map directory: signal only", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_NBANAVAR') + "gbt_out/maps/4_section_maps/old_cal/"
notes = 'calibration test maps: old calibration, old mapper'
dbcl.register_path('GBT_maps_Nidhi_oldcal', pathname,
                    "Nidhi's map directory", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_ESWITZER') + "GBT/maps/oldcal_plussim/"
notes = 'oldcal plus sim'
dbcl.register_path('GBT_maps_Eric_oldcal_plussim', pathname,
                    "Eric's map directory (oldcal plus sim)", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_KIYO') + "gbt_out/maps/"
notes = 'optimal maps, proposal-era cal'
dbcl.register_path('GBT_maps_Kiyo', pathname,
                    "Kiyo's map directory", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_KIYO') + "gbt_out/maps/jan16.2012/"
notes = 'optimal maps, proposal-era cal'
dbcl.register_path('GBT_maps_Kiyo_16jan2012', pathname,
                    "Kiyo's map directory", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_KIYO') + "gbt_out/maps/feb19.2012/"
notes = 'optimal maps, proposal-era cal, some bugs fixed following Jan 16 2012'
dbcl.register_path('GBT_maps_Kiyo_19feb2012', pathname,
                    "Kiyo's map directory", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_KIYO') + "gbt_out/maps/feb21.2012/"
notes = 'optimal maps, proposal-era cal, some bugs fixed following Jan 16 2012, fdgcal'
dbcl.register_path('GBT_maps_Kiyo_21feb2012', pathname,
                    "Kiyo's map directory", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_KIYO') + "gbt_out/maps/feb23.2012/"
notes = 'optimal maps, proposal-era cal, some bugs fixed following Jan 16 2012, fdgcal, TOD modes removed'
dbcl.register_path('GBT_maps_Kiyo_23feb2012', pathname,
                    "Kiyo's map directory", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_KIYO') + "gbt_out/maps/apr11.2012/"
notes = 'map where each session is calibrated wrt session 1 using zero mode'
dbcl.register_path('GBT_maps_Kiyo_11apr2012', pathname,
                    "Kiyo's map directory", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_KIYO') + "gbt_out/maps/apr16.2012/"
notes = 'individual sessions for map self-calibration; including combined map'
dbcl.register_path('GBT_maps_Kiyo_16apr2012', pathname,
                    "Kiyo's map directory", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_KIYO') + "gbt_out_new/maps/may04.2012/"
notes = '15hr optimal mapper, best band, self-calibration'
dbcl.register_path('GBT_maps_Kiyo_04may2012', pathname,
                    "Kiyo's map directory", notes=notes)

# Temporary (until noise_weight is included in the map directory)
pathname = dbcl.fetch_path('GBTDATA_ESWITZER') + "GBT/maps/optimal_15hr_selfcal/"
notes = '15hr optimal mapper, best band, self-calibration'
dbcl.register_path('GBT_maps_Eric_04may2012', pathname,
                    "Eric's version of Kiyo map with noise weight", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_KIYO') + "gbt_out_new/maps/may07.2012/"
notes = '15hr optimal mapper, best band, self-calibration'
dbcl.register_path('GBT_maps_Kiyo_07may2012', pathname,
                    "Kiyo's map directory", notes=notes)

#pathname = dbcl.fetch_path('GBTDATA_ESWITZER') + "GBT/maps/oldcal_plussim10pct/"
#notes = 'oldcal plus sim * 1.1'
#dbcl.register_path('GBT_maps_Eric_oldcal_plussim10pct', pathname,
#                    "Eric's map directory (oldcal plus sim)", notes=notes)

#pathname = dbcl.fetch_path('GBTDATA_NBANAVAR') + "gbt_out/maps/4_section_maps/flux_cal/"
#notes = 'calibration test maps: old calibration, old mapper'
#dbcl.register_path('GBT_maps_Nidhi_fluxcal', pathname,
#                    "Nidhi's map directory", notes=notes)

#-----------------------------------------------------------------------------
# register map data; 'register_maprun'
# unique information: map key, group key, mappath key, map field tag, map desc,
# map notes, map status
#-----------------------------------------------------------------------------
# register the GBT keys and paths
tcv_cal_note1 = "Calibration was done with the point source files that are \
found in `/mnt/raid-project/gmrt/tcv/mueller_params/` with the suffix \
`_mueller_from_inverted_params.txt` These point source files were generated \
using 3-5 sets of onoff scans from one point source (3C286, 3C48, or 3C67). \
For details of which scans were used for each file, check with Tabitha. For \
the 1hr maps, sessions 80-90 on the `10B_036` project and 1-13 on the \
`11B_055` project were used."

tcv_cal_note2 = "Calibration was done with the point source files that are \
found in `/mnt/raid-project/gmrt/tcv/mueller_params/` with the suffix \
`_mueller_from_inverted_params.txt` These point source files were generated \
using 3-5 sets of onoff scans from one point source (3C286, 3C48, or 3C67). \
For details of which scans were used for each file, check with Tabitha. For \
the 15hr and 22hr, sessions 41-90 in `10B_036` were used."

tcv_cal_note3 = "Calibration is done with a new scheme of using XX, YY \
separately from the calibration source to constrain the differential gain, \
then the flux calibration is applied from source observations."

group_key = 'GBTmaps'
status = 'static'

key = 'GBT_15hr_map_fluxpolcal'
parent = 'GBT_maps_Tabitha'
desc = "15hr maps with Oct. 10 2011 flux+pol calibration"
field_tag = '15hr_41-90'
dbcl.register_maprun(key, group_key, parent, field_tag, desc,
                        notes=tcv_cal_note2, status=status)

key = 'GBT_22hr_map_fluxpolcal'
parent = 'GBT_maps_Tabitha'
desc = "22hr maps with Oct. 10 2011 flux+pol calibration"
field_tag = '22hr_41-90'
dbcl.register_maprun(key, group_key, parent, field_tag, desc,
                        notes=tcv_cal_note2, status=status)

key = 'GBT_1hr_map_fluxpolcal'
parent = 'GBT_maps_Tabitha'
desc = "1hr maps with Oct. 10 2011 flux+pol calibration"
field_tag = '1hr_41-16'
dbcl.register_maprun(key, group_key, parent, field_tag, desc,
                        notes=tcv_cal_note1, status=status)

# register Tabitha's new calibration
key = 'GBT_15hr_map_fdgcal'
parent = 'GBT_maps_Tabitha'
desc = "15hr maps with Feb. 20 XXYY+flux calibration"
field_tag = '15hr_41-90_fdg'
status = 'in development'
dbcl.register_maprun(key, group_key, parent, field_tag, desc,
                        notes=tcv_cal_note3, status=status)

# register the new maps plus sim
key = 'GBT_15hr_map_fdgcal_plussim'
parent = 'GBT_maps_Eric_fdgcal_plussim'
desc = "15hr maps with Feb. 20 XXYY+flux calibration plus 15hr_oldmap_str sim_beam_001.npy"
field_tag = '15hr_41-90_fdg'
status = 'in development'
dbcl.register_maprun(key, group_key, parent, field_tag, desc,
                        notes=tcv_cal_note3, status=status)

# register the signal-only
key = 'GBT_15hr_map_signal_only'
parent = 'GBT_maps_Eric_signal_only'
desc = "15hr maps with sim_beam_001.npy in place of the real data"
field_tag = '15hr_41-90_fdg'
status = 'in development'
dbcl.register_maprun(key, group_key, parent, field_tag, desc,
                        notes=tcv_cal_note3, status=status)

#-----------------------------------------------------------------------------
# register optimal maps; 'register_optimalmap_section_run',
# 'register_optimalmap_glued_run'
# unique information: map key, group key, mappath key, map field tag, map desc,
# map notes, map status, map subsection
#-----------------------------------------------------------------------------
group_key = 'GBTmaps'
status = 'in development'

key = 'GBT_15hr_optimalmap_selfcal_762'
parent = 'GBT_maps_Eric_04may2012'
desc = "15hr optimal mapper, best band, self-calibration"
field_tag = '15hr_41-90'
dbcl.register_optimalmap_section_run(key, group_key, parent, field_tag, '762', desc,
                        status=status)

#key = 'GBT_15hr_optimalmap737'
#parent = 'GBT_maps_Kiyo'
#desc = "optimal (section) maps from Kiyo"
#field_tag = '15hr_41-90'
#dbcl.register_optimalmap_section_run(key, group_key, parent, field_tag, '737', desc,
#                        status=status)
#key = 'GBT_15hr_optimalmap799'
#parent = 'GBT_maps_Kiyo'
#desc = "optimal (section) maps from Kiyo"
#field_tag = '15hr_41-90'
#dbcl.register_optimalmap_section_run(key, group_key, parent, field_tag, '799', desc,
#                        status=status)
#key = 'GBT_15hr_optimalmap862'
#parent = 'GBT_maps_Kiyo'
#desc = "optimal (section) maps from Kiyo"
#field_tag = '15hr_41-90'
#dbcl.register_optimalmap_section_run(key, group_key, parent, field_tag, '862', desc,
#                        status=status)

key = 'GBT_15hr_optimalmap_fluxpolcal'
parent = 'GBT_maps_Kiyo_16jan2012'
desc = "optimal maps glued into one cube"
field_tag = '15hr_41-90'
dbcl.register_optimalmap_glued_run(key, group_key, parent, field_tag, desc,
                        status=status)

key = 'GBT_15hr_optimalmap_mapv2oldcal'
parent = 'GBT_maps_Kiyo_19feb2012'
desc = "optimal maps glued into one cube; some bugs fixed following 16Jan2012"
field_tag = '15hr_41-90'
dbcl.register_optimalmap_glued_run(key, group_key, parent, field_tag, desc,
                        status=status)

key = 'GBT_15hr_optimalmap_mapv2fdgcal'
parent = 'GBT_maps_Kiyo_21feb2012'
desc = "optimal maps glued into one cube; some bugs fixed following 16Jan2012, fdg cal"
field_tag = '15hr_41-90'
dbcl.register_optimalmap_glued_run(key, group_key, parent, field_tag, desc,
                        status=status)

key = 'GBT_15hr_optimalmap_mapv2fdgcalmoderm'
parent = 'GBT_maps_Kiyo_23feb2012'
desc = "optimal maps glued into one cube; some bugs fixed following 16Jan2012, fdg cal, modes removed in TOD"
field_tag = '15hr_41-90'
dbcl.register_optimalmap_glued_run(key, group_key, parent, field_tag, desc,
                        status=status)

key = 'GBT_15hr_map_mapcal'
parent = 'GBT_maps_Kiyo_11apr2012'
desc = 'map where each session is calibrated wrt session 1 using zero mode'
field_tag = 'calib_15hr_41-73'
dbcl.register_maprun(key, group_key, parent, field_tag, desc,
                     status=status, sectag="comb", skip_firstpass=True)

key = 'GBT_15hr_map_mapcal2'
parent = 'GBT_maps_Kiyo_16apr2012'
desc = 'map where each session is calibrated wrt stack w factorize noise'
field_tag = 'calib_15hr_41-73'
dbcl.register_maprun(key, group_key, parent, field_tag, desc,
                     status=status, sectag="comb", skip_firstpass=True)

#-----------------------------------------------------------------------------
# register alternate calibrations from Nidhi; 'register_maprun'
# unique information: map key, group key, mappath key, map field tag, map desc,
# map notes, map status
#-----------------------------------------------------------------------------
group_key = 'GBTmaps'
status = 'in development'

key = 'GBT_15hr_map_oldcal'
parent = 'GBT_maps_Nidhi_oldcal'
desc = "maps with old-style calibration for comparison"
field_tag = '15hr_41-90'
dbcl.register_maprun(key, group_key, parent, field_tag, desc,
                     status=status)

# register the oldcal maps plus sim
key = 'GBT_15hr_map_oldcal_plussim'
parent = 'GBT_maps_Eric_oldcal_plussim'
desc = "15hr maps with proposal-era calibration plus 15hr_oldmap_str sim_beam_001.npy"
field_tag = '15hr_41-90'
status = 'in development'
dbcl.register_maprun(key, group_key, parent, field_tag, desc,
                        notes=tcv_cal_note3, status=status)

# register the oldcal maps plus sim
#key = 'GBT_15hr_map_oldcal_plussim10pct'
#parent = 'GBT_maps_Eric_oldcal_plussim10pct'
#desc = "15hr maps with proposal-era calibration plus 15hr_oldmap_str sim_beam_001.npy"
#field_tag = '15hr_41-90'
#status = 'in development'
#dbcl.register_maprun(key, group_key, parent, field_tag, desc,
#                        notes=tcv_cal_note3, status=status)

#key = 'GBT_15hr_map_fluxcal'
#parent = 'GBT_maps_Nidhi_fluxcal'
#desc = "maps with fluxcal only"
#field_tag = '15hr_41-90'
#dbcl.register_maprun(key, group_key, parent, field_tag, desc,
#                     status=status)

#-----------------------------------------------------------------------------
# functions to register the cleaned maps (new style)
# here, the cases are
# map
# map_plussim_alt
# map_plussim_minussim_alt
# map_plussim_minusmap_alt
# where _alt is an additional tag like mult1p1 for simulation multiplied by 1.1
#
# this only need to be modified if the structure of the maps cleaning chain
# is modified. Different maps can be input below.
#-----------------------------------------------------------------------------
def mode_clean_run(source_key, username, modelist,
                   tag="", status=None, notes=None, simkey=None, simtag="",
                   alt="", extdesc=""):
    if simkey:
        key = '%s_cleaned%s%s' % (source_key, simtag, alt)
        if simtag == "_plussim":
            desc = '`%s` cleaned using map `%s`; %s' % (simkey, source_key, extdesc)

        if simtag == "_plussim_minussim":
            desc = '`%s` cleaned using map `%s`, cleaned sim subtracted; %s' % (simkey, source_key, extdesc)

        if simtag == "_plussim_minusmap":
            desc = '`%s` cleaned using map `%s`, cleaned map subtracted; %s' % (simkey, source_key, extdesc)
    else:
        key = '%s_cleaned%s' % (source_key, alt)
        desc = '`%s` cleaned ; %s' % (source_key, extdesc)

    combined_key = '%s_combined' % key
    group_key = 'GBTcleaned'
    parent_key = '%s_path_%s' % (key, username)

    pathname = '%s/GBT/cleaned_maps/%s%s%s/' % \
                (dbcl.fetch_path('GBTDATA_ESWITZER'), source_key, simtag, alt)

    dbcl.register_path(parent_key, pathname, desc, notes=notes)

    dbcl.register_fourway_list(key, group_key, parent_key, desc, modelist,
                 notes=notes, status=status, paramfile="params.ini", tag="",
                 register_modes=True)

    dbcl.register_combined_maprun(combined_key, group_key, parent_key, desc,
                                  modelist, notes=notes, status=status)

def mapsim_mode_clean_run(map_source_key, sim_source_key,
                   username, modelist,
                   tag="", status=None, notes=None,
                   alt="", extdesc=""):

    mode_clean_run(map_source_key, username, modelist,
                   tag=tag, status=status, notes=notes,
                   simkey=sim_source_key,
                   simtag="_plussim",
                   alt=alt, extdesc=extdesc)

    mode_clean_run(map_source_key, username, modelist,
                   tag=tag, status=status, notes=notes,
                   simkey=sim_source_key,
                   simtag="_plussim_minussim",
                   alt=alt, extdesc=extdesc)

    mode_clean_run(map_source_key, username, modelist,
                   tag=tag, status=status, notes=notes,
                   simkey=sim_source_key,
                   simtag="_plussim_minusmap",
                   alt=alt, extdesc=extdesc)

notes = "proposal era calibration and old mapmaker"
extdesc = "the mean is removed, radial modes subtracted (no common res conv.)"
status = 'static'
modelist = range(0, 105, 5)
#    mode_clean_run(map_source_key, username, modelist,
#                   tag=tag, status=status, notes=notes, sim=None,
#                   alt=alt, extdesc=extdesc)

# TODO: uncomment these: ONLY TEMPORARY
#mode_clean_run('GBT_15hr_map_oldcal', 'Eric', modelist,
#               status=status, notes=notes, simkey=None,
#               alt="", extdesc=extdesc)

#mapsim_mode_clean_run('GBT_15hr_map_oldcal', 'sim_15hr_oldmap_str_beam',
#                      "Eric", modelist, status=status, notes=notes,
#                      extdesc=extdesc)

#mapsim_mode_clean_run('GBT_15hr_map_oldcal', 'sim_15hr_oldmap_str_beam',
#                      "Eric", modelist, status=status, notes=notes,
#                      extdesc=extdesc, alt="_simx1p1")

notes = "map calibration round 1: old pixelization, 1x1, no factorization"
status = 'development'
mode_clean_run('GBT_15hr_map_mapcal', 'Eric', modelist,
               status=status, notes=notes, simkey=None,
               alt="", extdesc=extdesc)

mapsim_mode_clean_run('GBT_15hr_map_mapcal', 'sim_15hr_oldmap_str_beam',
                      "Eric", modelist, status=status, notes=notes,
                      extdesc=extdesc)

mapsim_mode_clean_run('GBT_15hr_map_mapcal', 'sim_15hr_oldmap_str_beam',
                      "Eric", modelist, status=status, notes=notes,
                      extdesc=extdesc, alt="_simx1p1")

notes = "map calibration round 1: old pixelization, 1xcomb, factorization"
status = 'development'
mode_clean_run('GBT_15hr_map_mapcal2', 'Eric', modelist,
               status=status, notes=notes, simkey=None,
               alt="", extdesc=extdesc)

mapsim_mode_clean_run('GBT_15hr_map_mapcal2', 'sim_15hr_oldmap_str_beam',
                      "Eric", modelist, status=status, notes=notes,
                      extdesc=extdesc)

mapsim_mode_clean_run('GBT_15hr_map_mapcal2', 'sim_15hr_oldmap_str_beam',
                      "Eric", modelist, status=status, notes=notes,
                      extdesc=extdesc, alt="_simx1p1")

#-----------------------------------------------------------------------------
# functions to register the cleaned maps
# this only need to be modified if the structure of the maps cleaning chain
# is modified. Different maps can be input below.
#-----------------------------------------------------------------------------
def mode_clean_run_old(source_key, username, modelist,
                   tag="", status=None, notes=None, sim=None,
                   alt="", extdesc=""):
    mapsim = ""
    if sim:
        mapsim = "_sims"

    # TODO: REVERT THIS TO OLDCLEAN
    key = '%s_cleaned%s%s' % (source_key, mapsim, alt)
    #key = '%s_oldcleaned%s%s' % (source_key, mapsim, alt)
    combined_key = '%s_combined' % key
    group_key = 'GBTcleaned'
    parent_key = '%s_path_%s' % (key, username)

    pathname = '%s/GBT/cleaned_maps/%s%s%s/' % \
                (dbcl.fetch_path('GBTDATA_ESWITZER'), source_key, mapsim, alt)
    if sim:
        desc = '`%s` (map index 000) cleaned using map `%s`; %s' % (sim, source_key, extdesc)
    else:
        desc = '`%s` cleaned ; %s' % (source_key, extdesc)

    dbcl.register_path(parent_key, pathname, desc, notes=notes)

    dbcl.register_fourway_list(key, group_key, parent_key, desc, modelist,
                 notes=notes, status=status, paramfile="params.ini", tag="",
                 register_modes=True)

    dbcl.register_combined_maprun(combined_key, group_key, parent_key, desc,
                                  modelist, notes=notes, status=status)

def mapsim_mode_clean_run_old(map_source_key, sim_source_key,
                   username, modelist,
                   tag="", status=None, notes=None,
                   alt="", extdesc=""):

    mode_clean_run_old(map_source_key, username, modelist,
                   tag=tag, status=status, notes=notes, sim=None,
                   alt=alt, extdesc=extdesc)

    mode_clean_run_old(map_source_key, username, modelist,
                   tag=tag, status=status, notes=notes,
                    sim=sim_source_key,
                   alt=alt, extdesc=extdesc)

# high level call to register a mode cleaning run + sims for conv, noconv case
def mapsimnoconv_mode_clean_run(map_source_key, sim_source_key, username="Eric",
                                tag="", status=None, notes=None, modelist=None):
    # For all maps, we clean 0, 5 ... 50 modes
    if modelist is None:
        modelist = range(0, 105, 5)

    extdesc = "the mean is removed, convolved to a common beam, radial modes subtracted"
    mapsim_mode_clean_run_old(map_source_key, sim_source_key, 'Eric', modelist,
                          status=status, notes=notes, extdesc=extdesc)

    extdesc = "the mean is removed, radial modes subtracted"
    mapsim_mode_clean_run_old(map_source_key, sim_source_key, 'Eric', modelist,
                          status=status, notes=notes, alt="_noconv", extdesc=extdesc)

#-----------------------------------------------------------------------------
# calls to register various complete cleaned map + sim runs
#-----------------------------------------------------------------------------
notes = "Tabitha flux+pol calibration and old mapmaker"
status = 'static'
mapsimnoconv_mode_clean_run('GBT_15hr_map_fluxpolcal', 'sim_15hr_oldmap_str_beam',
                             status=status, notes=notes)

notes = "Tabitha flux+pol calibration and old mapmaker"
status = 'static'
mapsimnoconv_mode_clean_run('GBT_22hr_map_fluxpolcal', 'sim_22hr_oldmap_str_beam',
                             status=status, notes=notes)

notes = "Tabitha flux+pol calibration and old mapmaker"
status = 'static'
mapsimnoconv_mode_clean_run('GBT_1hr_map_fluxpolcal', 'sim_1hr_oldmap_str_beam',
                             status=status, notes=notes)

notes = "proposal era calibration and old mapmaker"
status = 'static'
mapsimnoconv_mode_clean_run('GBT_15hr_map_oldcal', 'sim_15hr_oldmap_str_beam',
                             status=status, notes=notes)

notes = "Tabitha XX,YY + flux calibration calibration and old mapmaker"
status = 'active development'
mapsimnoconv_mode_clean_run('GBT_15hr_map_fdgcal', 'sim_15hr_oldmap_str_beam',
                             status=status, notes=notes)

notes = "Tabitha XX,YY + flux calibration calibration and old mapmaker plus sim"
status = 'active development'
mapsimnoconv_mode_clean_run('GBT_15hr_map_fdgcal_plussim', 'sim_15hr_oldmap_str_beam',
                             status=status, notes=notes)

notes = "Tabitha XX,YY + flux calibration calibration and old mapmaker cleaned with modes from the same map as in plussim"
status = 'active development'
mapsimnoconv_mode_clean_run('GBT_15hr_map_fdgcal_cleanedplussim', 'sim_15hr_oldmap_str_beam',
                             status=status, notes=notes)

notes = "proposal-era calibration calibration and old mapmaker plus sim"
status = 'active development'
mapsimnoconv_mode_clean_run('GBT_15hr_map_oldcal_plussim', 'sim_15hr_oldmap_str_beam',
                             status=status, notes=notes)

notes = "proposal-era calibration calibration and old mapmaker cleaned with modes from the same map as in plussim"
status = 'active development'
mapsimnoconv_mode_clean_run('GBT_15hr_map_oldcal_cleanedplussim', 'sim_15hr_oldmap_str_beam',
                             status=status, notes=notes)

notes = "signal-only"
status = 'active development'
mapsimnoconv_mode_clean_run('GBT_15hr_map_signal_only', 'sim_15hr_oldmap_str_beam',
                             status=status, notes=notes)

notes = "Tabitha flux+pol calibration calibration and new mapmaker"
status = 'active development'
mapsimnoconv_mode_clean_run('GBT_15hr_optimalmap_fluxpolcal', 'sim_15hr_optimalmap_str_beam',
                             status=status, notes=notes)

notes = "prop. era calibration calibration and new mapmaker, v2"
status = 'active development'
mapsimnoconv_mode_clean_run('GBT_15hr_optimalmap_mapv2oldcal', 'sim_15hr_optimalmap_str_beam',
                             status=status, notes=notes)

notes = "Tabitha fdg calibration calibration and new mapmaker, v2"
status = 'active development'
mapsimnoconv_mode_clean_run('GBT_15hr_optimalmap_mapv2fdgcal', 'sim_15hr_optimalmap_str_beam',
                             status=status, notes=notes)

notes = "Tabitha fdg calibration calibration and new mapmaker v2 with time modes removed"
status = 'active development'
mapsimnoconv_mode_clean_run('GBT_15hr_optimalmap_mapv2fdgcalmoderm', 'sim_15hr_optimalmap_str_beam',
                             status=status, notes=notes)

notes = "Kiyo optimal maps, map-domain calibration, funny border"
status = 'active development'
mapsimnoconv_mode_clean_run('GBT_15hr_optimalmap_selfcal_762', 'sim_15hr_optimalmap762_str_beam',
                             status=status, notes=notes)

#notes = "Tabitha flux-only calibration calibration and old mapmaker"
#status = 'UNRESOLVED FACTOR OF 2 ERROR'
#mapsimnoconv_mode_clean_run('GBT_15hr_map_fluxcal', 'sim_15hr_oldmap_str_beam',
#                             status=status, notes=notes)

#notes = "proposal-era calibration calibration and old mapmaker plus sim*1.1"
#status = 'active development'
#mapsimnoconv_mode_clean_run('GBT_15hr_map_oldcal_plussim10pct', 'sim_15hr_oldmap_str_beam',
#                             status=status, notes=notes)

#-----------------------------------------------------------------------------
# batch cleaning of simulations for a transfer function
#-----------------------------------------------------------------------------
pathname = dbcl.fetch_path('GBTDATA_ESWITZER') + "map_cleaning_cache/"
notes = 'map cleaning cache for batch runs'
dbcl.register_path('GBT_15hr_map_batch_sim_depot', pathname,
                    "Map cleaning cache", notes=notes)

key = 'batch_cleaned_cache'
group_key = 'GBTcleaned'
parent_key = 'GBT_15hr_map_batch_sim_depot'
desc = 'cache of files for cleaned maps (intermediate product)'
modelist = range(0, 105, 5)
simlist = range(0, 100)
notes = "this is temporary and these files are continuously overwritten"
status = "in development"
dbcl.register_fourway_list(key, group_key, parent_key, desc, modelist,
                 notes=notes, status=status, paramfile="params.ini", tag="",
                 register_modes=True)

parent_key = "GBT_15hr_oldcal_corrfg_path"
pathname = '%s/GBT/cleaned_maps/GBT_15hr_oldcal_corrfg/' % \
                dbcl.fetch_path('GBTDATA_ESWITZER')

dbcl.register_path(parent_key, pathname, desc, notes=notes)

desc = 'clean_{map+sim} (map+sim)'
notes = "cleaned foregrounds are correlated with the signal"
status = "in development"
dbcl.register_combined_batchsimrun("GBT_15hr_oldcal_corrfg", group_key,
                                   parent_key, desc,
                                   modelist, simlist, notes=notes, status=status)

#-----------------------------------------------------------------------------
# paths to WiggleZ data
#-----------------------------------------------------------------------------
def wigglez_paths(fieldname, rootdir, type):
    key = "WiggleZ_%s_%s_path" % (fieldname, type)
    desc = "path to %s wigglez `%s` data" % (fieldname, type)
    notes = "generated by map/optical_catalog.py"
    status = "in development"
    pathname = rootdir + fieldname + "/"
    #pathname = rootdir + "/"
    dbcl.register_path(key, pathname,
                    desc, notes=notes)

for fielditem in field_list:
    rootdir = dbcl.fetch_path('GBTDATA_ESWITZER') + 'wiggleZ/binned/'
    wigglez_paths(fielditem, rootdir, "binned")
    rootdir = dbcl.fetch_path('GBTDATA_ESWITZER') + 'wiggleZ/complete_binned/'
    wigglez_paths(fielditem, rootdir, "complete_binned")
    rootdir = dbcl.fetch_path('GBTDATA_ESWITZER') + 'wiggleZ/binned_delta/'
    wigglez_paths(fielditem, rootdir, "delta_binned")
    rootdir = dbcl.fetch_path('GBTDATA_ESWITZER') + 'wiggleZ/complete_binned_delta/'
    wigglez_paths(fielditem, rootdir, "complete_delta_binned")
    rootdir = dbcl.fetch_path('GBTDATA_ESWITZER') + 'wiggleZ/catalogs/'
    wigglez_paths(fielditem, rootdir, "catalog")

#-----------------------------------------------------------------------------
# register WiggleZ data sets
#-----------------------------------------------------------------------------
group_key = 'WiggleZ'
notes = 'catalog from C. Blake'
status = 'static'

for fielditem in field_list:
    parent = 'WiggleZ_%s_catalog_path' % fielditem
    hrbase = int(fielditem.split("hr")[0])

    key = 'WiggleZ_%s_catalog_data' % fielditem
    filename = 'reg%02ddata.dat' % hrbase
    desc = 'WiggleZ %s field catalog' % fielditem
    dbcl.register_file(key, group_key, parent, filename, desc,
                       notes=notes, status=status)

    key = 'WiggleZ_%s_mock_catalog' % fielditem
    prefix = "reg%02drand" % hrbase
    indices = range(0, 1000)
    desc = '%s WiggleZ mock catalog' % fielditem
    dbcl.register_file_set(key, group_key, parent, prefix, indices, desc,
                  notes=notes, status=status, suffix=".dat")

    key = 'WiggleZ_%s_priority_table' % fielditem
    filename = 'nzpri_reg%02d_tzuchingcats.dat' % hrbase
    desc = 'WiggleZ %s random catalog priority table' % fielditem
    dbcl.register_file(key, group_key, parent, filename, desc,
                       notes=notes, status=status)

def register_wigglez_products(fieldname, complete=False, delta=False):
    tdesc = ""
    if complete:
        ctag = "complete_"
        tdesc += " (complete survey region)"
    else:
        ctag = ""

    if delta:
        dtag = "delta_"
        tdesc += " (overdensity from sel. func)"
    else:
        dtag = ""

    notes = 'generated by `map/optical_catalog.py`; catalog C. Blake'
    status = 'in development'
    parent = 'WiggleZ_%s_%s%sbinned_path' % (fieldname, ctag, dtag)
    hrbase = int(fieldname.split("hr")[0])

    key = 'WiggleZ_%s_%s%sbinned_data' % (fieldname, ctag, dtag)
    filename = 'reg%02ddata.npy' % hrbase
    desc = 'binned %s WiggleZ data%s' % (fieldname, tdesc)
    dbcl.register_file(key, group_key, parent, filename, desc,
                       notes=notes, status=status)

    key = 'WiggleZ_%s_%s%smock' % (fieldname, ctag, dtag)
    prefix = "reg%02drand" % hrbase
    indices = range(0, 100)
    desc = '%s WiggleZ mock binned catalog%s' % (fieldname, tdesc)
    dbcl.register_file_set(key, group_key, parent, prefix, indices, desc,
                  notes=notes, status=status)

    if not delta:
        key = 'WiggleZ_%s_%s%sselection' % (fieldname, ctag, dtag)
        filename = 'reg%02dselection.npy' % hrbase
        desc = '%s WiggleZ data selection function, 1000 catalogs%s' % (fieldname, tdesc)
        dbcl.register_file(key, group_key, parent, filename, desc,
                           notes=notes, status=status)

        key = 'WiggleZ_%s_%s%smontecarlo' % (fieldname, ctag, dtag)
        filename = 'reg%02dmontecarlo.npy' % hrbase
        desc = '%s WiggleZ data selection function, Monte Carlo%s' % (fieldname, tdesc)
        dbcl.register_file(key, group_key, parent, filename, desc,
                           notes=notes, status=status)

        key = 'WiggleZ_%s_%s%sseparable_selection' % (fieldname, ctag, dtag)
        filename = 'reg%02dseparable.npy' % hrbase
        desc = '%s WiggleZ data selection function, separability%s' % (fieldname, tdesc)
        dbcl.register_file(key, group_key, parent, filename, desc,
                           notes=notes, status=status)

for fielditem in field_list:
    register_wigglez_products(fielditem)
    register_wigglez_products(fielditem, complete=True)
    register_wigglez_products(fielditem, delta=True)
    register_wigglez_products(fielditem, complete=True, delta=True)

#-----------------------------------------------------------------------------
# register the simulations
#-----------------------------------------------------------------------------
# register the various preparations of the simulation
def register_sim(fieldname, tag="", basedesc="sim.; WiggleZ pwrspec",
                 notes=None, status=None):
    # start by registering the parent directory
    group_key = "Simulations"
    parent = "sim_%s_%s_path" % (fieldname, tag)
    rootdir = dbcl.fetch_path('GBTDATA_ESWITZER') + 'GBT/simulations/'
    pathname = "%s%s_%s/" % (rootdir, fieldname, tag)
    pathdesc = "path to all %s %s simulations" % (fieldname, tag)
    pathnotes = 'generated by ERS using JRS code + wigglez'
    dbcl.register_path(parent, pathname, pathdesc, notes=pathnotes)

    # register 100 simulations per case
    indices = range(0, 100)

    simcases = {"_physical": "physical coordinates",
                "_temperature": "raw 21cm temperature map",
                "_delta": "delta overdensity map (21cm/T_b)",
                "_beam": "beam-convolved 21cm temperature map",
                "_beam_meansub": "beam-convolved 21cm temperature map, mean subtracted",
                "_beam_meansubconv": "beam-convolved 21cm temperature map, common-res. conv., mean subtracted",
                "_beam_conv": "beam-convolved 21cm temperature map, common-res. conv.",
                "_beam_plus_data": "beam-convolved 21cm temperature map plus real data",
                "_beam_plus_fg": "beam-convolved 21cm temperature map plus fg model"
               }

    for prep in simcases:
        prefix = "sim%s_" % prep
        key = "sim_%s_%s%s" % (fieldname, tag, prep)
        desc = "%s %s; %s" % (fieldname, basedesc, simcases[prep])
        dbcl.register_file_set(key, group_key, parent, prefix, indices, desc,
                  notes=notes, status=status)

# register sims of the various preparations above for three power spectra
def register_simset(treatment, fielditem, treatmentdesc):
    '''For each of the map treatments, register three classes of simulations:
    an ideal case with only P_dd, then P_ddvv and P_ddvvstreaming
    '''
    status = "code is stable"
    basedesc = "sim.; WiggleZ pwrspec (no mean, no evolution, density only)"
    basedesc += "; " + treatmentdesc
    register_sim(fielditem, status=status, basedesc=basedesc,
                 tag=treatment + "_ideal")

    basedesc = "sim.; WiggleZ pwrspec with dd and vv but no streaming"
    basedesc += "; " + treatmentdesc
    register_sim(fielditem, status=status, basedesc=basedesc,
                 tag=treatment + "_nostr")

    basedesc = "sim.; WiggleZ pwrspec standard + streaming velocities"
    basedesc += "; " + treatmentdesc
    register_sim(fielditem, status=status, basedesc=basedesc,
                 tag=treatment + "_str")

def register_strsim(treatment, fielditem, treatmentdesc):
    '''For each of the map treatments, register three classes of simulations:
    an ideal case with only P_dd, then P_ddvv and P_ddvvstreaming
    '''
    status = "code is stable"
    basedesc = "sim.; WiggleZ pwrspec standard + streaming velocities"
    basedesc += "; " + treatmentdesc
    register_sim(fielditem, status=status, basedesc=basedesc,
                 tag=treatment + "_str")

# highest level of registering simulations by field/pixelization and omega_HI
register_simset("oldmap", "15hr", "old-style map pixelization")
register_strsim("oldmap", "22hr", "old-style map pixelization")
register_strsim("oldmap", "1hr", "old-style map pixelization")
register_strsim("oldmap_HI5em4", "15hr", "old-style map pixelization, Omega_HI=5e-4")
register_strsim("optimalmap", "15hr", "optimal map pixelization")
register_strsim("optimalmap762", "15hr", "optimal map pixelization, best region")

#-----------------------------------------------------------------------------
# register quadratic product output and cache directories
#-----------------------------------------------------------------------------
pathname = dbcl.fetch_path('GBTDATA_ESWITZER') + "quadratic_products/data/"
notes = 'quadratic products of the data (batch runs)'
dbcl.register_path('quadratic_batch_data', pathname,
                    "Quadratic batch runs of the data", notes=notes)

#-----------------------------------------------------------------------------
# finalize and write out the database
#-----------------------------------------------------------------------------
db_out = shelve.open("path_database.shelve", 'n')
db_out['groups'] = dbcl.groups
db_out['group_order'] = dbcl.grouplist
db_out['pdb'] = dbcl.pdb
db_out['version_tag'] = version_tag
db_out.close()
