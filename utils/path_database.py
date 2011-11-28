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
version_tag = "October 31st, 2011"
pdb = {}
groups = {}
dbcl = path_forms.PathForms(verbose=False)
field_list = ['15hr', '22hr', '1hr']

dbcl.register_list_empty_groups(['GBTmaps', 'Correlation_products',
                           'GBTcleaned', 'WiggleZ', 'Simulations'])

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
# paths to maps
#-----------------------------------------------------------------------------
pathname = dbcl.fetch_path('GBTDATA_TCV') + "maps/"
notes = '15hr, 22hr, and 1hr maps made with the new calibration'
dbcl.register_path('GBT_maps_Tabitha', pathname,
                    "Tabitha's map directory", notes=notes)

pathname = dbcl.fetch_path('GBTDATA_KIYO') + 'wiggleZ/maps/'
notes = '15hr maps from the proposal era'
dbcl.register_path('GBT_maps_Kiyo', pathname,
                    "Kiyo's map directory", notes=notes)

#-----------------------------------------------------------------------------
# register map data
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

group_key = 'GBTmaps'
status = 'in development'

key = 'GBT_15hr_map'
parent = 'GBT_maps_Tabitha'
desc = "15hr maps with Oct. 10 2011 calibration"
field_tag = '15hr_41-90'
dbcl.register_maprun(key, group_key, parent, field_tag, desc,
                        notes=tcv_cal_note2, status=status)
key = 'GBT_22hr_map'
parent = 'GBT_maps_Tabitha'
desc = "22hr maps with Oct. 10 2011 calibration"
field_tag = '22hr_41-90'
dbcl.register_maprun(key, group_key, parent, field_tag, desc,
                        notes=tcv_cal_note2, status=status)
key = 'GBT_1hr_map'
parent = 'GBT_maps_Tabitha'
desc = "1hr maps with Oct. 10 2011 calibration"
#field_tag = '1hr_none'
field_tag = '1hr_41-16'
dbcl.register_maprun(key, group_key, parent, field_tag, desc,
                        notes=tcv_cal_note1, status=status)

# older maps
key = 'GBT_15hr_map_proposal'
parent = 'GBT_maps_Kiyo'
desc = "15 hr maps with proposal-era calibrations"
field_tag = '15hr_41-73'
status = 'in development'
notes = 'this dataset is broken because the noise_inv was over-written \
by something else at a later time'
dbcl.register_maprun(key, group_key, parent, field_tag, desc,
                        notes=notes, status=status)

#-----------------------------------------------------------------------------
# paths to cleaned maps
#-----------------------------------------------------------------------------
for fielditem in field_list:
    pathname = '%s/GBT/cleaned_maps/%s/' % \
                (dbcl.fetch_path('GBTDATA_ESWITZER'), fielditem)
    notes = 'active development stream'
    key = 'GBT_cleaned_%s_maps_Eric' % fielditem
    desc = 'ERS cleaned %s maps' % fielditem
    dbcl.register_path(key, pathname, desc, notes=notes)

# additional trials
pathname = dbcl.fetch_path('GBTDATA_ESWITZER') + \
           'GBT/cleaned_maps/15hr_nomeansub/'
notes = 'active development stream'
desc = 'cleaned 15hr maps without mean subtraction'
dbcl.register_path('GBT_cleaned_nomeansub_15hr_maps_Eric', pathname,
                    desc, notes=notes)

# proposal era
pathname = dbcl.fetch_path('GBTDATA_CALINLIV') + 'wiggleZ/corr/test/'
notes = 'made for the Aug. 1 proposal (Liviu)'
dbcl.register_path('GBT_15hr_Liviu', pathname,
                    "15hr maps with 15 modes removed", notes=notes)

#-----------------------------------------------------------------------------
# register the cleaned maps
#-----------------------------------------------------------------------------
def mode_clean_run(source_key, fieldname, mode_num, username,
                   tag="", status=None, notes=None):
    key = '%s_cleaned_%s%smode' % (source_key, tag, mode_num)
    group_key = 'GBTcleaned'
    parent = 'GBT_cleaned_%s_maps_%s' % (fieldname, username)
    desc = "`%s` maps with %s modes removed" % (source_key, mode_num)
    dbcl.register_fourway_list(key, group_key, parent, desc,
                 notes=None, status=status, paramfile="params.ini", tag="",
                 modenum=mode_num, register_modes=True,
                 register_pickles=False)


for fieldname in field_list:
    source_key = 'GBT_%s_map' % fieldname
    status = 'active development'
    for mode_num in range(0, 55, 5):
        mode_clean_run(source_key, fieldname, mode_num, 'Eric',
                       status=status, notes=notes)
# nomeansub

# proposal era
key = 'GBT_15hr_cleaned_Liviu_15mode'
group_key = 'GBTcleaned'
parent = 'GBT_15hr_Liviu'
status = 'frozen'
notes = 'made for the Aug. 1 proposal (Liviu)'
desc = '15hr maps with 15 modes removed'
tag = "15hr_41-73_"
dbcl.register_fourway_list(key, group_key, parent, desc,
                 notes=notes, status=status, paramfile=None, tag=tag,
                 register_pickles=False, register_modes=False, modenum=None,
                 register_corrsvd=False)


#-----------------------------------------------------------------------------
# paths to combined maps
#-----------------------------------------------------------------------------
for fielditem in field_list:
    pathname = '%s/GBT/combined_maps/%s/' % \
                (dbcl.fetch_path('GBTDATA_ESWITZER'), fielditem)
    notes = 'active development stream'
    key = 'GBT_combined_%s_maps_Eric' % fielditem
    desc = 'ERS combined %s maps' % fielditem
    dbcl.register_path(key, pathname, desc, notes=notes)

#-----------------------------------------------------------------------------
# register the combined datasets
#-----------------------------------------------------------------------------
def register_moderemoved_run(field, num_modes):
    key = 'GBT_%s_combined_cleaned_%smode_map' % (field, modenum)
    group_key = 'GBTcleaned'
    parent = 'GBT_combined_%s_maps_Eric' % field
    filename = 'GBT_%s_map_cleaned_%smode_signal.npy' % (field, modenum)
    desc = '%s field data with %d modes removed, combined' % (field, modenum)
    notes = 'active development stream'
    status = "in development"
    dbcl.register_file(key, group_key, parent, filename, desc,
                       notes=notes, status=status)

    key = 'GBT_%s_combined_cleaned_%smode_product' % (field, modenum)
    filename = 'GBT_%s_map_cleaned_%smode_product.npy' % (field, modenum)
    desc = 'map times weights for %s field data with %d modes removed, combined' % \
           (field, modenum)
    dbcl.register_file(key, group_key, parent, filename, desc,
                       notes=notes, status=status)

    key = 'GBT_%s_combined_cleaned_%smode_weight' % (field, modenum)
    filename = 'GBT_%s_map_cleaned_%smode_weight.npy' % (field, modenum)
    desc = 'weights for %s field data with %d modes removed, combined' % \
           (field, modenum)
    dbcl.register_file(key, group_key, parent, filename, desc,
                       notes=notes, status=status)

for fielditem in field_list:
    for modenum in range(0, 55, 5):
        register_moderemoved_run(fielditem, modenum)

#-----------------------------------------------------------------------------
# paths to WiggleZ data
#-----------------------------------------------------------------------------
# TODO: change desc to reflect binned, etc.
def wigglez_paths(fieldname, rootdir, type):
    key = "WiggleZ_%s_%s" % (fieldname, type)
    desc = "path to all %s wigglez data" % fieldname
    notes = "generated by map/optical_catalog.py"
    status = "in development"
    #pathname = rootdir + fieldname + "/"
    pathname = rootdir + "/"
    dbcl.register_path(key, pathname,
                    desc, notes=notes)

for fielditem in field_list:
    rootdir = dbcl.fetch_path('GBTDATA_ESWITZER') + 'wiggleZ/binned_wiggleZ/'
    wigglez_paths(fielditem, rootdir, "binned_path")
    rootdir = dbcl.fetch_path('GBTDATA_ESWITZER') + 'wiggleZ/wiggleZ_catalogs/'
    wigglez_paths(fielditem, rootdir, "catalog_path")

#-----------------------------------------------------------------------------
# register WiggleZ data sets
#-----------------------------------------------------------------------------
notes = 'generated by `map/optical_catalog.py`; catalog C. Blake'
status = 'in development'
group_key = 'WiggleZ'

for fielditem in field_list:
    parent = 'WiggleZ_%s_binned_path' % fielditem
    hrbase = fielditem.split("hr")[0]

    key = 'WiggleZ_%s_binned_data' % fielditem
    filename = 'reg%sdata.npy' % hrbase
    desc = 'binned %s WiggleZ data' % fielditem
    dbcl.register_file(key, group_key, parent, filename, desc,
                       notes=notes, status=status)

    key = 'WiggleZ_%s_selection' % fielditem
    filename = 'reg%smontecarlo.npy' % hrbase
    desc = '%s WiggleZ data selection function (Monte Carlo)' % fielditem
    dbcl.register_file(key, group_key, parent, filename, desc,
                       notes=notes, status=status)

    key = 'WiggleZ_%s_separable_selection' % fielditem
    filename = 'reg%sseparable.npy' % hrbase
    desc = '%s WiggleZ data selection function (separability)' % fielditem
    dbcl.register_file(key, group_key, parent, filename, desc,
                       notes=notes, status=status)

    key = 'WiggleZ_%s_mock' % fielditem
    prefix = "reg%srand" % hrbase
    indices = range(0, 100)
    desc = '%s WiggleZ mock binned catalog' % fielditem
    dbcl.register_file_set(key, group_key, parent, prefix, indices, desc,
                  notes=notes, status=status)


#-----------------------------------------------------------------------------
# paths to simulations
#-----------------------------------------------------------------------------
def sim_path(fieldname, rootdir, type):
    key = "sim_%s_%s" % (fieldname, type)
    desc = "path to all %s simulations" % fieldname
    notes = 'generated by ERS using JRS code + wigglez'
    status = "in development"
    pathname = rootdir + fieldname + "/"
    dbcl.register_path(key, pathname, desc, notes=notes)

for fielditem in field_list:
    rootdir = dbcl.fetch_path('GBTDATA_ESWITZER') + 'GBT/simulations/'
    sim_path(fielditem, rootdir, "path")


#-----------------------------------------------------------------------------
# register simulations
#-----------------------------------------------------------------------------
# file series: 100 15hr simulations
def register_sim(fieldname, tag="", basedesc="sim.; WiggleZ pwrspec",
                 notes=None, status=None):
    parent = "sim_%s_path" % fieldname
    group_key = "Simulations"
    indices = range(0, 100)

    prefix = "sim%s_" % tag
    key = "sim%s_%s" % (tag, fieldname)
    desc = "%s %s" % (fieldname, basedesc)
    dbcl.register_file_set(key, group_key, parent, prefix, indices, desc,
                  notes=notes, status=status)

    prefix = "sim%s_beam_" % tag
    key = "sim%s_%s_beam" % (tag, fieldname)
    desc = "%s %s; beam convolved" % (fieldname, basedesc)
    dbcl.register_file_set(key, group_key, parent, prefix, indices, desc,
                  notes=notes, status=status)

    prefix = "sim%s_beam_plus_data_" % tag
    key = "sim%s_%s_beam_plus_data" % (tag, fieldname)
    desc = "%s %s; beam convolved plus data" % (fieldname, basedesc)
    dbcl.register_file_set(key, group_key, parent, prefix, indices, desc,
                  notes=notes, status=status)

for fielditem in field_list:
    status = "in development"
    register_sim(fielditem, status=status)

for fielditem in field_list:
    status = "in development"
    basedesc = "sim.; WiggleZ pwrspec with streaming velocities"
    register_sim(fielditem, status=status, basedesc=basedesc, tag="vel")

# register the intermediate correlation product directories
#groups['Paths'].extend(['intermediate_products'])
#groups['Correlation_products'].extend(['GBT_15hr_intemediate_correlations'])
#groups['Correlation_products'].extend(['GBT_22hr_intemediate_correlations'])
#groups['Correlation_products'].extend(['GBT_1hr_intemediate_correlations'])
# foreground_correlations (6 corrs), SVD_set (6 SVDs)
# final correlation products (6 corrs) for fixed # modes removed

db_out = shelve.open("path_database.shelve", 'n')
db_out['groups'] = dbcl.groups
db_out['group_order'] = dbcl.grouplist
db_out['pdb'] = dbcl.pdb
db_out['version_tag'] = version_tag
db_out.close()
