"""Functions for getting data file names."""

import os
import re

_data_dir = os.getenv('GBT10B_DATA')
_kiyo_data_dir = os.getenv('GBT10B_KM')+'data/'
# This line is cheating, will break for other people.
_guppi_data_dir = os.getenv('GBT10B_KM') + 'guppi_data/'

def get_data_files(session_list, field, type=None) :
    """Gets a list of the data file names for each session and field.

    Can pass a list of sessions (as integers) but only a single field.
    """

    out_files = []
    if type is None :
        all_files = os.listdir(_data_dir)
        for session in session_list :
            match_str = ('(%02d'%session + '_wigglez' + field + 
                           "_azel_.*\.raw.acs.fits)")
            for file_name in all_files :
                # See if the file matches.
                if re.match(match_str, file_name) :
                    # Chop off extension and append to list to be returned.
                    root = file_name.split('.')[0]
                    out_files.append(root)
        # Do the same for kiyo data dir.
    elif type == 'kiyo' :
        all_files = os.listdir(_kiyo_data_dir)
        for session in session_list :
            match_str = ('(%02d'%session + '_wigglez' + field + 
                           "_azel_.*\.raw.acs.fits)")
            for file_name in all_files :
                # See if the file matches.
                if re.match(match_str, file_name) :
                    # Chop off extension and append to list to be returned.
                    root = file_name.split('.')[0]
                    out_files.append(root)
        # And finally for guppi data.
    elif type == 'guppi' :
        all_files = os.listdir(_guppi_data_dir)
        for session in session_list :
            match_str = ('(%02d'%session + '_wigglez' + field + 
                           "st_ralongmap_.*\.fits)")
            for file_name in all_files :
                # See if the file matches.
                if re.match(match_str, file_name) :
                    # Chop off extension and append to list to be returned.
                    root = file_name.split('.')[0]
                    out_files.append(root)

    return out_files

def get_cal_files(session_list) :
    """Given a list of data files, returns a list of lists of cal files for
    each data file.
    """
    all_files = os.listdir(_data_dir)
    cal_files = []

    for session in session_list:
        match_str = ('(%02d'%session+'_3c'+'.*_onoff_'+".*\.raw.acs.fits)")
        for file_name in all_files:
            if re.match(match_str, file_name) :
                root = file_name.split('.')[0]
                cal_files.append(root)

    return cal_files

def get_kiyo_data_files(session_list, field) :
    """Gets a list of the data file names for each session and field.

    Can pass a list of sessions (as integers) but only a single field.
    """

    all_files_kiyo = os.listdir(_kiyo_data_dir)
    out_files_kiyo = []
    for session in session_list :
        match_str = ('(%02d'%session + '_wigglez' + field +
                       "_azel_stepping_topo_.*\.raw.acs.fits)")
        for file_name in all_files_kiyo :
            # See if the file matches.
            if re.match(match_str, file_name) :
                # Chop off extension and append to list to be returned.
                root = file_name.split('.')[0] 
                out_files_kiyo.append(root)

    return out_files_kiyo


