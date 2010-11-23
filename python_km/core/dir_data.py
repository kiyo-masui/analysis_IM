"""Functions for getting data file names."""

import os
import re

_data_dir = os.getenv('GBT10B_DATA')

def get_data_files(session_list, field) :
    """Gets a list of the data file names for each session and field.

    Can pass a list of sessions (as integers) but only a single field.
    """

    all_files = os.listdir(_data_dir)
    out_files = []
    for session in session_list :
        match_str = ('(%02d'%session + '_wigglez' + field + 
                       "_azel_.*\.raw.acs.fits)")
        for file_name in all_files :
            # See if the file matches.
            if re.match(match_str, file_name) :
                # Chop off extension and append to list to be returned.
                root = file_name.split('.')[0]
                out_files.append(root)

    return out_files


def get_cal_files(data_file_names) :
    """Given a list of data files, returns a list of lists of cal files for
    each data file.

    Not written yet.
    """
    pass
