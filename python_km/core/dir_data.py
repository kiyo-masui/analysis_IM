"""Functions for getting data file names."""

import os
import re

data_dir = os.getenv('GBT10B_DATA')

def get_data_files(session, field) :
    """Gets a list of the data file names for each session and field.
    """
    
    session = int(session) # In case passed as a string.
    all_files = os.listdir(data_dir)
    out_files = []
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
