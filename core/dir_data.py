"""Functions for getting data file names."""

import os
import re
import glob

_data_dir = os.getenv('GBT10B_DATA')
_kiyo_data_dir = os.getenv('GBT10B_KM')+'data/'
# This line is cheating, will break for other people.
_guppi_data_dir = os.getenv('GBT10B_KM') + 'guppi_data/'

field_name_versions = {
    'wigglez15hr' : ['wigglez15hr', '15hr', 'wigglez15hrst'],
    'wigglez22hr' : ['wigglez22hr', '22hr', 'wigglez22hrst'],
    'wigglez1hr' : ['wigglez1hr', '1hr', 'wigglez1hrst', 'wigglez1hr_centre'],
    '3C286' : ['3C286', '3c286'],
    '3C348' : ['3C348', '3c348'],
    '3C48' : ['3C48', '3c48'],
    '3C67' : ['3C67', '3c67'],
}

data_dir = os.getenv('GBT_DATA')


def get_data_files(session_list, project='GBT10B_036', field='15hr',
                   type='ralongmap'):
    """Gets a list of data file names matching specified criteria.
    
    File names are retrieved for based on speified criteria.  The names include
    one level of directory name (this is the project name) but have the
    extension chopped off.

    Parameters
    ----------
    session_list : list of integers
        What sessions to use.
    project : string
        What project to use.
    field : string
        What source to use.  This can be specified in a variety of ways ('15hr,
        'wigglez15hr', 'wigglez1hr_centre', '3c67', '3C67' should all work.)
    type : string
        Scan type to use.

    Returns
    -------
    file_list : list of strings
        List of file names matching the criteria.  In includes the project name
        directory, and thus uniquly identifies the data files.

    Examples
    --------
    >>> get_data_files([89], project='GBT10B_036', field='1hr',
                       type='ralongmap')
    ["GBT10B_036/89_wigglez1hr_centre_ralongmap_82-91", ...]
    >>> get_data_files([72], project='GBT10B_036',field='3C48',
                       type='onoff')
    ["GBT10B_036/72_3C48_onoff_8-9", ...]

    Notes
    -----
    You can use unix wildcards in the `type` argument, but nowhere else.
    """
    
    for this_field, name_versions in field_name_versions.iteritems():
        if field in name_versions:
            field = this_field
            break
    else:
        raise ValueError("Invalid field name.")
    out_file_name_list = []
    for name_version in field_name_versions[field]:
        template = project + "/%02d_" + name_version + "_" + type + "_*.fits"
        template = data_dir + '/' + template
        for session in session_list:
            this_session_template = template%session
            files_this_session = glob.glob(this_session_template)
            for file_name in files_this_session:
                start_ind = file_name.find(project)
                formatted = file_name[start_ind:]
                formatted = formatted.split('.')[0]
                out_file_name_list.append(formatted)
    return out_file_name_list




#### All depricated. ####

def get_data_files_depricated(session_list, field, type=None) :
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
                           ".*_ralongmap_.*\.fits)")
            for file_name in all_files :
                # See if the file matches.
                if re.match(match_str, file_name) :
                    # Chop off extension and append to list to be returned.
                    root = file_name.split('.')[0]
                    out_files.append(root)
    
    # Remove files that don't have exactly 8 or 10 scans in them.
    good_scans = []
    for middle_string in out_files:
        rlong_index = middle_string.find('ralongmap')
        underscore_index = middle_string.find('_', rlong_index)
        dash_index = middle_string.find('-', underscore_index)
        first_scan = int(middle_string[underscore_index+1:dash_index])
        final_scan = int(middle_string[dash_index+1:])
        if (((final_scan - first_scan) == 7) 
            or ((final_scan - first_scan) == 9)) :
            good_scans.append(middle_string)

    return good_scans

def get_cal_files(session_list,calibrator,type=None) :
    """Gets a list of the cal file names for each session and calibrator
    calibrator=xxx (three numbers following 3c).
    """
    cal_files = []

    if type is None :
        all_files = os.listdir(_data_dir)
        for session in session_list:
            match_str = ('(%02d'%session+'_3c'+calibrator+'_onoff_'+".*\.raw.acs.fits)")
            for file_name in all_files:
                if re.match(match_str, file_name) :
                    root = file_name.split('.')[0]
                    cal_files.append(root)

    elif type == 'kiyo':
        all_files = os.listdir(_kiyo_data_dir)
        for session in session_list:
            match_str = ('(%02d'%session+'_3c'+calibrator+'_onoff_'+".*\.raw.acs.fits)")
            for file_name in all_files:
                if re.match(match_str, file_name) :
                    root = file_name.split('.')[0]
                    cal_files.append(root)

    elif type == 'guppi' :
        all_files = os.listdir(_guppi_data_dir) 
        for session in session_list: 
            match_str = ('(%02d'%session+'_3C'+calibrator+'_onoff_'+".*\.fits)")
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


