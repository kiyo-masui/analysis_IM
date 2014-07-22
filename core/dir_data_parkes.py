import os
import glob

def get_data_files(base_dir, source):
    """Gets a list of sdfits files in base_dir that contain the string source.
Removes the .sdfits ending.  Returned file list also removes base_dir from
start of string."""
    files = []
    dir_len = len(base_dir)
    for el in os.walk(base_dir):
        dir_files = glob.glob(el[0] + '/*' + source + '*.sdfits')
        for file in dir_files:
            files.append(file.split('.')[0][dir_len:])
    return files
