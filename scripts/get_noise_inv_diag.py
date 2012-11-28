"""Program to extract and save noise inverse diagonals.

Supply and input directory to search for noise_inv file and an output directory
to write noise_inv_diag files.

If only one directory is supplied, it is used for both input and output.
"""

import sys
import glob

import numpy as np

from core import algebra as al

def extract(in_dir, out_dir) :
    """Searches for noise_inv files, extracts the diagonal and writes it out.
    """
    
    files = glob.glob(in_dir + '/*noise_inv*.npy')
    for file_path in files:
        if 'noise_inv_diag' in file_path:
            continue
        file_name = file_path[len(in_dir):]
        parts = file_name.split('noise_inv')
        if len(parts) != 2:
            raise RuntimeError("'noise_inv' appears in file name more than"
                               " once.  Wasn't prepared for this.")
        out_path = out_dir + '/' + parts[0] + 'noise_inv_diag' + parts[1]
        mat = al.open_memmap(file_path, 'r')
        mat = al.make_mat(mat)
        mat_diag = mat.mat_diag()
        al.save(out_path, mat_diag)

if __name__ == "__main__" :
    if len(sys.argv) == 2:
        # Argument is the input directory
        directory = sys.argv[1]
        extract(directory, directory)
    elif len(sys.argv) == 3:
        extract(sys.argv[1], sys.argv[2])
    else :
        print ("Usage : python get_noise_inv_diag.py input_directory"
               " output_directory or"
               " python get_noise_inv_diag.py directory")
