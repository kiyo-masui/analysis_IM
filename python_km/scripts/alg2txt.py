"""Program to convert algebra objects to text files.

This is mainly written to convert maps and noise matricies to text files for
other peoples pipelines.  This program on treats a limited number of cases like
3D vects, and the diagonals of mat objects.
"""

import sys

import numpy as np

from core import algebra

def tofile(fname, data) :
    """Writes a 3D array to a text file.

    Writes a 3D array to a text file using comment lines to separate slices of
    the array.  The array dimensions are written in a comment in the first line 
    of the array.

    Parameters
    ----------
    fname : string
        File name to which to write the array.
    array : array like
        Data to write to file.

    See Also
    --------
    np.savetxt : Only works for < 2D arrays.
    """

    if array.ndim != 3 :
        raise ValueError("Array must be 3D.")
    # Write the array to disk
    with file(fname, 'w') as outfile:
        # I'm writing a header here just for the sake of readability
        # Any line starting with "#" will be ignored by numpy.loadtxt
        outfile.write('# Array shape: {0}\n'.format(data.shape))

        # Iterating through a ndimensional array produces slices along
        # the last axis. This is equivalent to data[i,:,:] in this case
        for data_slice in data:

            # The formatting string indicates that I'm writing out
            # the values in left-justified columns 7 characters in width
            # with 2 decimal places.
            np.savetxt(outfile, data_slice)

            # Writing out a break to indicate different slices...
            outfile.write('# New slice\n')

if __name__ == "__main__" :
    if len(sys.argv) == 2 :
        # Argument should just be a .npy file.
        array = algebra.load(sys.argv[1])
        out_fname = sys.argv[1].split('/')[-1][:-4] + '.txt'
        tofile(out_fname, array)
    elif len(sys.argv) == 3 and sys.argv[1] == str("diag") :
        # Second argument should be a .npy file that should be interpreted as a
        # matrix and we want to save the diagonal.
        mat = algebra.open_memmap(sys.argv[2])
        mat = algebra.make_mat(mat)
        array = mat.mat_diag()
        out_fname = sys.argv[2].split('/')[-1][:-4] + '.txt'
        tofile(out_fname, array)
    else :
        print ("Usage : python alg2txt.py [input file] or"
               " python alg2txt.py diag [input file]")
