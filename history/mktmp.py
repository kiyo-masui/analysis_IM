import numpy as np
import core.algebra as algebra
"""Trivial code to generate template files.

arr = sp.zeros((256,64,32))

# the axis info for all stuff to be done on the 15hr field.
info = {'ra_delta': -0.075045715822411624,
        'dec_delta': 0.074999999999999997,
        'dec_centre': 2.0,
        'axes': ('freq', 'ra', 'dec'),
        'ra_centre': 217.90000000000001,
        'freq_centre': 799609375.0,
        'freq_delta': -781250.0,
        'type': 'vect'}

map_template = algebra.make_vect(arr, axis_names=('freq', 'ra', 'dec'))
map_template.info = info

0. get pep8 and pylint; python setup.py install --user; pep8 document
1. given left/right bounds on RA, top/bottom bounds on Dec return template
for which the survey map is a submap (full GBT, wiggleZ)
2. simulate signal realizations on the volume (ERS to give WiggleZ survey dim)
3. hammer out ~eswitzer/code/analysis_IM/map/hitmap.py -- make very fine
3.' boolean array of where the subregions overlap

A. Lewis
"""

print "1"


def mktmp(rgn_i,rgn_j,rgn_k,srgn_i1,srgn_i2,srgn_j1,srgn_j2,srgn_k1,srgn_k2,outfilename):
    """Write to disk a file representing an empty matrix of given dimensions. Also write an identically
    shaped array of booleans, which are true if the index points to the subregion.
    rgn_i/j/k  : the dimensions of the full region to be simulated
        srgn_i/j/k : the dimensions of the deep integration subregion
    outfilename: the name of the file to be created
    """


    regiontype = np.zeros((rgn_i,rgn_j,rgn_k), bool)

    array = np.zeros((rgn_i,rgn_j,rgn_k))

    for i in range(0,rgn_i):
        for j in range(0,rgn_j):
            for k in range(0,rgn_k):
                if (i>=(srgn_i1-1) and i<=(srgn_i2-1)):
                    if (j>=(srgn_j1-1) and j<=(srgn_j2-1)):
                        if (k>=(srgn_k1-1) and k<=(srgn_k2-1)):
                            regiontype[i,j,k]=True
            else:
                regiontype[i,j,k]=False

    region=algebra.info_array(array)
    regiontypename = 'bool' + outfilename
    np.save(regiontypename, regiontype)
    algebra.save(outfilename,region)
    print "done"
    template_map = algebra.make_vect(algebra.load(outfilename))
