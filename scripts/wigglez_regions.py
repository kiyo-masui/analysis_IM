from utils import data_paths
import numpy as np
from core import constants as cc

def minmax(vector, previous):
    vmax = vector.max()
    vmin = vector.min()
    #print previous, vmin, vmax

    if previous[1] != None:
        if (vmax > previous[1]):
            previous[1] = vmax
    else:
        previous[1] = vmax

    if previous[0] != None:
        if (vmin < previous[0]):
            previous[0] = vmin
    else:
        previous[0] = vmin

    return previous

def find_wigglez_region(fieldname):
    """Find the RA/Dec edges of a WiggleZ field by exhaustion
    """
    datapath_db = data_paths.DataPath()

    db_key = "WiggleZ_%s_mock_catalog" % fieldname
    infile_mock = datapath_db.fetch(db_key, intend_read=True,
                             purpose="WiggleZ real data catalog", silent=True)

    n_rand_cats = len(infile_mock[0])

    ra_minmax = [None, None]
    dec_minmax = [None, None]
    freq_minmax = [None, None]

    ndtype = [('RA', float), ('Dec', float), ('z', float),
              ('r-mag', float), ('ijack', int), ('sec', int)]
    for index in infile_mock[0]:
        filename = infile_mock[1][index]
        #print "loading: " + filename
        catalog = np.genfromtxt(filename, dtype=ndtype,
                                skiprows=1)
        freq_vec = cc.freq_21cm_MHz * 1.e6 / (1 + catalog['z'])
        ra_vec = catalog['RA']
        dec_vec = catalog['Dec']

        ra_minmax = minmax(ra_vec, ra_minmax)
        dec_minmax = minmax(dec_vec, dec_minmax)
        freq_minmax = minmax(freq_vec, freq_minmax)
        #print ra_vec.min(), ra_vec.max(), ra_minmax
        #print dec_vec.min(), dec_vec.max(), dec_minmax
        #print freq_vec.min(), freq_vec.max(), freq_minmax

    return ra_minmax, dec_minmax, freq_minmax

if __name__ == "__main__":
    field_15hr_bounds = find_wigglez_region('15hr')
    print field_15hr_bounds
    field_22hr_bounds = find_wigglez_region('22hr')
    print field_22hr_bounds
    field_1hr_bounds = find_wigglez_region('1hr')
    print field_1hr_bounds
