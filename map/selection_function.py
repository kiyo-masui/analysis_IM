"""find the selection function using a monte carlo over catalogs"""
import numpy as np
import shelve
import random
import multiprocessing as mp
import copy
from core import algebra
from core import constants as cc
from utils import binning
from utils import data_paths

def randomize_catalog_redshifts(catalog, nzilookupfile):
    """
    re-draw the redshifts of a catalog from N(z) according to observation
    priority -- from Chris Blake
    """
    ndtype = [('CDF', float), ("z0", float), ("z1", float), ("z2", float),
              ("z3", float),  ("z4", float), ("z5", float)]
    nzilookup = np.genfromtxt(nzilookupfile, dtype=ndtype, skiprows=1)

    # find the priority zones of each random source
    whpriority = [np.where(catalog["sec"] == ipri)[0] for ipri in range(1, 7)]
    numpriority = [(index, entry.size) for (index, entry)
                    in zip(range(1, 7), whpriority)]
    randsample = [(index, np.random.random_sample(ssize)) for \
                  (index, ssize) in numpriority]
    randz = [np.interp(draw, nzilookup["CDF"],
                       nzilookup["z" + repr(index - 1)])
                    for (index, draw) in randsample]
    for zone in range(0, 6):
        catalog["z"][whpriority[zone]] = randz[zone]

    #print catalog["z"]
    #sys.exit()

    return catalog


def wrap_bin_catalog_data(entry):
    """hand some tuple of information to the binning code
    """
    (ranindex, rancat, freq_axis, ra_axis, dec_axis) = entry
    #print ranindex
    return bin_catalog_data(rancat, freq_axis, ra_axis, dec_axis)


def template_map_axes(filename):
    """Open a numpy array map and extract its axis/etc. information
    """
    print "using the volume template file: " + filename
    template_map = algebra.make_vect(algebra.load(filename))
    freq_axis = template_map.get_axis('freq')
    ra_axis = template_map.get_axis('ra')
    dec_axis = template_map.get_axis('dec')
    return (freq_axis, ra_axis, dec_axis, template_map.shape, template_map)


# TODO: remove A/B split in averaging? probably not needed
def estimate_selection_function(fieldname, template_file,
                                catalog_shelvefile="./randcatdata.shelve",
                                generate_rancat_shelve=True,
                                complete=False):
    """Estimate the selection function using random draws in C. Blake N(z)
    PDFs.

    Parameters:
    -----------
    fieldname: string
        name of field, e.g. "15hr", "22hr" or "1hr"
    template_file: string
        file (not db key) to use to grab dimensions for the binnins
    catalog_shevefile: string keyword
        the name of a intermediate (binary) catalog format
        this speeds the computation.
    generate_rancat_shelve: boolean
        if the catalog shelve is already generated, set
        this to False and it will use the pre-existing database
    complete: boolean keyword
        use the complete WiggleZ regions
    """
    datapath_db = data_paths.DataPath()
    if complete:
        ctag = "complete_"
    else:
        ctag = ""

    db_key = "WiggleZ_%s_%smontecarlo" % (fieldname, ctag)
    outfile_selection = datapath_db.fetch(db_key, intend_write=True,
                           purpose="WiggleZ MC selection function")

    db_key = "WiggleZ_%s_mock_catalog" % fieldname
    infile_mock = datapath_db.fetch(db_key, intend_read=True,
                             purpose="WiggleZ real data catalog", silent=True)

    db_key = "WiggleZ_%s_priority_table" % fieldname
    priority_file = datapath_db.fetch(db_key, intend_read=True,
                        purpose="WiggleZ mock cat priority file", silent=True)

    n_rand_cats = len(infile_mock[0])
    chunking_size = 10  # break the averaging into pooled multiprocess jobs
    num_chunks = 100000 # 9000 for testing, 100000 for production

    # read the WiggleZ catalog and convert redshift axis to frequency
    randdata = shelve.open(catalog_shelvefile)
    if generate_rancat_shelve:
        ndtype = [('RA', float), ('Dec', float), ('z', float),
                  ('r-mag', float), ('ijack', int), ('sec', int)]
        for index in infile_mock[0]:
            filename = infile_mock[1][index]
            print "loading: " + filename
            randdata[index] = np.genfromtxt(filename, dtype=ndtype,
                                               skiprows=1)
        print "done generating random catalog shelve lookup"

    (freq_axis, ra_axis, dec_axis, template_shape, template_map) = \
        template_map_axes(filename=template_file)

    selection_function = np.zeros(template_shape)
    sel_func_a = np.zeros(template_shape)
    sel_func_b = np.zeros(template_shape)
    for iternum in range(0, num_chunks):
        runlist_a = []
        runlist_b = []
        for count in range(0, chunking_size):
            ran_a = random.randint(0, n_rand_cats - 1)
            ran_b = random.randint(0, n_rand_cats - 1)
            # TODO: do we need deep copy here, or paranoid?
            rancat_a = randomize_catalog_redshifts(
                            copy.deepcopy(randdata[repr(ran_a)]),
                            priority_file)
            rancat_b = randomize_catalog_redshifts(
                            copy.deepcopy(randdata[repr(ran_b)]),
                            priority_file)
            runlist_a.append((ran_a, rancat_a, freq_axis, ra_axis, dec_axis))
            runlist_b.append((ran_b, rancat_b, freq_axis, ra_axis, dec_axis))

        chunk_sel_func_a = np.zeros(template_shape)
        chunk_sel_func_b = np.zeros(template_shape)

        # leave 3 processors free to be nice
        pool = mp.Pool(processes=(mp.cpu_count() - 3))
        results_a = pool.map(wrap_bin_catalog_data, runlist_a)
        results_b = pool.map(wrap_bin_catalog_data, runlist_b)
        pool.close()

        for resultitem in results_a:
            chunk_sel_func_a += resultitem

        for resultitem in results_b:
            chunk_sel_func_b += resultitem

        chunk_sel_func_a /= float(len(results_a))
        chunk_sel_func_b /= float(len(results_b))

        sel_func_a += chunk_sel_func_a
        sel_func_b += chunk_sel_func_b

        print np.std((sel_func_a - sel_func_b) / float(iternum + 1)), \
              np.mean(sel_func_a) / float(iternum + 1), \
              np.mean(sel_func_b) / float(iternum + 1), iternum

    selection_function = (sel_func_a +
                          sel_func_b) / 2. / float(num_chunks)
    randdata.close()

    map_wigglez_selection = algebra.make_vect(selection_function,
                                              axis_names=('freq', 'ra', 'dec'))
    map_wigglez_selection.copy_axis_info(template_map)
    print np.min(map_wigglez_selection), np.max(map_wigglez_selection)
    algebra.save(outfile_selection, map_wigglez_selection)


def montecarlo_selection():
    map_root = '/mnt/raid-project/gmrt/tcv/maps/'
    template_15hr = map_root + 'sec_A_15hr_41-90_clean_map_I.npy'
    template_22hr = map_root + 'sec_A_22hr_41-90_clean_map_I.npy'
    template_1hr = map_root + 'sec_A_1hr_41-16_clean_map_I.npy'

    estimate_selection_function('15hr', template_15hr,
                                catalog_shelvefile="./rand15hrcatdata.shelve",
                                generate_rancat_shelve=True)

    #estimate_selection_function('22hr', template_22hr,
    #                            catalog_shelvefile="./rand22hrcatdata.shelve",
    #                            generate_rancat_shelve=True)


if __name__ == '__main__':
    montecarlo_selection()
