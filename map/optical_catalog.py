"""
A set of functions to bin optical catalog data into cubes

    In the 15h field,
    ra range is 214. to 223.
    dec range is 0. to 4.
    freq range is 676. to 947.
"""
import numpy as np
import unittest
import time
import shelve
import random
import multiprocessing
import copy
from core import algebra
from core import constants as cc
# TODO: make better parameter passing for catalog binning
# TODO: put ranges of other fields in docstring


def find_edges(axis, delta=False):
    """
    service function for bin_catalog_data which
    finds the bin edges for the histogram
    """
    if not delta:
        delta = axis[1] - axis[0]

    edges = np.array(axis) - delta / 2.
    return np.append(edges, edges[-1] + delta)


def print_edges(sample, edges, name):
    """print bin edges for a catalog
    """
    print "Binning %s from range (%5.3g, %5.3g) into (%5.3g, %5.3g)" % (
           name, min(sample), max(sample), min(edges), max(edges))


def histogram3d(sample, xedges, yedges, zedges):
    """Make a 3D histogram from the sample and edge specification
    indices in the sample: 0=x, 1=y, 2=z;
    histogramdd was having problems with the galaxy catalogs
    """
    numcatalog = sample.size
    x_size = xedges.size - 1
    y_size = yedges.size - 1
    z_size = zedges.size - 1
    box_index = np.zeros(numcatalog)
    count_array = np.zeros((x_size + 1) * (y_size + 1) * (z_size + 1))
    # the final array to return is the value within the bin
    count_cube = np.zeros((x_size, y_size, z_size))

    # find which bin each galaxies lies in
    x_index = np.digitize(sample[:, 0], xedges)
    y_index = np.digitize(sample[:, 1], yedges)
    z_index = np.digitize(sample[:, 2], zedges)

    # digitize puts values outside of the bins either to 0 or len(bins)
    x_out = np.logical_or((x_index == 0), (x_index == (x_size + 1)))
    y_out = np.logical_or((y_index == 0), (y_index == (y_size + 1)))
    z_out = np.logical_or((z_index == 0), (z_index == (z_size + 1)))
    # now flag all those point which are inside the region
    box_in = np.logical_not(np.logical_or(np.logical_or(x_out, y_out), z_out))

    # the 0th bin center is recorded in the digitized index=1, so shift
    # also throw out points that are not in the volume
    x_index = x_index[box_in] - 1
    y_index = y_index[box_in] - 1
    z_index = z_index[box_in] - 1

    box_index = x_index + y_index * x_size + z_index * x_size * y_size

    # note that bincount will only count up to the largest object in the list,
    # which may be smaller than the dimension of the full count cube
    try:
        count_array[0:max(box_index) + 1] = np.bincount(box_index)

        # make the x y and z axes which index the bincount output
        count_index = np.arange(x_size * y_size * z_size)
        zind = count_index / (x_size * y_size)
        yind = (count_index - x_size * y_size * zind) / x_size
        xind = count_index - x_size * y_size * zind - x_size * yind

        #count_cube[xind, yind, zind] = count_array[xind + yind * x_size +
        #                                           zind * x_size * y_size]
        count_cube[xind, yind, zind] = count_array[count_index]
        #split_indices = cartesian((np.arange(z_size),
        #                           np.arange(y_size),
        #                           np.arange(x_size)))
        #count_cube[split_indices] = count_array[count_index]
    except MemoryError:
        print "histogram3d: all points out of the volume"

    return count_cube


def bin_catalog_data(catalog, freq_axis, ra_axis,
                     dec_axis, verbose=False):
    """
    bin catalog data onto a grid in RA, Dec, and frequency
    This currently assumes that all of the axes are uniformly spaced
    """
    catalog_frequencies = cc.freq_21cm_MHz * 1.e6 / (1 + catalog['z'])
    num_catalog = catalog.size
    sample = np.zeros((num_catalog, 3))
    sample[:, 0] = catalog_frequencies
    sample[:, 1] = catalog['RA']
    sample[:, 2] = catalog['Dec']

    freq_edges = find_edges(freq_axis)
    ra_edges = find_edges(ra_axis)
    dec_edges = find_edges(dec_axis)

    if verbose:
        #print len(freq_axis), len(ra_axis), len(dec_axis)
        #print len(freq_edges), len(ra_edges), len(dec_edges)
        print_edges(sample[:, 0], freq_edges, "frequency")
        print_edges(sample[:, 1], ra_edges, "RA")
        print_edges(sample[:, 2], dec_edges, "Dec")

    #print sample, freq_edges, ra_edges, dec_edges

    count_cube = histogram3d(sample, freq_edges, ra_edges, dec_edges)
    #count_cube, edges = np.histogramdd(sample,
    #                                    bins=[freq_edges,
    #                                          ra_edges, dec_edges])
    #print edges
    return count_cube


def bin_catalog_file(filename, freq_axis, ra_axis,
                     dec_axis, skip_header=None, verbose=True):
    """Bin the catalog given in `filename` using the given freq, ra, and dec
    axes.
    """

    # read the WiggleZ catalog and convert redshift axis to frequency
    ndtype = [('RA', float), ('Dec', float), ('z', float),
              ('r-mag', float), ('ijack', int)]
    # TODO: numpy seems to be an old version that does not have the skip_header
    # argument here! skiprows is identical
    output = np.genfromtxt(filename, dtype=ndtype, skiprows=skip_header)

    if verbose:
        print filename + ": " + repr(output.dtype.names) + \
              ", n_records = " + repr(output.size)

    return bin_catalog_data(output, freq_axis, ra_axis, dec_axis,
                            verbose=verbose)


def wrap_bin_catalog_data(entry):
    """hand some tuple of information to the binning code
    """
    (ranindex, rancat, freq_axis, ra_axis, dec_axis) = entry
    #print ranindex
    return bin_catalog_data(rancat, freq_axis, ra_axis, dec_axis)


def template_map_axes(filename):
    """Open a numpy array map and extract its axis/etc. information
    """
    #if filename is None:
    #    root_template = '/mnt/raid-project/gmrt/calinliv/wiggleZ/corr/test/'
    #    filename = root_template + \
    #               'sec_A_15hr_41-73_cleaned_clean_map_I_with_B.npy'

    print "using the volume template file: " + filename
    template_map = algebra.make_vect(algebra.load(filename))
    freq_axis = template_map.get_axis('freq')
    ra_axis = template_map.get_axis('ra')
    dec_axis = template_map.get_axis('dec')
    return (freq_axis, ra_axis, dec_axis, template_map.shape, template_map)


# TODO: need printoptions?
def bin_wigglez(fieldname, template_file):
    """process the WiggleZ optical catalog, binning against a template map
    given by the template_file.
    """

    #np.set_printoptions(threshold=np.nan)
    root_data = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/"
    binned_data = root_data + "binned_wiggleZ/"
    root_catalogs = root_data + "wiggleZ_catalogs/"
    catalog_subdir = fieldname + "hr/"
    n_rand_cats = 1000
    n_to_save = 100

    outfile_data = binned_data + "reg" + fieldname + "data.npy"
    outfile_selection = binned_data + "reg" + fieldname + "selection.npy"
    outfile_separable = binned_data + "reg" + fieldname + "separable.npy"
    outfile_rand = binned_data + "reg" + fieldname + "rand"

    infile_data = root_catalogs + catalog_subdir + "reg" + fieldname + "data.dat"
    rootrand = root_catalogs + catalog_subdir + "reg" + fieldname + "rand"

    randlist = [(index, (rootrand + "%03d.dat" % index))
                for index in range(0, n_rand_cats)]

    (freq_axis, ra_axis, dec_axis, template_shape, template_map) = \
        template_map_axes(filename=template_file)

    realmap_binning = bin_catalog_file(infile_data, freq_axis,
                                          ra_axis, dec_axis, skip_header=1)
    map_wigglez = algebra.make_vect(realmap_binning,
                                    axis_names=('freq', 'ra', 'dec'))
    map_wigglez.copy_axis_info(template_map)
    algebra.save(outfile_data, map_wigglez)

    selection_function = np.zeros(template_map.shape)
    for (randindex, randfile) in randlist:
        random_binning = bin_catalog_file(randfile, freq_axis, ra_axis,
                                             dec_axis, skip_header=1)
        selection_function += random_binning
        if randindex < n_to_save:
            map_wigglez_random = algebra.make_vect(random_binning,
                                         axis_names=('freq', 'ra', 'dec'))
            map_wigglez_random.copy_axis_info(template_map)
            algebra.save(outfile_rand + str(randindex).zfill(3) + ".npy",
                         map_wigglez_random)

    # adding the real map back to the selection function is a kludge which
    # ensures the selection function is not zero where there is real data; this
    # should only be used in the limit of a small number of realizations of
    # random catalogs. (note: n_random + 1)
    selection_function += realmap_binning
    selection_function /= float(n_rand_cats + 1)

    map_wigglez_selection = algebra.make_vect(selection_function,
                                              axis_names=('freq', 'ra', 'dec'))
    map_wigglez_selection.copy_axis_info(template_map)
    algebra.save(outfile_selection, map_wigglez_selection)

    # now assume separability of the selection function
    spatial_selection = np.sum(selection_function, axis=0)  # 2D array
    freq_selection = np.apply_over_axes(np.sum,
                                        selection_function, [1, 2])  # 1D
    separable_selection = (freq_selection * spatial_selection)  # back to 3D
    separable_selection /= np.sum(freq_selection.flatten())

    map_wigglez_separable = algebra.make_vect(separable_selection,
                                              axis_names=('freq', 'ra', 'dec'))
    map_wigglez_separable.copy_axis_info(template_map)
    algebra.save(outfile_separable, map_wigglez_separable)


def randomize_catalog_redshifts(catalog, fieldname):
    """
    re-draw the redshifts of a catalog from N(z) according to observation
    priority -- from Chris Blake
    """
    nzilookupfile = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/wiggleZ_catalogs/"
    nzilookupfile += fieldname + "hr/"
    nzilookupfile += "nzpri_reg" + fieldname + "_tzuchingcats.dat"
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


# TODO: remove A/B split in averaging? probably not needed
def estimate_selection_function(fieldname, template_file,
                                catalog_shelvefile="randcatdata.shelve",
                                generate_rancat_shelve=True):
    """Estimate the selection function using random draws in C. Blake N(z)
    PDFs.
    fieldname is the identifier for the field; "01", "15", or "22" hr
    template_file -- numpy template to use for the axes of the output
    catalog_shevefile is a keyword giving the name of a intermediate (binary)
        catalog format -- this speeds the computation.
    generate_rancat_shelve -- if the catalog shelve is already generated, set
        this to False and it will use the pre-existing database
    """
    root_data = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/"
    binned_data = root_data + "binned_wiggleZ/"
    root_catalogs = root_data + "wiggleZ_catalogs/"
    catalog_subdir = fieldname + "hr/"

    outfile_selection = binned_data + "reg" + fieldname + "montecarlo.npy"
    n_rand_cats = 1000
    chunking_size = 10  # break the averaging into pooled multiprocess jobs
    #num_chunks = 9000
    #num_chunks = 20000
    num_chunks = 60000

    rootrand = root_catalogs + catalog_subdir + "reg" + fieldname + "rand"
    randlist = [(repr(index), (rootrand + "%03d.dat" % index))
                for index in range(0, n_rand_cats)]

    # read the WiggleZ catalog and convert redshift axis to frequency
    randdata = shelve.open(binned_data + catalog_shelvefile)
    if generate_rancat_shelve:
        ndtype = [('RA', float), ('Dec', float), ('z', float),
                  ('r-mag', float), ('ijack', int), ('sec', int)]
        #randdata = {}
        for entry in randlist:
            print "loading: " + entry[1]
            randdata[entry[0]] = np.genfromtxt(entry[1], dtype=ndtype,
                                               skiprows=1)
        randdata.close()

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
                            fieldname)
            rancat_b = randomize_catalog_redshifts(
                            copy.deepcopy(randdata[repr(ran_b)]),
                            fieldname)
            runlist_a.append((ran_a, rancat_a, freq_axis, ra_axis, dec_axis))
            runlist_b.append((ran_b, rancat_b, freq_axis, ra_axis, dec_axis))

        chunk_sel_func_a = np.zeros(template_shape)
        chunk_sel_func_b = np.zeros(template_shape)

        # leave 3 processors free to be nice
        pool = multiprocessing.Pool(processes=(multiprocessing.cpu_count()-3))
        resultsA = pool.map(wrap_bin_catalog_data, runlist_a)
        resultsB = pool.map(wrap_bin_catalog_data, runlist_b)

        for resultitem in resultsA:
            chunk_sel_func_a += resultitem

        for resultitem in resultsB:
            chunk_sel_func_b += resultitem

        chunk_sel_func_a /= float(len(resultsA))
        chunk_sel_func_b /= float(len(resultsB))

        sel_func_a += chunk_sel_func_a
        sel_func_b += chunk_sel_func_b

        print np.std((sel_func_a - sel_func_b) / float(iternum + 1)), \
              np.mean(sel_func_a) / float(iternum + 1), \
              np.mean(sel_func_b) / float(iternum + 1)

    selection_function = (sel_func_a +
                          sel_func_b) / 2. / float(num_chunks)
    print np.mean(selection_function)

    map_wigglez_selection = algebra.make_vect(selection_function,
                                              axis_names=('freq', 'ra', 'dec'))
    map_wigglez_selection.copy_axis_info(template_map)
    print np.min(map_wigglez_selection), np.max(map_wigglez_selection)
    algebra.save(outfile_selection, map_wigglez_selection)


class CatalogGriddingTest(unittest.TestCase):
    """Unit test class for catalog gridding
    """

    def test_simple(self):
        """bin a simple 3x3x3 array
        """

        parent_axis = np.array([0.25, 0.75, 1.25])
        edges = find_edges(parent_axis)
        self.assertTrue(np.array_equal(edges, [0., 0.5, 1., 1.5]))

        # test a sample (with some outliers)
        sample = np.array([[0., 0., 0.],
                           [0.75, 0., 0.],
                           [1.25, 0., 0.],
                           [1.75, 0., 0.],
                           [0., 0., 0.],
                           [0.75, 0.75, 0.75],
                           [1.25, 1.25, 1.25]])

        result = histogram3d(sample, edges, edges, edges)
        alternate, histo_edges = np.histogramdd(sample,
                                                bins=[edges, edges, edges])

        answer = np.array([[[2,  0,  0],
                            [0,  0,  0],
                            [0,  0,  0]],
                           [[1,  0,  0],
                            [0,  1,  0],
                            [0,  0,  0]],
                           [[1,  0,  0],
                            [0,  0,  0],
                            [0,  0,  1]]])

        self.assertTrue(np.array_equal(answer, result))
        self.assertTrue(np.array_equal(alternate, result))

        # test the case where no points are in the volume
        sample2 = np.array([[-1., -1., -1.]])
        result2 = histogram3d(sample2, edges, edges, edges)
        alternate2, histo_edges = np.histogramdd(sample2,
                                                bins=[edges, edges, edges])

        answer2 = np.zeros((3, 3, 3), dtype=int)
        self.assertTrue(np.array_equal(answer2, result2))
        self.assertTrue(np.array_equal(alternate2, result2))

    def test_timing(self):
        """compare the timing of histogram3d and histogramdd"""
        # TODO: compare sum of two histogram methods;
        # edge cases do not seem to match
        # TODO: speed up histogram3d class
        edges = np.array([0., 0.25, 0.75, 1.])
        sample = np.random.rand(1e7, 3)

        # profiling tools do not seem to work well with numpy
        start = time.clock()
        result = histogram3d(sample, edges, edges, edges)
        end = time.clock()
        print (end - start) / 1000.
        alternate, histo_edges = np.histogramdd(sample,
                                                bins=[edges, edges, edges])
        endalt = time.clock()
        print (endalt - end) / 1000.
        print result - alternate


def generate_sel_15hr_22hr_gbtregion():
    """generate selection functions for the 15hr and 22hr field based on GBT
    map template for those regions
    """
    calinliv_dir = '/mnt/raid-project/gmrt/calinliv/wiggleZ/'

    root_template = calinliv_dir + "corr/test/"
    template_mapname = root_template + \
                    'sec_A_15hr_41-73_cleaned_clean_map_I_with_B.npy'

    bin_wigglez("15", template_mapname)

    root_template = calinliv_dir + "corr/84_ABCD_all_15_modes/"
    template_mapname = root_template + \
                    'sec_A_22hr_41-84_cleaned_clean_map_I_with_C.npy'

    bin_wigglez("22", template_mapname)


def generate_MCsel_15hr_gbtregion():
    """generate selection functions for the 22hr field using a Monte Carlo
    over random source positions, consistently with priority zones from
    C. Blake
    """

    calinliv_dir = '/mnt/raid-project/gmrt/calinliv/wiggleZ/'
    root_template = calinliv_dir + "corr/test/"
    template_mapname = root_template + \
                    'sec_A_15hr_41-73_cleaned_clean_map_I_with_B.npy'

    estimate_selection_function('15', template_mapname,
                                catalog_shelvefile="randcatdata15.shelve",
                                generate_rancat_shelve=True)


def generate_MCsel_22hr_gbtregion():
    """generate selection functions for the 22hr field using a Monte Carlo
    over random source positions, consistently with priority zones from
    C. Blake
    """

    calinliv_dir = '/mnt/raid-project/gmrt/calinliv/wiggleZ/'
    root_template = calinliv_dir + "corr/84_ABCD_all_15_modes/"
    template_mapname = root_template + \
                    'sec_A_22hr_41-84_cleaned_clean_map_I_with_C.npy'

    estimate_selection_function('22', template_mapname,
                                catalog_shelvefile="randcatdata.shelve",
                                generate_rancat_shelve=False)


if __name__ == '__main__':
    #bin_wigglez()
    generate_MCsel_15hr_gbtregion()
    #generate_MCsel_22hr_gbtregion()
    #generate_GBT_sel_15hr_22hr_gbtregion()
