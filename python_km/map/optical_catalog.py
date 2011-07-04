import numpy as np
import scipy as sp
import unittest
import profile
import time
import sys
from core import algebra 
from optparse import OptionParser
# TODO: make better parameter passing for catalog binning

def find_edges(axis, delta=False):
    """
    service function for bin_catalog_data which
    finds the bin edges for the histogram
    """
    if not delta:
        delta = axis[1] - axis[0]

    edges = np.array(axis) - delta/2.
    return np.append(edges, edges[-1]+delta)


def print_edges(sample, edges, name):
    print "Binning %s from range (%5.3g, %5.3g) into (%5.3g, %5.3g)" % (
           name, min(sample), max(sample), min(edges), max(edges))


def histogram3d(sample, xedges, yedges, zedges):
    """Make a 3D histogram from the sample and edge specification
    indices in the sample: 0=x, 1=y, 2=z; 
    histogramdd was having problems with the galaxy catalogs
    """
    numcatalog = sample.size
    x_size = xedges.size-1
    y_size = yedges.size-1
    z_size = zedges.size-1
    box_index = np.zeros(numcatalog)
    count_array = np.zeros((x_size+1)*(y_size+1)*(z_size+1))
    # the final array to return is the value within the bin 
    count_cube = np.zeros((x_size, y_size, z_size))

    # find which bin each galaxies lies in
    x_index = np.digitize(sample[:,0], xedges)
    y_index = np.digitize(sample[:,1], yedges)
    z_index = np.digitize(sample[:,2], zedges)

    # digitize puts values outside of the bins either to 0 or len(bins)
    x_out = np.logical_or((x_index == 0), (x_index == (x_size+1)))
    y_out = np.logical_or((y_index == 0), (y_index == (y_size+1)))
    z_out = np.logical_or((z_index == 0), (z_index == (z_size+1)))
    # now flag all those point which are inside the region
    box_in = np.logical_not(np.logical_or(np.logical_or(x_out, y_out), z_out))

    # the 0th bin center is recorded in the digitized index=1, so shift
    # also throw out points that are not in the volume
    x_index = x_index[box_in] - 1
    y_index = y_index[box_in] - 1
    z_index = z_index[box_in] - 1

    box_index = x_index + y_index*x_size + z_index*x_size*y_size

    # note that bincount will only count up to the largest object in the list,
    # which may be smaller than the dimension of the full count cube
    try:
        count_array[0:max(box_index)+1] = np.bincount(box_index)

        # make the x y and z axes which index the bincount output
        count_index = np.arange(x_size*y_size*z_size)
        z = count_index/(x_size*y_size)
        y = (count_index-x_size*y_size*z)/x_size
        x = count_index-x_size*y_size*z-x_size*y

        #count_cube[x,y,z] = count_array[x + y*x_size + z*x_size*y_size]
        count_cube[x,y,z] = count_array[count_index]
        #split_indices = cartesian((np.arange(z_size), np.arange(y_size), np.arange(x_size)))
        #count_cube[split_indices] = count_array[count_index]
    except MemoryError:
        print "histogram3d: all points out of the volume" 

    return count_cube


def bin_catalog_data(filename, freq_axis, ra_axis, dec_axis, skip_header=None, verbose=True):
    """
    bin catalog data onto a grid in RA, Dec, and frequency
    This currently assumes that all of the axes are uniformly spaced
    """
    # TODO: move this to a constants file
    nu_21cm_MHz = 1420.40575177

    # read the WiggleZ catalog and convert redshift axis to frequency
    ndtype=[('RA', float), ('Dec', float), ('z', float), ('r-mag', float), ('ijack', int)]
    # TODO: numpy seems to be an old version that does not have the skip_header
    # argument here! skiprows is identical
    output = np.genfromtxt(filename, dtype=ndtype, skiprows=skip_header)

    if verbose:
        print filename+": "+repr(output.dtype.names)+", n_records = "+repr(output.size)

    catalog_frequencies = nu_21cm_MHz*1.e6/(1+output['z'])
    num_catalog = output.size
    sample = np.zeros((num_catalog,3))
    sample[:,0] = catalog_frequencies
    sample[:,1] = output['RA']
    sample[:,2] = output['Dec']

    freq_edges = find_edges(freq_axis)
    ra_edges = find_edges(ra_axis)
    dec_edges = find_edges(dec_axis)

    if verbose:
        #print len(freq_axis), len(ra_axis), len(dec_axis)
        #print len(freq_edges), len(ra_edges), len(dec_edges)
        print_edges(sample[:,0], freq_edges, "frequency")
        print_edges(sample[:,1], ra_edges, "RA")
        print_edges(sample[:,2], dec_edges, "Dec")

    #print sample, freq_edges, ra_edges, dec_edges

    count_cube = histogram3d(sample, freq_edges, ra_edges, dec_edges)
    #count_cube, edges = np.histogramdd(sample, bins=[freq_edges, ra_edges, dec_edges])
    #print edges

    return count_cube


def bin_wigglez(filename=None):
    """process the WiggleZ optical catalog, binning against a template map
    given by the filename.
    In the 15h field, 
    ra range is 214. to 223.
    dec range is 0. to 4.
    freq range is 676. to 947.
    """

    np.set_printoptions(threshold=np.nan)
    if filename is None: 
        #filename = '/mnt/raid-project/gmrt/kiyo/wiggleZ/corr/sec_A_15hr_41-69_cleaned_clean_map_I.npy'
        filename ='/mnt/raid-project/gmrt/calinliv/wiggleZ/corr/test/sec_A_15hr_41-73_cleaned_clean_map_I_with_B.npy'

    root_data = "/cita/h/home-2/eswitzer/data/wigglez/"
    binned_data = "/cita/h/home-2/eswitzer/data/binned_wigglez/"
    n_random = 1000
    n_to_save = 100

    Template_map = algebra.make_vect(algebra.load(filename))
    freq_axis = Template_map.get_axis('freq')
    ra_axis = Template_map.get_axis('ra')
    dec_axis = Template_map.get_axis('dec')

    realmap_binning = bin_catalog_data(root_data+"reg15data.dat", freq_axis,
                                          ra_axis, dec_axis, skip_header=1)
    Map_wigglez = algebra.make_vect(realmap_binning, axis_names=('freq', 'ra', 'dec'))
    Map_wigglez.copy_axis_info(Template_map)
    algebra.save(binned_data+"reg15data.npy", Map_wigglez)

    selection_function = np.zeros(Template_map.shape)
    for i in range(n_random):
        wfile = root_data+"reg15rand"+str(i).zfill(3)+".dat"
        random_binning = bin_catalog_data(wfile, freq_axis, ra_axis,
                                             dec_axis, skip_header=1) 
        selection_function += random_binning
        if i < n_to_save:
            Map_wigglez_random = algebra.make_vect(random_binning, axis_names=('freq', 'ra', 'dec'))
            Map_wigglez_random.copy_axis_info(Template_map)
            algebra.save(binned_data+"reg15rand"+str(i).zfill(3)+".npy", Map_wigglez_random)

    # adding the real map back to the selection function is a kludge which
    # ensures the selection function is not zero where there is real data; this
    # should only be used in the limit of a small number of realizations of
    # random catalogs. (note: n_random+1)
    selection_function += realmap_binning
    selection_function /= float(n_random+1)

    Map_wigglez_selection = algebra.make_vect(selection_function, axis_names=('freq', 'ra', 'dec'))
    Map_wigglez_selection.copy_axis_info(Template_map)
    algebra.save(binned_data+"reg15selection.npy", Map_wigglez_selection)

    # now assume separability of the selection function    
    spatial_selection = np.sum(selection_function, axis=0) # 2D array 
    freq_selection = np.apply_over_axes(np.sum, selection_function, [1,2]) # 1D 
    separable_selection = (freq_selection * spatial_selection) # back to 3D
    separable_selection /= np.sum(freq_selection.flatten())

    Map_wigglez_separable = algebra.make_vect(separable_selection, axis_names=('freq', 'ra', 'dec'))
    Map_wigglez_separable.copy_axis_info(Template_map)
    algebra.save(binned_data+"reg15separable.npy", Map_wigglez_separable)


class CatalogGriddingTest(unittest.TestCase):
    """Unit test class for catalog gridding"""

    def test_simple(self):
        """bin a simple 3x3x3 array"""

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
        sample2 = np.array([[-1.,-1.,-1.]])
        result2 = histogram3d(sample2, edges, edges, edges)
        alternate2, histo_edges = np.histogramdd(sample2, 
                                                bins=[edges, edges, edges])

        answer2 = np.zeros((3,3,3), dtype=int)
        self.assertTrue(np.array_equal(answer2, result2))
        self.assertTrue(np.array_equal(alternate2, result2))

    def test_timing(self):
        """compare the timing of histogram3d and histogramdd"""
        # TODO: compare sum of two histogram methods; edge cases do not seem to match
        # TODO: speed up histogram3d class
        edges = np.array([0.,0.25,0.75,1.])
        sample = np.random.rand(1e7, 3)

        # profiling tools do not seem to work well with numpy
        start = time.clock()
        result = histogram3d(sample, edges, edges, edges)
        end = time.clock()
        print (end - start)/1000.
        alternate, histo_edges = np.histogramdd(sample,
                                                bins=[edges, edges, edges])
        endalt = time.clock()
        print (endalt - end)/1000.
        print result - alternate


if __name__ == '__main__':
    parser = OptionParser(usage="usage: %prog [options]",
                          version="%prog 1.0")
    parser.add_option("-t", "--test",
                      action="store_true",
                      dest="run_unittest",
                      default=False,
                      help="run some units tests on the catalog binning")
    parser.add_option("-f", "--file",
                      action="store", 
                      dest="template_file",
                      default=None,
                      help="File for radio cube to use as template",)
    (options, args) = parser.parse_args()

    #if len(args) != 1:
    #    parser.error("wrong number of arguments")
    print options
    print args

    if options.__dict__["run_unittest"]:
        del sys.argv[1:] # unittest would try to interpret these
        unittest.main()
    else:
        bin_wigglez()
