"""
A set of functions to bin optical catalog data into cubes
http://astronomy.swin.edu.au/~cblake/tzuching.html
    In the 15h field,
    ra range is 214. to 223.
    dec range is 0. to 4.
    freq range is 676. to 947.

    map_root = '/mnt/raid-project/gmrt/tcv/maps/'
    template_15hr = map_root + 'sec_A_15hr_41-90_clean_map_I.npy'
    template_22hr = map_root + 'sec_A_22hr_41-90_clean_map_I.npy'
    template_1hr = map_root + 'sec_A_1hr_41-16_clean_map_I.npy'
    template_15hr = 'templates/wigglez_15hr_complete.npy'
    template_22hr = 'templates/wigglez_22hr_complete.npy'
    template_1hr = 'templates/wigglez_1hr_complete.npy'
"""
import numpy as np
import shelve
import random
import copy
from core import algebra
from core import constants as cc
from utils import binning
from utils import data_paths
from kiyopy import parse_ini
# TODO: make better parameter passing for catalog binning


def bin_catalog_data(catalog, freq_axis, ra_axis,
                     dec_axis, debug=False, use_histogramdd=False):
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

    freq_edges = binning.find_edges(freq_axis)
    ra_edges = binning.find_edges(ra_axis)
    dec_edges = binning.find_edges(dec_axis)

    if debug:
        binning.print_edges(sample[:, 0], freq_edges, "frequency")
        binning.print_edges(sample[:, 1], ra_edges, "RA")
        binning.print_edges(sample[:, 2], dec_edges, "Dec")
        print sample, freq_edges, ra_edges, dec_edges

    if use_histogramdd:
        count_cube, edges = np.histogramdd(sample,
                                            bins=[freq_edges,
                                                  ra_edges, dec_edges])
        print edges
    else:
        count_cube = binning.histogram3d(sample, freq_edges,
                                         ra_edges, dec_edges)

    return count_cube


def bin_catalog_file(filename, freq_axis, ra_axis,
                     dec_axis, skip_header=None, debug=False, mock=False):
    """Bin the catalog given in `filename` using the given freq, ra, and dec
    axes.
    """

    # read the WiggleZ catalog and convert redshift axis to frequency
    if mock:
        ndtype = [('RA', float), ('Dec', float), ('z', float),
                  ('r-mag', float), ('ijack', int), ('ipri', int)]
    else:
        ndtype = [('RA', float), ('Dec', float), ('z', float),
                  ('r-mag', float), ('ijack', int)]

    # TODO: numpy seems to be an old version that does not have the skip_header
    # argument here! skiprows is identical
    catalog = np.genfromtxt(filename, dtype=ndtype, skiprows=skip_header)

    if debug:
        print filename + ": " + repr(catalog.dtype.names) + \
              ", n_records = " + repr(catalog.size)

    return bin_catalog_data(catalog, freq_axis, ra_axis, dec_axis,
                            debug=debug)


binwigglezparams_init = {
        "infile_data": "WiggleZ_B_catalog_data",
        "infile_mock": "WiggleZ_B_mock_catalog",
        "outfile_data": "WiggleZ_B_Bbinned_data",
        "outfile_mock": "WiggleZ_B_Bmock",
        "outfile_deltadata": "WiggleZ_B_Bdelta_binned_data",
        "outfile_deltamock": "WiggleZ_B_Bdelta_mock",
        "outfile_selection": "WiggleZ_B_Bselection",
        "outfile_separable": "WiggleZ_B_Bseparable_selection",
        "template_file": "file"
               }
binwigglezprefix = 'binwigglez_'


class BinWigglez(object):
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        self.datapath_db = data_paths.DataPath()

        if parameter_file:
            self.params = parse_ini.parse(parameter_file, binwigglezparams_init,
                                          prefix=binwigglezprefix)

        # gather names of the input catalogs
        self.infile_data = self.datapath_db.fetch(self.params['infile_data'],
                                             intend_read=True,
                                             purpose="WiggleZ data catalog")

        self.infile_mock = self.datapath_db.fetch(self.params['infile_mock'],
                                             intend_read=True,
                                             purpose="WiggleZ mock catalog",
                                             silent=True)

        # gather names of all the output files
        self.outfile_data = self.datapath_db.fetch(self.params['outfile_data'],
                                              intend_write=True,
                                              purpose="binned WiggleZ data")

        self.outfile_delta_data = self.datapath_db.fetch(\
                                        self.params['outfile_deltadata'],
                                        intend_write=True,
                                        purpose="binned WiggleZ data delta")

        self.outfile_selection = \
                     self.datapath_db.fetch(self.params['outfile_selection'],
                                        intend_write=True,
                                        purpose="WiggleZ selection function")

        self.outfile_separable = \
                     self.datapath_db.fetch(self.params['outfile_separable'],
                                        intend_write=True,
                                        purpose="WiggleZ separable sel func")

        self.outfile_mock = self.datapath_db.fetch(self.params['outfile_mock'],
                                        intend_write=True,
                                        purpose="WiggleZ binned mock",
                                        silent=True)

        self.outfile_delta_mock = self.datapath_db.fetch(\
                                        self.params['outfile_deltamock'],
                                        intend_write=True,
                                        purpose="WiggleZ binned delta mock",
                                        silent=True)

        # gather axis information from the template file
        self.template_map = \
            algebra.make_vect(algebra.load(self.params['template_file']))

        self.freq_axis = self.template_map.get_axis('freq')
        self.ra_axis = self.template_map.get_axis('ra')
        self.dec_axis = self.template_map.get_axis('dec')

        # placeholders for data products
        self.realmap_binning = None
        self.selection_function = None
        self.separable_selection = None

    def execute(self, processes):
        #np.set_printoptions(threshold=np.nan)
        print "finding the binned data"
        self.realmap()
        print "finding the binned mock and selection function"
        self.selection()
        print "finding the separable form of the selection"
        self.separable()
        print "finding the optical overdensity"
        self.delta()

    def realmap(self):
        """bin the real WiggleZ catalog"""
        self.realmap_binning = bin_catalog_file(self.infile_data,
                                                self.freq_axis,
                                                self.ra_axis, self.dec_axis,
                                                skip_header=1)

        map_wigglez = algebra.make_vect(self.realmap_binning,
                                        axis_names=('freq', 'ra', 'dec'))

        map_wigglez.copy_axis_info(self.template_map)
        algebra.save(self.outfile_data, map_wigglez)

        return

    def selection(self):
        """bin the mock catalogs"""
        self.selection_function = np.zeros(self.template_map.shape)

        for mockindex in self.infile_mock[0]:
            print mockindex
            mockfile = self.infile_mock[1][mockindex]
            mock_binning = bin_catalog_file(mockfile, self.freq_axis,
                                            self.ra_axis, self.dec_axis,
                                            skip_header=1, mock=True)

            self.selection_function += mock_binning

            # if this binned mock catalog should be saved
            if mockindex in self.outfile_mock[0]:
                print "mock", self.outfile_mock[1][mockindex]
                map_wigglez_mock = algebra.make_vect(mock_binning,
                                    axis_names=('freq', 'ra', 'dec'))

                map_wigglez_mock.copy_axis_info(self.template_map)

                algebra.save(self.outfile_mock[1][mockindex], map_wigglez_mock)

        # adding the real map back to the selection function is a kludge which
        # ensures the selection function is not zero where there is real data
        # (limit of insufficient mocks)
        self.selection_function += self.realmap_binning
        self.selection_function /= float(len(self.infile_mock[0]) + 1)
        print np.mean(self.selection_function)

        map_wigglez_selection = algebra.make_vect(self.selection_function,
                                                  axis_names=('freq', 'ra', 'dec'))

        map_wigglez_selection.copy_axis_info(self.template_map)

        algebra.save(self.outfile_selection, map_wigglez_selection)

        return

    def separable(self):
        # now assume separability of the selection function
        spatial_selection = np.sum(self.selection_function, axis=0)

        freq_selection = np.apply_over_axes(np.sum,
                                            self.selection_function, [1, 2])

        self.separable_selection = (freq_selection * spatial_selection)

        self.separable_selection /= np.sum(freq_selection.flatten())

        map_wigglez_separable = algebra.make_vect(self.separable_selection,
                                                  axis_names=('freq', 'ra', 'dec'))

        map_wigglez_separable.copy_axis_info(self.template_map)

        algebra.save(self.outfile_separable, map_wigglez_separable)

        return


    def produce_delta_map(self, optical_file, optical_selection_file):
        map_optical = algebra.make_vect(algebra.load(optical_file))
        map_nbar = algebra.make_vect(algebra.load(optical_selection_file))

        old_settings = np.seterr(invalid="ignore", under="ignore")
        map_delta = map_optical / map_nbar - 1.
        np.seterr(**old_settings)

        # TODO: also consider setting the nbar to zero outside of galaxies?
        map_delta[np.isinf(map_delta)] = 0.
        map_delta[np.isnan(map_delta)] = 0.
        # if e.g. nbar is zero, then set the point as if there were no galaxies
        # downstream, nbar=0 should coincide with zero weight anyway
        #map_delta[np.isinf(map_delta)] = -1.
        #map_delta[np.isnan(map_delta)] = -1.

        return map_delta

    def delta(self):
        """find the overdensity using a separable selection function"""
        delta_data = self.produce_delta_map(self.outfile_data,
                                            self.outfile_separable)

        algebra.save(self.outfile_delta_data, delta_data)

        for mockindex in self.outfile_mock[0]:
            print "mock delta", mockindex
            mockinfile = self.outfile_mock[1][mockindex]
            mockoutfile = self.outfile_delta_mock[1][mockindex]

            delta_mock = self.produce_delta_map(mockinfile,
                                                self.outfile_separable)

            algebra.save(mockoutfile, delta_mock)
