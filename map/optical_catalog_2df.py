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
import scipy as sp
import shelve
import random
import copy
import ephem
from core import algebra
from core import constants as cc
from utils import binning
from utils import data_paths
from utils import cosmology as cosmo
from utils import units
from kiyopy import parse_ini
from map import physical_gridding
# TODO: make better parameter passing for catalog binning

def estimate_bias(apparent_mag, redshift, type='blue', absolute_mag_star=-19.66):
    '''
        absolute_mag_star comes from paper Hawkins 2003 MNRAS 346, 78
        bias calculation comes from paper Cole 2005 MNRAS 362, 505
    '''

    cosmology = cosmo.Cosmology()
    distance = cosmology.luminosity_distance(redshift)
    absolute_mag = apparent_mag - 5. * np.log10(distance) + 2.5
    if type=='red':
        bias = 1.3 * ( 0.85 + 0.15 * (absolute_mag - absolute_mag_star))
    else: 
        bias = 0.9 * ( 0.85 + 0.15 * (absolute_mag - absolute_mag_star))
    return bias


def bin_catalog_data(catalog, freq_axis, ra_axis,
                     dec_axis, debug=False, use_histogramdd=False, 
                     get_bias = False):
    """
    bin catalog data onto a grid in RA, Dec, and frequency
    This currently assumes that all of the axes are uniformly spaced
    """
    #catalog_frequencies = cc.freq_21cm_MHz * 1.e6 / (1 + catalog['z'])
    num_catalog = catalog.size
    sample = np.zeros((num_catalog, 3))
    sample[:, 0] = catalog['z']
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
        print np.max(sample, axis=0)
        print np.min(sample, axis=0)

    if use_histogramdd:
        count_cube, edges = np.histogramdd(sample, bins=[freq_edges,
                                                         ra_edges, dec_edges])
        print edges
    else:
        count_cube = binning.histogram3d(sample, freq_edges, ra_edges, dec_edges)
        if get_bias:
            bias = estimate_bias(catalog['mag'], catalog['z'])
            print bias
            bias_cube = binning.histogram3d(sample, freq_edges, ra_edges, 
                                            dec_edges, weight=bias)
            old_settings = np.seterr(invalid='ignore', divide='ignore')
            bias_cube /= count_cube
            bias_cube[np.isnan(bias_cube)] = 0.
            bias_cube[np.isinf(bias_cube)] = 0.
            np.seterr(**old_settings)

    if get_bias:
        return count_cube, bias_cube
    else:
        return count_cube
def convert_to_physical_coordinate(freq, ra, dec, ra_fact):
    
    ra = np.arcsin(ra_fact*np.sin(ra*np.pi/180.))*180./np.pi

    freq = freq/1.e6
    cosmology = cosmo.Cosmology()
    z = units.nu21 / freq - 1.
    c = cosmology.comoving_distance(z)
    d = cosmology.proper_distance(z)
    y = d * ra * np.pi / 180.
    z = d * dec * np.pi / 180.
    #x = np.sqrt(c*c - y*y - z*z)
    x = c * np.cos(ra*np.pi/180.)

    return x, y, z

def convert_B1950_to_J2000(ra, dec, z, degree_in=False, degree_out=True):

    if degree_in:
        ra = ephem.hours(ra*np.pi/180.)
        dec = ephem.degrees(dec*np.pi/180.)
    else:
        ra = ephem.hours(ra)
        dec = ephem.degrees(dec)
    coor = ephem.Equatorial(ra, dec, epoch=ephem.B1950)
    coor = ephem.Equatorial(coor, epoch=ephem.J2000)

    ra2000 = coor.ra
    dec2000 = coor.dec

    # for test if the 2df have B1950
    #ra2000 = ra
    #dec2000 = dec

    #print z, 
    #print " | ra: ",ra," --> ",ra2000," [",(ra-ra2000)*180./np.pi,"] | ",
    #print "dec: ",dec," --> ",dec2000," [",(dec-dec2000)*180./np.pi,"]"

    if degree_out:
        ra2000 = ra2000*180./np.pi
        dec2000 = dec2000*180./np.pi

    return ra2000, dec2000, z

def bin_catalog_file(filename, freq_axis, ra_axis, dec_axis, 
                     physical_cube=False, ra_centre=None, dec_centre=None,
                     skip_header=None, debug=False, mock=False, 
                     Negative_RA=False):
    """Bin the catalog given in `filename` using the given freq, ra, and dec
    axes.
    """

    # read the WiggleZ catalog and convert redshift axis to frequency
    if mock:
        ndtype = [('RA', float), ('Dec', float), ('z', float),
                  ('mag', float), ('wsel', float), ('compl', float),
                  ('selz', float), ('selz_mu', float), ('bjlim', float)]
    else:
        ndtype = [('RA', float), ('Dec', float), ('z', float),
                  ('mag', float), ('wsel', float), ('compl', float),
                  ('selz', float), ('selz_mu', float), ('bjlim', float), 
                  ('serial', int)]
        #ndtype = [('RA', float), ('Dec', float), ('mag', float), ('z', float)]

    # TODO: numpy seems to be an old version that does not have the skip_header
    # argument here! skiprows is identical
    catalog = np.genfromtxt(filename, dtype=ndtype, skiprows=skip_header)

    convert_B1950_to_J2000_vect = np.vectorize(convert_B1950_to_J2000)
    catalog['RA'], catalog['Dec'], catalog['z'] =\
        convert_B1950_to_J2000_vect(catalog['RA'], catalog['Dec'], catalog['z'])

    #catalog['RA']  = catalog['RA']*180./np.pi
    #catalog['Dec'] = catalog['Dec']*180./np.pi
    catalog['z']   = cc.freq_21cm_MHz * 1.e6 / (1. + catalog['z'])

    # change the RA range to -180 ~ 180
    if np.any(ra_axis < 0) and not Negative_RA:
        print "Should set Negative_RA = True"
        exit()
    if Negative_RA:
        catalog['RA'][catalog['RA']>180.] -= 360.

    if ra_centre != None:
        catalog['RA'] -= ra_centre
    if dec_centre != None:
        catalog['Dec'] -= dec_centre
    if physical_cube:
        catalog['z'], catalog['RA'], catalog['Dec'] = convert_to_physical_coordinate(
            catalog['z'], catalog['RA'], catalog['Dec'], sp.cos(sp.pi*dec_centre/180.0))
        print catalog['z'].min(), catalog['z'].max()

    if debug:
        print filename + ": " + repr(catalog.dtype.names) + \
              ", n_records = " + repr(catalog.size)

    if mock:
        return bin_catalog_data(catalog, freq_axis, ra_axis, dec_axis, 
                                debug=debug)
    else:
        return bin_catalog_data(catalog, freq_axis, ra_axis, dec_axis, 
                                debug=debug, get_bias=False)

bin2dfparams_init = {
        "infile_data": "/Users/ycli/DATA/2df/catalogue/real_catalogue_2df.out",
        "infile_mock": "/Users/ycli/DATA/2df/catalogue/mock_catalogue_2df_%03d.out",
        "outfile_data": "/Users/ycli/DATA/2df/map/real_map_2df.npy",
        "outfile_bias": "/Users/ycli/DATA/2df/map/bias_map_2df.npy",
        "outfile_mock": "/Users/ycli/DATA/2df/map/mock_map_2df_%03d.npy",
        "outfile_deltadata": "/Users/ycli/DATA/2df/map/real_map_2df_delta.npy",
        "outfile_deltamock": "/Users/ycli/DATA/2df/map/mock_map_2df_delta_%03d.npy",
        "outfile_selection": "/Users/ycli/DATA/2df/map/sele_map_2df.npy",
        "outfile_separable": "/Users/ycli/DATA/2df/map/sele_map_2df_separable.npy",
        "template_file": "/Users/ycli/DATA/2df/tempfile",
        "physical_cube": False,
        "mock_number": 10,
        "mock_alpha": 1,
        "negative_RA" : False,
        }
bin2dfprefix = 'bin2df_'


class Bin2dF(object):
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        #self.datapath_db = data_paths.DataPath()

        if parameter_file:
            self.params = parse_ini.parse(parameter_file, bin2dfparams_init,
                                          prefix=bin2dfprefix)

        # gather names of the input catalogs
        self.infile_data = self.params['infile_data']

        self.infile_mock = self.params['infile_mock']

        # gather names of all the output files
        self.outfile_data = self.params['outfile_data']

        self.outfile_bias = self.params['outfile_bias']

        self.outfile_delta_data = self.params['outfile_deltadata']

        self.outfile_selection  = self.params['outfile_selection']

        self.outfile_separable  = self.params['outfile_separable']

        self.outfile_mock       = self.params['outfile_mock']

        self.outfile_delta_mock = self.params['outfile_deltamock']

        # gather axis information from the template file
        self.template_map = \
            algebra.make_vect(algebra.load(self.params['template_file']))

        if self.params['physical_cube']:

            self.ra_centre = self.template_map.info['ra_centre']
            self.dec_centre = self.template_map.info['dec_centre']
            self.physical_cube = True

            #self.template_map = physical_gridding.physical_grid_largeangle(
            #    self.template_map, refinement=1, pad=5, order=2, infor_only=True)
            #self.template_map = physical_gridding.physical_grid(
            #    self.template_map, refinement=1, pad=5, order=2, info_only=True)
            self.template_map = physical_gridding.physical_grid_largeangle(
                self.template_map, refinement=1, pad=5, order=2, info_only=True, 
                cube_force=1.)
        else:
            self.ra_centre = None
            self.dec_centre = None
            self.physical_cube = False


        self.freq_axis = self.template_map.get_axis('freq')
        self.ra_axis = self.template_map.get_axis('ra')
        self.dec_axis = self.template_map.get_axis('dec')

        # placeholders for data products
        self.realmap_binning = None
        self.selection_function = None
        self.separable_selection = None

    def execute(self, processes):
        pass
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
        #self.realmap_binning, self.biasmap_binning\
        self.realmap_binning\
            = bin_catalog_file(self.infile_data, self.freq_axis, self.ra_axis,
                               self.dec_axis, skip_header=1, debug=False,
                               physical_cube=self.physical_cube, 
                               ra_centre=self.ra_centre, 
                               dec_centre=self.dec_centre,
                               Negative_RA=self.params['negative_RA'])

        map_2df = algebra.make_vect(self.realmap_binning,
                                        axis_names=('freq', 'ra', 'dec'))
        map_2df.copy_axis_info(self.template_map)
        algebra.save(self.outfile_data, map_2df)

        #bias_2df = algebra.make_vect(self.biasmap_binning,
        #                                axis_names=('freq', 'ra', 'dec'))
        #bias_2df.copy_axis_info(self.template_map)
        #algebra.save(self.outfile_bias, bias_2df)

        #old_setting = np.seterr(divide='ignore')
        #bias_x_map = self.realmap_binning / self.biasmap_binning
        #bias_x_map[np.isnan(bias_x_map)] = 0.
        #bias_x_map[np.isinf(bias_x_map)] = 0.
        #np.seterr(**old_setting)
        #bias_x_map = algebra.make_vect(bias_x_map, axis_names=('freq', 'ra', 'dec'))
        #bias_x_map.copy_axis_info(self.template_map)
        #algebra.save(self.outfile_data, bias_x_map)

        return

    def selection(self):
        """bin the mock catalogs"""
        self.selection_function = np.zeros(self.template_map.shape)

        for mockindex in range(self.params['mock_number']):
            print mockindex
            mockfile = self.infile_mock%mockindex
            mock_binning = bin_catalog_file(mockfile, self.freq_axis,
                                            self.ra_axis, self.dec_axis,
                                            skip_header=1, mock=True,
                                            physical_cube=self.physical_cube, 
                                            ra_centre=self.ra_centre, 
                                            dec_centre=self.dec_centre, 
                                            Negative_RA=self.params['negative_RA'])

            self.selection_function += mock_binning

            # if this binned mock catalog should be saved
            #if mockindex in self.outfile_mock[0]:

            print "mock", self.outfile_mock%mockindex
            map_2df_mock = algebra.make_vect(mock_binning,
                                             axis_names=('freq', 'ra', 'dec'))

            map_2df_mock.copy_axis_info(self.template_map)

            algebra.save(self.outfile_mock%mockindex, map_2df_mock)

        #print "normalized the selection function"
        #self.selection_function *= self.params['mock_alpha']
        # adding the real map back to the selection function is a kludge which
        # ensures the selection function is not zero where there is real data
        # (limit of insufficient mocks)
        self.selection_function += self.realmap_binning
        self.selection_function /= float(self.params['mock_number'] + 1)
        print np.mean(self.selection_function)

        map_2df_selection = algebra.make_vect(self.selection_function,
                                              axis_names=('freq', 'ra', 'dec'))

        map_2df_selection.copy_axis_info(self.template_map)

        algebra.save(self.outfile_selection, map_2df_selection)

        return

    def separable(self):
        # now assume separability of the selection function
        spatial_selection = np.sum(self.selection_function, axis=0)

        freq_selection = np.apply_over_axes(np.sum, self.selection_function, [1, 2])

        self.separable_selection = (freq_selection * spatial_selection)

        self.separable_selection /= np.sum(freq_selection.flatten())

        map_2df_separable = algebra.make_vect(self.separable_selection,
                                                  axis_names=('freq', 'ra', 'dec'))

        map_2df_separable.copy_axis_info(self.template_map)

        algebra.save(self.outfile_separable, map_2df_separable)

        return


    def produce_delta_map(self, optical_file, optical_selection_file, mock=False):
        map_optical = algebra.make_vect(algebra.load(optical_file))
        #if mock:
        #    map_optical *= self.params['mock_alpha']
        map_nbar = algebra.make_vect(algebra.load(optical_selection_file))
        if not mock:
            map_nbar *= self.params['mock_alpha']
            print "normalize selectin fuctino by %f"%self.params['mock_alpha']

        old_settings = np.seterr(invalid="ignore", under="ignore")
        map_delta = map_optical / map_nbar - 1.
        np.seterr(**old_settings)

        # TODO: also consider setting the nbar to zero outside of galaxies?
        #map_delta[np.isinf(map_delta)] = 0.
        #map_delta[np.isnan(map_delta)] = 0.
        # if e.g. nbar is zero, then set the point as if there were no galaxies
        # downstream, nbar=0 should coincide with zero weight anyway
        map_delta[np.isinf(map_delta)] = -1.
        map_delta[np.isnan(map_delta)] = -1.

        return map_delta

    def delta(self):
        #selection_file = self.outfile_separable
        selection_file = self.outfile_selection
        """find the overdensity using a separable selection function"""
        delta_data = self.produce_delta_map(self.outfile_data, selection_file)

        algebra.save(self.outfile_delta_data, delta_data)

        for mockindex in range(self.params['mock_number']):
            print "mock delta", mockindex
            mockinfile = self.outfile_mock%mockindex
            mockoutfile = self.outfile_delta_mock%mockindex

            delta_mock = self.produce_delta_map(mockinfile, selection_file, mock=True)

            algebra.save(mockoutfile, delta_mock)

if __name__=="__main__":
    
    import os

    tempfile = algebra.make_vect(np.ones(shape=(256,256,32)), 
                                 axis_names=('freq', 'ra', 'dec'))
    tempfile.info['ra_delta']  = -0.35
    tempfile.info['dec_delta'] = 0.35
    tempfile.info['ra_centre'] = 7.00
    tempfile.info['dec_centre'] = -29.5
    tempfile.info['freq_delta'] = -1000000.0
    tempfile.info['freq_centre'] = 1221900000.0

    #tempfile = algebra.make_vect(np.ones(shape=(256,128,128)), 
    #                             axis_names=('freq', 'ra', 'dec'))
    #tempfile.info['ra_delta']  = -0.079
    #tempfile.info['dec_delta'] = 0.079
    #tempfile.info['ra_centre'] = 29.0
    #tempfile.info['dec_centre'] = -29.5
    #tempfile.info['freq_delta'] = -1000000.0
    #tempfile.info['freq_centre'] = 1221900000.0

    algebra.save('/mnt/scratch-gl/ycli/2df_catalog/temp/tempfile', tempfile)

    #map_dir = '/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_oneseed_separable/'
    #map_dir = '/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_full_selection_10000mock/'
    map_dir = '/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_parkes_selection_1000mock/'
    if not os.path.exists(map_dir):
        os.makedirs(map_dir)

    bin2dfparams_init = {
        "infile_data": "/mnt/scratch-gl/ycli/2df_catalog/catalog/real_catalogue_2df.out",
        "infile_mock": "/mnt/scratch-gl/ycli/2df_catalog/catalog/mock_catalogue_2df_%04d.out",
        "outfile_data": map_dir + "real_map_2df.npy",
        "outfile_bias": map_dir + "bias_map_2df.npy",
        "outfile_mock": map_dir + "mock_map_2df_%04d.npy",
        "outfile_deltadata": map_dir + "real_map_2df_delta.npy",
        "outfile_deltamock": map_dir + "mock_map_2df_delta_%03d.npy",
        "outfile_selection": map_dir + "sele_map_2df.npy",
        "outfile_separable": map_dir + "sele_map_2df_separable.npy",
        #"template_file": "/mnt/scratch-gl/ycli/2df_catalog/temp/tempfile",
        "template_file": "/mnt/raid-project/gmrt/anderson/first_parkes_pipe/maps/test_allbeams_27n30_10by7_clean_map_I_1315.npy",
        "mock_number": 100,
        "mock_alpha" : 1.,
        }
    
    Bin2dF(params_dict=bin2dfparams_init).execute(2)
    
