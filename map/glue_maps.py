"""Pipeline module for gluing multiple maps for adjacent bands together."""

import glob

import numpy as np

import core.algebra as al
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce

params_init = {'input_root' : './',
               'polarizations' : ('I',),
               'output_root' : './',
               'bands' : ()
               }
prefix = 'gm_'

class GlueMaps(object):

    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix=prefix, feedback=feedback)
        self.feedback = feedback

    def execute(self, nprocesses):
        params = self.params
        # Make parent directory and write parameter file.
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        in_root = params['input_root']        
        # Figure out what the band names are.
        bands = params['bands']
        if not bands:
            map_files = glob.glob(in_root + pol_str + "_*.npy")
            bands = []
            root_len = len(in_root)
            for file_name in map_files:
                bands.append(file_name[root_len:-4])
        # Loop over polarizations.
        for pol_str in params['polarizations']:
            # Read in all the maps to be glued.
            maps = []
            for band in bands:
                band_map_fname = (in_root + pol_str + "_" +
                              repr(band) + '.npy')
                if self.feedback > 1:
                    print "Read using map: " + band_map_fname
                band_map = al.load(band_map_fname)
                band_map = al.make_vect(band_map)
                if band_map.axes != ('freq', 'ra', 'dec') :
                    msg = ("Expeced maps to have axes ('freq',"
                           "'ra', 'dec'), but it has axes: "
                           + str(band_map.axes))
                    raise ce.DataError(msg)
                maps.append(band_map)
            # Now glue them together.
            out_map = glue(maps)
            out_fname = (params['output_root']
                         + pol_str + "_" + "all" + '.npy')
            if self.feedback > 1:
                print "Writing glued map to: " + out_fname
            al.save(out_fname, out_map)


def glue(maps):
    """Function that glues a set of maps together into one map.

    Gluing is done along frequency axis
    """
    
    # Check the non-frequency axis to make sure all the maps are the same.
    angular_shape = maps[0].shape[1:]
    ra_centre = maps[0].info["ra_centre"]
    dec_centre = maps[0].info["dec_centre"]
    ra_delta = maps[0].info["ra_delta"]
    dec_delta = maps[0].info["dec_delta"]
    for map in maps:
        if map.shape[1:] != angular_shape:
            msg = ("Angular shapes don't match up. " + str(map.shape[1:])
                   + " vs " + str(angular_shape))
            raise ce.DataError(msg)
        if (not np.allclose(map.info["ra_centre"], ra_centre)
            or not np.allclose(map.info["dec_centre"], dec_centre)
            or not np.allclose(map.info["ra_delta"], ra_delta)
            or not np.allclose(map.info["dec_delta"], dec_delta)):
            msg = "Angular information doesn't match up."
            raise ce.DataError(msg)
    # Now check that the frequecy axes all have the same spacing and get the
    # centres.
    delta = maps[0].info["freq_delta"]
    centres = []   # In units of delta.
    # Size of the frequency axis.
    size = 0
    for map in maps:
        if not np.allclose(map.info["freq_delta"], delta):
            msg = "Frequency spacing not all the same."
            raise ce.DataError(msg)
        centres.append(map.info["freq_centre"] / delta)
        size += map.shape[0]
    # Sort the map by thier centre.
    inds = np.argsort(centres)
    # Allocate memory for the output map.
    out_map = np.empty((size,) + angular_shape, dtype=float)
    out_map = al.make_vect(out_map, axis_names=('freq', 'ra', 'dec'))
    out_map.copy_axis_info(maps[0])
    # Loop over the maps and put them in the larger map.
    so_far = 0  # frequency axis used so far.
    for index in inds:
        map = maps[index]
        s = map.shape[0]
        freq = map.get_axis('freq')
        # Check if the centre of the total band is in this map and if so, set
        # the centre.
        if size//2 >= so_far and size//2 < so_far + s:
            centre = freq[size//2 - so_far]
            out_map.set_axis_info('freq', centre, delta)
        # Check that the space between the last frequency of the last map and
        # the first frequency of the first map is delta.
        if so_far > 0:
            if not np.allclose(freq[0], last_freq + delta):
                raise ce.DataError("Bands don't line up correctly.")
        last_freq = freq[-1]
        # Finally copy the map into the conglomerate map.
        out_map[so_far:so_far + s,:,:] = map
        so_far += s
    if so_far != size:
        raise RuntimeError("Something went horrible wrong.")
    return out_map





