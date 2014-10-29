import os
import numpy as np
from scipy import interpolate


point_source_data_root = os.getenv('RAID_PRO') + 'kiyo/beam_fits/'
beam_file = '3c147/full_beam_3c147.npy'

pulsar_data_root = os.getenv('RAID_PRO') + 'kiyo/beam_fits/pulsar_grid/'

coords_fname = 'coords.txt'

jones_data_root = 'sgrid_'
jones_data_end = '.txt'

n_chan = 256

#### Get all the pulsar data, first the coordinates, then the jones matrices.
# Read in the coordinate data.
f = open(pulsar_data_root + coords_fname, 'r')
lines = f.readlines()
f.close()
# Reorganize it.
points = {}
for line in lines:
    # Chop off the newline character.
    entries = line.split()
    points[entries[0]] = (float(entries[1]), float(entries[2]))

# Figure out the field center.
center = points['C']
az_factor = np.cos(center[1] * np.pi / 180)
points['C'] = (0., 0.)
# Adjust all the azimuths for the cos(el) factor.
new_points = {}
for name, coord in points.iteritems():
    coord = (coord[0] * az_factor, coord[1])
    new_points[name] = coord
points = new_points

# Now read in the Jones matrix data.
jones_points = {}
for name in points.iterkeys():
    fname = pulsar_data_root + jones_data_root + name + jones_data_end
    data = np.loadtxt(fname)
    n_entries = data.shape[0]
    jones = np.zeros((n_chan, 2, 2), dtype=np.complex128)
    for ii in xrange(n_entries):
        jones[data[ii,0],0,0] = data[ii,1] + data[ii,2]*1j
        jones[data[ii,0],0,1] = data[ii,3] + data[ii,4]*1j
        jones[data[ii,0],1,0] = data[ii,5] + data[ii,6]*1j
        jones[data[ii,0],1,1] = data[ii,7] + data[ii,8]*1j
    jones_points[name] = jones

#### Now read in the fits to the point source beam maps.
# Map has 256 spectral channels and is the XX, XY, YX, YY response to a pure I
# source.
beam_map = np.load(point_source_data_root + beam_file)
map_side = 1.0  # Map dimension in degrees.
n_side = beam_map.shape[2]
#if not n_side % 2:
#    # Want an odd number of points so the n_side // 2 point is at the exact
#    # centre.
#    print "Warning, beam map has an even number of points."

#### Now do the comparison by interpolating the beam maps to pulsar pointings.
for location in points.iter_keys():
    coord = points[location]
    jones = jones_points[location]
    beam_skewer = interpolate.interp2d(
