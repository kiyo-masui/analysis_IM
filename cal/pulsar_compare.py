import os
import numpy as np


pulsar_data_root = os.getenv('RAID_PRO') + 'kiyo/beam_fits/pulsar_grid/'

coords_fname = 'coords.txt'

jones_data_root = 'sgrid_'
jones_data_end = '.txt'

n_chan = 256

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
    print name
    print coord
    coord = (coord[0] * az_factor, coord[1])
    new_points[name] = coord
    print coord
    print np.sqrt(coord[0]**2 + coord[1]**2)
    print

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


