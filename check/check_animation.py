#! /usr/bin/env python 

import numpy as np
from core import algebra
from utils import binning
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

map_root = '/cita/h/home-2/ycli/'
#map_file = 'real_map_2df'
map_file = 'mock_map_2df_097'
#map_file = 'fir_15hr_41-80_avg_fdgp_new_clean_map_I_800'

map = algebra.make_vect(algebra.load(map_root + map_file + '.npy'))
ra = map.get_axis('ra')
dec = map.get_axis('dec')
freq = map.get_axis('freq')
ra_edges = binning.find_edges(ra)
dec_edges = binning.find_edges(dec)

map = np.ma.array(map)
map[map==0] = np.ma.masked

# Set up formatting for the movie files
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=18)

fig = plt.figure(figsize=(10,7))
ims = []
for i in range(len(freq)):
    ims.append((plt.pcolormesh(ra_edges, dec_edges, map[i], 
                label='freq: %7.3f MHz'%(freq[i]/1.e6)),
                plt.text(ra[-1], dec[-4], 
                '%s 2df freq: %7.3f MHz'%(map_file, freq[i]/1.e6)),
                #plt.title('2df freq: %7.3f MHz'%(freq[i]/1.e6)),
                ))
    #plt.legend(loc='upper centre', frameon=False, bbox_to_anchor = (0.5, 0.5))

plt.xlim(xmin=ra_edges.min(), xmax=ra_edges.max())
plt.ylim(ymin=dec_edges.min(), ymax=dec_edges.max())
plt.xlabel('RA [degree]')
plt.ylabel('DEc [degree]')
#plt.title('2df freq: %7.3f MHz'%(freq[i]/1.e6))

im_ani = animation.ArtistAnimation(fig, ims,interval=500)

im_ani.save('./png/%s.mp4'%map_file)
#plt.show()
