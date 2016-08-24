#! /usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import math

WORKROOT = '/mnt/raid-project/gmrt/ycli/ps_result/'
#CASE = ('cros', 'auto')
CASE = ('auto',)

HOUR = 15

#FILELIST = (5, 15, 25, 35, 45, 55, 65, 80)
#SAVENAME = 'bias_%dhr_5-80mode'%(HOUR, )
#MAP = (
#       'GBT_15hr_map_oldcal_legendre_modes_15gwj',
#      )

FILELIST = (80, )
SAVENAME = 'bias_%dhr_80mode_0-75'%(HOUR, )
MAP = (
       'GBT_15hr_map_oldcal_legendre_modes_0gwj',
       'GBT_15hr_map_oldcal_legendre_modes_5gwj',
#      'GBT_15hr_map_oldcal_legendre_modes_10gwj',
        'GBT_15hr_map_oldcal_legendre_modes_15gwj',
#      'GBT_15hr_map_oldcal_legendre_modes_20gwj',
        'GBT_15hr_map_oldcal_legendre_modes_25gwj',
#      'GBT_15hr_map_oldcal_legendre_modes_30gwj',
        'GBT_15hr_map_oldcal_legendre_modes_35gwj',
#      'GBT_15hr_map_oldcal_legendre_modes_40gwj',
        'GBT_15hr_map_oldcal_legendre_modes_45gwj',
#      'GBT_15hr_map_oldcal_legendre_modes_50gwj',
        'GBT_15hr_map_oldcal_legendre_modes_55gwj',
      'GBT_15hr_map_oldcal_legendre_modes_65gwj',
#        'GBT_15hr_map_oldcal_legendre_modes_70gwj',
       'GBT_15hr_map_oldcal_legendre_modes_75gwj',
      )

#FILELIST = (15, 25, 50)
#SAVENAME = 'bias_%dhr_15_25_50mode%s%s'%(HOUR, FLAG, BEAM)

linestyle = ('--', '-')
plt.figure(figsize=(8,8))
for map in MAP:
    for case in CASE:
    	for i in FILELIST:
    		cros_root = WORKROOT + 'bias/%s_%s_%d/'%(case, map, i)
    		print cros_root
    	
    		B = np.load(cros_root+'b_bias.npy')
    		B = (1/B)
    		if case=='cros':
    			B = B**2
    		k = np.load(cros_root+'k_bias.npy')
    	
    		dk = math.sqrt(k[1]/k[0])
    		k = k*dk
    	
    		if case=='cros':
    			lt='--'
    		else: 
    			lt='-'
    		plt.plot(k, B, label='%s %s %02dmodes'%(case, map, i), 
                linewidth=2, linestyle=lt)
	
plt.axhline(y=1, linewidth=2, linestyle='--', color='k')

#plt.semilogx()
plt.loglog()
#plt.xlim(xmin=k.min(), xmax=k.max())
plt.xlim(xmin=0.025, xmax=1.5)
plt.ylim(ymin=1.e-6,ymax=100)
#legendtitle = 'Transfer function\n'
legendtitle = '(P[clean_{map+sim}(map+sim)]-P[clean_{map}(map)])/P[sim]'
plt.legend(loc=1, frameon=False, title=legendtitle)

plt.tick_params(length=6, width=1.)
plt.tick_params(which='minor', length=3, width=1.)

plt.xlabel('k [h/Mpc]')
plt.ylabel('Transfer Function [$T^2$]')

plt.savefig('./png/'+SAVENAME+'.png', formate='png')
plt.savefig('./eps/'+SAVENAME+'.eps', formate='eps')

#plt.show()
