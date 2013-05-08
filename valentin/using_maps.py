import numpy as np
import matplotlib.pyplot as plt

import core.algebra as al

plt.ion()

file_name = "/mnt/raid-project/gmrt/eswitzer/GBT/cleaned_maps/GBT_15hr_map_autopaper/sec_C_cleaned_clean_map_I_with_B_20modes.npy"
#file_name = "/mnt/raid-project/gmrt/eswitzer/GBT/cleaned_maps/GBT_15hr_map_autopaper/combined_clean_map_20modes.npy"

map = al.load(file_name)
map = al.make_vect(map)



plt.imshow(map[50,:,:].T)

print map.axes
print map.shape
print len(map.get_axis('freq'))

#sec = ['A','B','C','D']
#map_moy = np.zeros(map[50,:,:].shape)
#for X in sec :
#	for Y in sec :
#		if X!=Y:
#			file_name_tmp = "/mnt/raid-project/gmrt/eswitzer/GBT/cleaned_maps/GBT_15hr_map_autopaper/" \
#			"sec_"+ X + "_cleaned_clean_map_I_with_"+ Y +"_10modes.npy"
#			map_tmp = al.load(file_name_tmp)
#			map_tmp = al.make_vect(map_tmp)
#			map_moy[:] += map_tmp[50,:,:]
#map_moy[:] /= 8
	    
def affiche(n) :
	plt.clf()
	plt.imshow(map[n,8:70,8:35].T) #, vmin = -0.65, vmax = 0.65)
	plt.colorbar()        
