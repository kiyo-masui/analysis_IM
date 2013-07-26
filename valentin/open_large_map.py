import numpy as np
import array
import time

#Load Joachim's simulation

t = time.time()

path = '/mnt/raid-project/gmrt/goblot/'
file = open(path+"recon67_0.054den.dat","rb")
data = file.read()

file.close()

a = array.array('f')
a.fromstring(data)

del data

b = np.array(a)
del a
c = b.reshape(1024,1024,1026)
del b

d = np.zeros((1024,1024,1024))
d = c[:,:,:1024]

t = time.time() -t
print "Map loaded in %.2f sec" %t


def low_resolution(map, k=2):
	l,m,n = map.shape
	map_red = map.reshape(l/k,k,m/k,k,n/k,k)
	map_red = np.sum(map_red, axis=5)
	map_red = np.sum(map_red, axis=3)
	map_red = np.sum(map_red, axis=1)
	return map_red/k**3


#Load 2df data

def load_2df(filename=None):
	if filename==None:
		path = '/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5/'
		name = 'real_map_2df_delta.npy'
		filename = path+name
	map = al.load(filename)
	map = al.make_vect(map)
	map_xyz,info = physical_grid(map,refinement=1,order=1)
#	sample = map_xyz[183:255,400:912,22:150]
	sample = map_xyz
	sample = al.make_vect(sample,('z','x','y'))
	sample.set_axis_info('z',info['freq_centre'],info['freq_delta'])
	sample.set_axis_info('x',info['ra_centre'],info['ra_delta'])
	sample.set_axis_info('y',info['dec_centre'],info['dec_delta'])
	return sample