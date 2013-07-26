import numpy as np
import h5py
import core.algebra as al
from utils import binning

def JD_cat_binning(path,n_pix,size,calc_delta=False,data='num_galaxies'):
	z_edges = np.arange(n_pix[0]+1)*size[0]/(n_pix[0])
	x_edges = np.arange(n_pix[1]+1)*size[1]/(n_pix[1])
	y_edges = np.arange(n_pix[2]+1)*size[2]/(n_pix[2])
#	return z_edges,x_edges,y_edges
	map = np.zeros(n_pix)
	for i in range(8):
		print 'Step',i,'...',
		f = h5py.File(path+'0.800halo_catalog000%s.hdf5'%i,'r')
		pos = f['Positions']
		num = f[data][data][:]
		if (pos['z'][:].min())<size[0]:
			sample = np.ones((len(num),3))
			sample[:,0] = pos['z'][:]
			sample[:,1] = pos['x'][:]
			sample[:,2] = pos['y'][:]
			f.close()
			map += binning.histogram3d(sample, z_edges, x_edges, y_edges, weight=num)
		else: 
			f.close()
		print 'done'
	if calc_delta:
		mn = map.mean()
		map = map/mn - 1.
	map = al.make_vect(map,('z','x','y'))
	map.set_axis_info('z',size[0]/2.,size[0]/n_pix[0])
	map.set_axis_info('x',size[1]/2.,size[1]/n_pix[1])
	map.set_axis_info('y',size[2]/2.,size[2]/n_pix[2])
	return map