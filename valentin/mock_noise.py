import numpy as np
import core.algebra as al
import valentin.map_build as mb
import valentin.aggregate_outputs as ao
from valentin.my_utils import own_fftconvolve
from map.physical_gridding import *

import time
import shelve
import random

def load_2df(filename=None):
	t = time.time()
	if filename==None:
		path = '/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5/'
		name = 'real_map_2df_delta.npy'
		filename = path+name
	map = al.load(filename)
	map = al.make_vect(map)
	map_xyz,info = physical_grid_largeangle(map,refinement=1,order=1)
#	sample = map_xyz[183:255,400:912,22:150]
	sample = map_xyz
	sample = al.make_vect(sample,('z','x','y'))
	sample.set_axis_info('z',info['freq_centre'],info['freq_delta'])
	sample.set_axis_info('x',info['ra_centre'],info['ra_delta'])
	sample.set_axis_info('y',info['dec_centre'],info['dec_delta'])
	t = time.time() - t
	print t
	return sample
	
def rescale(map, epsilon=1e-5):
	return map*(-1+epsilon)/map.min()	
	
def init_w():
	delta = al.load('/mnt/raid-project/gmrt/goblot/2df_delta/large_angle_ref_1.npy')
	delta = al.make_vect(delta)
#	data = shelve.open('/mnt/raid-project/gmrt/goblot/2df_mock/mock_0.hdf5','r',protocol=-1)
#	delta = data['map']
#	data.close()
	tmp = mb.MapBuilder(delta,sigma_gaussian=1.25,full=False)
	return tmp.w/tmp.w.sum()

def log_smoothing(ind,w):
	path = '/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_full_selection_10000mock/'
	if ind<10:
		name = 'mock_map_2df_delta_00%s.npy'%ind
	elif ind<100:
		name = 'mock_map_2df_delta_0%s.npy'%ind
	else:
		name = 'mock_map_2df_delta_%s.npy'%ind
	
	delta = load_2df(path+name)
#	filename = '/mnt/raid-project/gmrt/goblot/2df_mock/rs_mock_%s.hdf5'%ind
#	data = shelve.open(filename,'r',protocol=-1)
#	delta = data['map']
#	data.close()
	delta = rescale(delta)
#	del delta
	delta = own_fftconvolve(delta,w)
	m = delta.min()
	assert m>-1, "Error in delta min value: %f"%m
	del w
	delta = np.log(1+delta)
	
	filename = '/mnt/raid-project/gmrt/goblot/2df_mock/mock_%s.hdf5'%ind
	output = shelve.open(filename,'n',protocol=-1)
	output['map'] = delta
	output.close()
	return filename



def save_step(delta,ind):
	filename = '/mnt/raid-project/gmrt/goblot/2df_mock/mock_%s.hdf5'%ind
	output = shelve.open(filename,'n',protocol=-1)
	output['map'] = delta
	output.close()
	return filename
	
def load_step(ind, prefix=''):
	filename = '/mnt/raid-project/gmrt/goblot/2df_mock/'+prefix+'mock_%s.hdf5'%ind
	data = shelve.open(filename,'r',protocol=-1)
	delta = data['map']
	data.close()
	return delta
	
#Convert from freq,ra,dec to z,x,y, max 20 cpu
def step0(ind):
	path = '/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_full_selection_10000mock/'
	if ind<10:
		name = 'mock_map_2df_delta_00%s.npy'%ind
	elif ind<100:
		name = 'mock_map_2df_delta_0%s.npy'%ind
	else:
		name = 'mock_map_2df_delta_%s.npy'%ind
	delta = load_2df(path+name)
	return save_step(delta,ind)

#Rescale, max 20 cpu
def step1(ind):
	delta = load_step(ind)
	delta = rescale(delta)
	return save_step(delta,ind)

#Convolve with gaussian window, max 8 cpu
def step2(ind,w):
#	if ind<=3 : time.sleep(30)
	delta = load_step(ind,prefix = 'rs_')
#	filename = '/mnt/raid-project/gmrt/goblot/2df_mock/rs_mock_%s.hdf5'%ind
#	data = shelve.open(filename,'r',protocol=-1)
#	delta = data['map']
#	data.close()
	print "Mock",ind,"loaded"	
	delta = own_fftconvolve(delta, w)
	return save_step(delta,ind)

#Do the log smoothing, max 20 cpu
def step3(ind):
	delta = load_step(ind)
	print "Mock",ind,"loaded"
	delta = np.log(1+delta)
	return save_step(delta,ind)
	
	
def run():
	w = None
#	w = np.ones((261, 1384, 177))
#	w = init_w()
	t = time.time()
#	mock_agg = ao.AggregateOutputs('valentin.mock_noise.log_smoothing')
	mock_agg = {
		0 : ao.AggregateOutputs('valentin.mock_noise.step0'),
		1 : ao.AggregateOutputs('valentin.mock_noise.step1'),
#		2 : ao.AggregateOutputs('valentin.mock_noise.step2'),
#		3 : ao.AggregateOutputs('valentin.mock_noise.step3')
	}
	for s in mock_agg:
		print "----------------------------------------------------"
		print "Step",s
		for i in range(1000):
			execute_key = i
#		mock_agg.execute(i,w,execute_key=execute_key)
			if s==2: mock_agg[s].execute(i,w,execute_key=execute_key)
			else: mock_agg[s].execute(i,execute_key=execute_key)
#	mock_agg.multiprocess_stack(filename='/mnt/raid-project/gmrt/goblot/mock_catalogue.hdf5',save_cpu=None,ncpu=8)
		ncpu = 7 if s==2 else 15
		mock_agg[s].multiprocess_stack(filename='/mnt/raid-project/gmrt/goblot/mock_catalogue.hdf5',save_cpu=None,ncpu=ncpu)
	t = time.time() - t
	print 'Done :',t

	
def part_std(ind):
	t = time.time()
#	i = ind%2
#	k = (ind>>1)%2
#	j = ind>>2
#	l_min = 261/2*i
#	m_min = 1384/16*j
#	n_min = 177/2*k
#	l_max = 261/2+1 if i==1 else 261/2
#	m_max = 1384/16+1384%16 if j==15 else 1384/16
#	n_max = 177/2+1 if k==1 else 177/2
#	bloc = np.zeros((100,l_max,m_max,n_max))
	m_min = 1384/16*ind
	m = 1384/16+1384%16 if ind==15 else 1384/16
	m_max = m_min + m
	bloc = np.zeros((100,261,m,177))
	for mock in range(100):
		data = shelve.open('/mnt/raid-project/gmrt/goblot/2df_mock/mock_%s.hdf5'%mock,'r',protocol=-1)
#		bloc[mock] = data['map'][l_min:l_min+l_max,m_min:m_min+m_max,n_min:n_min+n_max]
		bloc[mock] = data['map'][:,m_min:m_max,:]
		data.close()
		print "Bloc",ind,"of mock",mock,"extracted"
	t = time.time() - t
	print t
	return np.std(bloc,axis=0)


def full_std():
	t = time.time()
	noise = ao.AggregateOutputs('valentin.mock_noise.part_std')
	for bloc in range(16):
		execute_key = bloc
		noise.execute(bloc,execute_key=execute_key)
	noise.multiprocess_stack(filename='/mnt/raid-project/gmrt/goblot/test_mock_noise.hdf5',save_cpu=None,ncpu=4)
	t = time.time() - t
	print 'Done :',t

		
def recombine():
	blocs = shelve.open('/mnt/raid-project/gmrt/goblot/mock_noise.hdf5','r',protocol=-1)
	noise = np.zeros((261,1384,177))
	for b in range(64):
		i = b%2
		k = (b>>1)%2
		j = b>>2
		l_min = 261/2*i
		m_min = 1384/16*j
		n_min = 177/2*k
		l_max = 261/2+1 if i==1 else 261/2
		m_max = 1384/16+1384%16 if j==15 else 1384/16
		n_max = 177/2+1 if k==1 else 177/2
		noise[l_min:l_min+l_max,m_min:m_min+m_max,n_min:n_min+n_max] = blocs[b]
	return noise

def my_shuffle(a):
	random.shuffle(a)
	return a

if __name__=='__main__':
	full_std()