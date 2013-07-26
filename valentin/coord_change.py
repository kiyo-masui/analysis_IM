import numpy as np
import scipy as sp
import scipy.ndimage


def convert_info(info):
	l,m,n = info['shape']
	z_delta = info['z_full']/l
	x_delta = info['x_full']/m
	y_delta = info['y_full']/n
	info['z_delta'] = z_delta
	info['x_delta'] = x_delta
	info['y_delta'] = y_delta
	return info

def ind2coord(i,j,k, info_in):
	l,m,n = info_in['shape']
	z_i = (i - l/2)*info_in['z_delta'] + info_in['z_centre']
	x_j = (j - m/2)*info_in['x_delta'] + info_in['x_centre']
	y_k = (k - n/2)*info_in['y_delta'] + info_in['y_centre']
	return (z_i,x_j,y_k)

def change_box(map_in, info_in, info_out,order=1):
	l,m,n = info_out['shape']
	li,mi,ni = map_in.shape
	
	z_axis_out = (np.arange(l)-l/2)*info_out['z_delta'] + info_out['z_centre'] - info_in['z_centre']
	x_axis_out = (np.arange(m)-m/2)*info_out['x_delta'] + info_out['x_centre'] - info_in['x_centre']
	y_axis_out = (np.arange(n)-n/2)*info_out['y_delta'] + info_out['y_centre'] - info_in['y_centre']
	
	grid_z = z_axis_out[:,None,None]*np.ones((1,m,n))
	grid_x = x_axis_out[None,:,None]*np.ones((l,1,n))
	grid_y = y_axis_out[None,None,:]*np.ones((l,m,1))
	
	interpol_grid = np.zeros((3,l,m,n))
	interpol_grid[0] = grid_z/info_in['z_delta'] + li/2
	interpol_grid[1] = grid_x/info_in['x_delta'] + mi/2
	interpol_grid[2] = grid_y/info_in['y_delta'] + ni/2

	map_out = sp.ndimage.map_coordinates(map_in, interpol_grid, order=order)
	return map_out