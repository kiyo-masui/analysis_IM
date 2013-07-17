import numpy as np
import core.algebra as al
import copy as cp

def preprocess(map, info=None):
	try:
		if info==None:
			info = cp.deepcopy(map.info)
		
		if not info['axes']==('z','x','y'):
			valid_info = {'axes' : ('z','x','y'), 'type' : 'vect'}
			for i in range(3):
				axis = info['axes'][i]
				valid_axis = valid_info['axes'][i]
				valid_info[valid_axis+'_centre'] = info[axis+'_centre']
				valid_info[valid_axis+'_delta'] = info[axis+'_delta']
			info = valid_info
		map_out = al.make_vect(map)
		map_out.info = info
		return map_out
		
	except KeyError as e :
		print "Error : axis names are not defined"
	except TypeError as e :
		print "Error : axis names are None"
	except AttributeError as e :
		print "Error : no information provided"
		
def set_info(map, size):
	map_out = al.make_vect(map_out,axes=('z','x','y'))
	map_out.set_axis_info('z', 0., size[0]/map.shape[0])
	map_out.set_axis_info('x', 0., size[1]/map.shape[1])
	map_out.set_axis_info('y', 0., size[2]/map.shape[2])
	return map_out