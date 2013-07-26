import numpy as np
import core.algebra as al
import copy as cp


def preprocess(map, info=None):
    r"""Ensures that the input map has the correct structure and information
    for the kappa estimator
    
    Input parameters :
        - map : ndarray or vect, input overdensity mapfor the estimator
        - info : dictionnary, information on box axes. 
                 must to contain the same keys as vect info
    Output :
        - map_out : vect, structure adapted to the estimator
    """
    try:
        if info is None:
            info = cp.deepcopy(map.info)

        if not info['axes'] == ('z', 'x', 'y'):
            valid_info = {'axes': ('z', 'x', 'y'), 'type': 'vect'}
            for i in range(3):
                axis = info['axes'][i]
                valid_axis = valid_info['axes'][i]
                valid_info[valid_axis+'_centre'] = info[axis+'_centre']
                valid_info[valid_axis+'_delta'] = info[axis+'_delta']
            info = valid_info
        map_out = al.make_vect(map)
        map_out.info = info
        return map_out

    except KeyError as e:
        print "Error : axis names are not defined"
    except TypeError as e:
        print "Error : axis names are None"
    except AttributeError as e:
        print "Error : no information provided"


def set_info(map, size):
    r"""This function sets the information of a box 
    so that is can be used by the estimator
    
    Input parameters:
        - map : ndarray, input overdensity mapfor the estimator
        - size : tuple of length 3, physical length of each size of the box
    Output :
        - map_out : vect, structure adapted to the estimator
    """
    map_out = al.make_vect(map, ('z', 'x', 'y'))
    map_out.set_axis_info('z', 0., size[0]/map.shape[0])
    map_out.set_axis_info('x', 0., size[1]/map.shape[1])
    map_out.set_axis_info('y', 0., size[2]/map.shape[2])
    return map_out
