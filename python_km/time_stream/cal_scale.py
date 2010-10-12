#!/usr/bin/python
"""Module that put data in units of cal temperture and subtracts median."""

from core import data_block, fitsGBT

def scale_by_cal(Data) :
    """Puts all data in units of the cal temperature."""
    
    # Here we check the polarizations
    xx_ind = 0
    yy_ind = 3
    xy_inds = [1,2]
    on_ind = 0
    off_ind = 1

    #cal_med_xx = 
