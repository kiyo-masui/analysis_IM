#! /usr/bin/env python 

import numpy as np
import numpy.ma as ma
import copy
import matplotlib.pyplot as plt

class MultibeamCal(object):
    '''
        This mudule is used to read and analysis the mbcal file
    '''

    def __init__(self, mbcal_file, calibrator_flux=14.9, full_mbcal_file=None, 
            moving_beam_list = []):

        self.cal_inf = np.loadtxt(mbcal_file)
        self.calibrator_flux = calibrator_flux

        if full_mbcal_file != None:
            self.read_full_mbcal(full_mbcal_file)
            self.factor = self.get_calibrator_flux(moving_beam_list)/self.cal_inf[:2,:]
        else:
            self.factor = 1.
            

    def get_cal_factor(self):

        self.Jy2K  = np.array([self.cal_inf[6]/self.cal_inf[12],
                               self.cal_inf[7]/self.cal_inf[13]])
        self.Ginv  = np.array([self.calibrator_flux/self.cal_inf[12],
                               self.calibrator_flux/self.cal_inf[13]])
        self.Ginv /= self.factor

    def read_full_mbcal(self, full_mbcal_file):
	    f = open(full_mbcal_file, 'r')
	    data = []
	    bat_n = 0
	    for line in f:
	    	line = line.split()
	    	if line[0]=='##':
	    		continue
	    	elif line[0]=='$$':
	    		if line[1]=='beam':
	    			#print bat_n
	    			#print line 
	    			bat_n = 0
	    	elif line[0]=='BAT':
	    		bat_n += 1
	    	else:
	    		data.append([float(x) for x in line])
	    #print bat_n
	    f.close()
	    
	    self.data = ma.array(data).reshape(-1,13,7)

    def get_calibrator_flux(self, moving_beam_list = []):

        data = copy.deepcopy(self.data)

        for i in moving_beam_list:
            data = self.mask_moving_time(data, i[0], i[1])
    
    	flux = []
    
    	# noise cal off
    	A = data[...,4]/(data[...,3] - data[...,4])
    	B = data[...,6]/(data[...,5] - data[...,6])
    
    	A = self.mask_peak(A)
    	B = self.mask_peak(B)
    
    	for i in range(13):
    		st = i*20
    		et = st + 20
    
    		on  = copy.deepcopy(A[:,i])
    		off = copy.deepcopy(A[:,i])
    		on[:st]   = ma.masked
    		on[et:]   = ma.masked
    		off[st:et]= ma.masked
    		flux.append(ma.mean(on) - ma.mean(off))
    
    		on  = copy.deepcopy(B[:,i])
    		off = copy.deepcopy(B[:,i])
    		on[:st]   = ma.masked
    		on[et:]   = ma.masked
    		off[st:et]= ma.masked
    		flux.append(ma.mean(on) - ma.mean(off))
    
    	# noise cal on
    	A = data[...,3]/(data[...,3] - data[...,4])
    	B = data[...,5]/(data[...,5] - data[...,6])
    
    	A = self.mask_peak(A)
    	B = self.mask_peak(B)
    
    	for i in range(13):
    		st = i*20
    		et = st + 20
    
    		on  = copy.deepcopy(A[:,i])
    		off = copy.deepcopy(A[:,i])
    		on[:st]   = ma.masked
    		on[et:]   = ma.masked
    		off[st:et]= ma.masked
    		flux.append(ma.mean(on) - ma.mean(off))
    
    		on  = copy.deepcopy(B[:,i])
    		off = copy.deepcopy(B[:,i])
    		on[:st]   = ma.masked
    		on[et:]   = ma.masked
    		off[st:et]= ma.masked
    		flux.append(ma.mean(on) - ma.mean(off))
    
    	flux = ma.array(flux).reshape(2,13,2)
    
    	flux = ma.mean(flux, axis=0)
    
    	return flux.T


    def mask_moving_time(self, data, moving_beam, length=6):
        print moving_beam, length
    	st = moving_beam*5*4
    	et = st + length*4
    	data[st:et,...] = ma.masked
    	data = data[~data.mask].reshape(-1,13,7)
    	return data
    
    def mask_peak(self, data, sig=3):
    	shape = data.shape
        print shape
    	data = data.reshape((-1,20)+shape[1:])
    	data_mean = ma.mean(data, axis=1)
    	data_std = ma.std(data, axis=1)
    	data[data>data_mean[:,None,...]+sig*data_std[:,None,...]] = ma.masked
    	data[data<data_mean[:,None,...]-sig*data_std[:,None,...]] = ma.masked
    
    	data = data.reshape(shape)
    	return data

