#! /usr/bin/env python

import pyfits
import numpy as np
import scipy as sp
from math import *

# Modified module
import data_block
import fitsGBT

from kiyopy import parse_ini

params_init = {
    'parkesroot' : '/mnt/raid-project/gmrt/raid-pen/pen/Parkes/2dF/DATA/p641/sdfits/rawdata/sept12/west/',
    'parkesfile' : [
                    '2008-09-12_1530_west1_1315_P641.sdfits',
                    '2008-09-12_1534_west2_1315_P641.sdfits',
                   ],
    'outputname' : 'test',
    'outputroot' : './',
    }

prefix = 'rfparkes_'

class ReformatParkes(object):
    
    def __init__(self, parameter_file_or_dict=None, feedback=2):
        
        self.params = parse_ini.parse(parameter_file_or_dict,
            params_init, prefix=prefix, feedback=feedback)

        self.feedback=feedback

        self.parkes_fieldlabel = []
        self.parkes_data = None

    def execute(self, nprocesses=1):
        self.readparkesdata()
        self.reformat()
        self.save()
        #self.check_data()

    def reformat(self):

        self.datablock = data_block.DataBlock()
        
        data_array, scanlist, band_cf, band_wd,\
        chan_cf,    chan_wd,  object,  exposure,\
        date_obs,   elevation,azimuth, beam, pol, cal\
            = self.readoutinfo(self.parkes_data)

        if self.feedback>2:
            print 'Dim:'
            print 'data:     ', data_array.shape
            print 'scan:     ', scanlist.shape
            print 'band_cf:  ', band_cf.shape
            print 'band_wd:  ', band_wd.shape
            print 'chan_cf:  ', band_cf.shape
            print 'chan_wd:  ', band_wd.shape
            print 'object:   ', object.shape
            print 'exposure: ', exposure.shape
            print 'data_obs: ', date_obs.shape
            print 'elevation:', elevation.shape
            print 'azimuth:  ', azimuth.shape
            print 'beam:     ', beam.shape
            print 'pol:      ', pol.shape
            print 'cal:      ', cal.shape

        self.datablock.set_data(data_array)
        self.datablock.set_field('SCAN',    tuple(scanlist), ('time',), '1I')
        self.datablock.set_field('CRVAL1',  tuple(band_cf),  (), '1D')
        self.datablock.set_field('BANDWID', tuple(band_wd),  (), '1D')
        self.datablock.set_field('CRPIX1',  tuple(chan_cf),  (), '1E')
        self.datablock.set_field('CDELT1',  tuple(chan_wd),  (), '1D')
        self.datablock.set_field('OBJECT',  tuple(object),   ('time',), '16A')
        self.datablock.set_field('EXPOSURE',tuple(exposure), ('time',), '1D')
        self.datablock.set_field('DATE-OBS',tuple(date_obs), ('time',), '22A')
        self.datablock.set_field('CRVAL3',  tuple(elevation),('time','beam',),'1D')
        self.datablock.set_field('CRVAL2',  tuple(azimuth),  ('time','beam',),'1D')
        self.datablock.set_field('BEAM',    tuple(beam),     ('beam',), '1I')
        self.datablock.set_field('CRVAL4',  tuple(pol),      ('pol', ), '1E')
        self.datablock.set_field('CAL',     tuple(cal),      ('cal', ), '1A')

        self.datablock.verify()

    def save(self):
        par = self.params
        print self.datablock.dims
        fitsfile = fitsGBT.Writer(self.datablock)
        fitsfile.write( par['outputroot'] + par['outputname'] + '.fits' )



    def readoutinfo(self, parkes_data):
        data = []
        band_cf = []
        band_wd = []
        chan_cf = []
        chan_wd = []
        object = []
        exposure = []
        date = []
        time = []
        elevation = []
        azimuth = []
        pol = [1,2]
        cal = ['F',]
        for i in range(len(parkes_data)):
            data.append(parkes_data[i].field('DATA'))
            band_cf.append(parkes_data[i].field('CRVAL1'))
            band_wd.append(parkes_data[i].field('BANDWID'))
            chan_cf.append(parkes_data[i].field('CRPIX1'))
            chan_wd.append(parkes_data[i].field('CDELT1'))
            object.append(parkes_data[i].field('OBJECT'))
            exposure.append(parkes_data[i].field('EXPOSURE'))
            date.append(parkes_data[i].field('DATE-OBS'))
            time.append(parkes_data[i].field('TIME'))
            elevation.append(parkes_data[i].field('ELEVATIO'))
            azimuth.append(parkes_data[i].field('AZIMUTH'))
        data = np.array(data)
        band_cf = np.array(band_cf)
        band_wd = np.array(band_wd)
        chan_cf = np.array(chan_cf)
        chan_wd = np.array(chan_wd)
        object = np.array(object)
        exposure = np.array(exposure)
        date = np.array(date)
        time = np.array(time)
        elevation = np.array(elevation)
        azimuth = np.array(azimuth)

        if data.ndim != 6:
            print 'error: data dim : '
            exit()
        #data = data[:,beam-1::13,:,:,:,:]
        # Note:  The raw Parkes data have 5 dims. We
        #        combined several files, so 6 dims:
        #        they are:
        #        shape[0] = file
        #        shape[1] = time*beam
        #        shape[2] = 1 ? not sure about it
        #        shape[3] = cal
        #        shape[4] = pol
        #        shape[5] = freq
        shape = data.shape
        if shape[1]%13 != 0:
            print 'error: data beam number error'
            exit()
        # Note: new data[file*time, beam, pol, cal, frequency]
        data = data.reshape(shape[0]*shape[1]/13,
                            13,
                            shape[4]*shape[2],
                            shape[3],
                            shape[5])
        #data_new = np.zeros([shape[0]*shape[1],4,shape[3],shape[5]])
        #data_new[:,::3,:,:] = data

        scannum = shape[0]
        scanlist = np.ones([shape[1]/13*shape[2]*shape[3], shape[0]])
        scanlist = scanlist * np.array(range(scannum))
        scanlist = scanlist.T.flatten()

        band_cf = np.array([band_cf[:,::13].flatten()[0]]) 
        band_wd = np.array([band_wd[:,::13].flatten()[0]]) 
        chan_cf = np.array([chan_cf[:,::13].flatten()[0]]) 
        chan_wd = np.array([chan_wd[:,::13].flatten()[0]]) 
        object = object[:,::13].flatten()
        exposure = exposure[:,::13].flatten()
        date_obs = self.getobstime(date[:,::13].flatten(), 
                                   time[:,::13].flatten())
        elevation = elevation[:,::13].flatten()
        azimuth = azimuth[:,::13].flatten()
        elevation, azimuth = self.elaz_multibeam(elevation, azimuth)
        pol = np.array(pol)
        cal = np.array(cal)
        beam= np.array(range(13))
        return data, scanlist, band_cf, band_wd, chan_cf, chan_wd, object, exposure,\
               date_obs, elevation, azimuth, beam, pol, cal

    def elaz_multibeam(self, elevation, azimuth):
        angle = 15. # Scan angle
        r_in = 29.1 # arc minutes
        r_out= 50.8 # arc minutes
        beamoffsetxy = np.array((np.arange(13.), np.arange(13.))).T
        # for beam 2 - 7 (index: 1 - 6)
        beamoffsetxy[1:7] = ((3 - beamoffsetxy[1:7])*60.+30.-angle)*pi/180.
        beamoffsetxy[1:7,0] = np.cos(beamoffsetxy[1:7,0])*r_in
        beamoffsetxy[1:7,1] = np.sin(beamoffsetxy[1:7,1])*r_in
        # for beam 8 - 13 (index: 7 - 12)
        beamoffsetxy[7:13] = ((10 - beamoffsetxy[7:13])*60.-angle)*pi/180.
        beamoffsetxy[7:13,0] = np.cos(beamoffsetxy[7:13,0])*r_out
        beamoffsetxy[7:13,1] = np.sin(beamoffsetxy[7:13,1])*r_out

        elevation = elevation.repeat(13).reshape(-1,13)
        elevation[:] += beamoffsetxy[:,1]/60.

        delta_azimuth = beamoffsetxy[:,0]/np.cos(elevation[:]*pi/180.)
        azimuth = azimuth.repeat(13).reshape(-1,13)
        azimuth[:] -= delta_azimuth[:]/60

        return elevation, azimuth

        
    def getobstime(self, date, time):
        if date.shape != time.shape:
            print 'error: time date not match'
            exit()
    
        time = time.astype(float)
        h = (time/60./60.).astype(int)
        m = ((time - h*60*60)/60.).astype(int)
        s = time - h*60*60 - m*60

        date_obs = []
        for i in range(date.shape[0]):
            date_obs.append(date[i] + 'T%d:%d:%f'%(h[i], m[i], s[i]))
        return np.array(date_obs)


    def readparkesdata(self):
        
        par = self.params
        parkes_root = par[ 'parkesroot' ]
        parkes_file = par[ 'parkesfile' ]

        first_parkesfile = True
        self.parkes_data = []

        for (i, file) in enumerate(parkes_file):
            hdulist = pyfits.open( parkes_root + file )

            # get all the label name
            for j in range(hdulist[1].header['TFIELDS']):
                if first_parkesfile == True:
                    self.parkes_fieldlabel.append(hdulist[1].header['TTYPE%d'%(j+1)])
                else:
                    if self.parkes_fieldlabel[j] != hdulist[1].header['TTYPE%d'%(j+1)]:
                        print 'fits label not match'
                        exit()

            # get the data
            self.parkes_data.append(hdulist[1].data)

            if self.feedback>3:
                print '-'*50
                print 'sdfits label of : [%s]' %file
                print hdulist[1].header
                print
                for j in range(hdulist[1].header['TFIELDS']):
                    print hdulist[1].header['TTYPE%d'%(j+1)],
                    print sp.unique(self.parkes_data[i].field(
                                    self.parkes_fieldlabel[j])).shape,
                    print self.parkes_data[i].field(self.parkes_fieldlabel[j]).shape
                    print self.parkes_data[i].field(self.parkes_fieldlabel[j])
                    print

            first_parkesfile = False

    def check_data(self):
        self.datablock.calc_pointing()
        self.datablock.calc_freq()
        print self.datablock.ra.tolist()
        print self.datablock.freq/1.e6


if __name__ == '__main__':
    import sys
    if len(sys.argv)==2 :
        ReformatParkes(str(sys.argv[1]), feedback=3).execute()
    elif len(sys.argv)>2 :
        print 'Maximun one argument, a parameter file name.'
    else:
        ReformatParkes(feedback=3).execute()
