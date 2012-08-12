import pyfits
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

from parkes import fitsGBT

#rawdatapath = ('/mnt/raid-project/gmrt/raid-pen/pen/Parkes/2dF/DATA/p641/sdfits/rawdata/sept11/east/2008-09-11_1217_east1_1392_P641.sdfits',)
#rawdatapath = ('/mnt/raid-project/gmrt/raid-pen/pen/Parkes/2dF/DATA/p641/sdfits/rawdata/sept12/west/2008-09-12_1534_west2_1315_P641.sdfits',)
#rawdatapath = ('/mnt/raid-project/gmrt/raid-pen/pen/Parkes/2dF/DATA/p641/sdfits/rawdata/sept11/west/2008-09-11_1647_west1_1290_drift_P641.sdfits',)

#rawdatapath = ('/mnt/raid-project/gmrt/ycli/86_wigglez1hr_centre_ralongmap_19-28.fits',)
#rawdatapath = ('/mnt/raid-project/gmrt/ycli/55_wigglez15hrst_ralongmap_272-279.fits',)
#rawdatapath = ('/home/ycli/workspace/map_result/parkes/parkes_2008_09_12_west_P641.fits',)
#rawdatapath = ('/mnt/data-pen3/ycli/map_result/parkes/parkes_2008_09_12_west_P641.fits',)
#rawdatapath = ('/mnt/data-pen3/ycli/map_result/flagged/parkes_2008_09_12_west_P641.fits',)
#rawdatapath = ('/mnt/data-pen3/ycli/map_result/rebinned/parkes_2008_09_12_west_P641.fits',)
rawdatapath = ('/mnt/data-pen3/ycli/map_result/pol_selected/parkes_2008_09_12_west_P641.fits',)

class CheckFitsFile(object):

    def __init__(self, datapath):
        self.hdulist = pyfits.open(datapath)
        
        self.tbdata = self.hdulist[1].data
        
        self.hdulist.close()

        self.reader = fitsGBT.Reader(datapath)

    def printhead(self):
        print self.hdulist[1].header
        
        
    def printlabel(self):
        self.fieldlabel = []
        
        for i in range(self.hdulist[1].header['TFIELDS']):
            self.fieldlabel.append(self.hdulist[1].header['TTYPE%d'%(i+1)])
        
        #for i in range(len(tbdata)):
            #print tbdata[i][fieldlabel[3]], 

        for i in range(self.hdulist[1].header['TFIELDS']):
            print self.hdulist[1].header['TTYPE%d'%(i+1)]
            print sp.unique(self.tbdata.field(self.fieldlabel[i])).shape
            print self.tbdata.field(self.fieldlabel[i]).shape
            print self.tbdata.field(self.fieldlabel[i])
            print
        
        #print self.tbdata.field('AZIMUTH')[:13]
        #print self.tbdata.field(fieldlabel[6])[:200]
        #print self.tbdata.field(fieldlabel[6])[200:1000]

    def plotfreq(self):
        spectrum_xx = self.tbdata.field('DATA')[0::2,:][0::13,:]
        spectrum_yy = self.tbdata.field('DATA')[1::2,:][0::13,:]
        x = range(spectrum_xx.shape[1])
        print spectrum_xx.shape
        print spectrum_yy.shape
        plt.figure(figsize=(8,5))
        #for i in range(100):
        for i in range(spectrum_xx.shape[0]):
            plt.plot(x, spectrum_xx[i], c='0.6')

        plt.ylim(ymax=350,ymin=0)
        plt.savefig('./png/parkes_test_pol.png', format='png')

    def plotradec(self):
        scan_inds = self.reader.scan_set
        plt.figure(figsize=(33,4))
        ax = plt.gca()
        beamfwhp = 14.0/60.
        for thisscan in scan_inds:
            block = self.reader.read(thisscan)
            block.calc_pointing('W')
            for i in range(block.ra.shape[0]):
                for j in range(block.ra.shape[1]):
                    cir = plt.Circle((block.ra[i,j], 
                                      block.dec[i,j]),
                                      radius=beamfwhp/2., 
                                      fc='b', 
                                      ec="none",
                                      alpha=0.1)
                    ax.add_patch(cir)
            ax.autoscale_view()
        plt.xlabel('RA [deg]')
        plt.ylabel('DEC [deg]')
        plt.tick_params(length=6, width=1.)
        plt.tick_params(which='minor', length=3, width=1.)
        plt.savefig('./png/parkes_test_radec.png', format='png')

    def plotradec_one(self):
        scan_inds = self.reader.scan_set
        plt.figure(figsize=(6,6))
        beamfwhp = 14.0/60.
        ax = plt.gca()
        for thisscan in scan_inds[0:1]:
            block = self.reader.read(thisscan)
            block.calc_pointing()
            print block.ra.shape
            print block.dec.shape
            for i in range(13):
                cir = plt.Circle((block.ra[0,i],
                                  block.dec[0,i]),
                                  radius=beamfwhp/2., 
                                  fc="none", 
                                  edgecolor='k')
                ax.add_patch(cir)
            ax.autoscale_view()
        plt.xlabel('RA [deg]')
        plt.ylabel('DEC [deg]')
        plt.xlim(xmin=315.5, xmax=318.5)
        plt.ylim(ymin=-37.5, ymax=-34.5)
        plt.tick_params(length=6, width=1.)
        plt.tick_params(which='minor', length=3, width=1.)
        plt.savefig('./png/parkes_test_radec_one.png', format='png')

    def plotelaz(self):
        scan_inds = self.reader.scan_set
        plt.figure(figsize=(15,3))
        ax = plt.gca()
        beamfwhp = 14.0/60.
        for thisscan in scan_inds:
            block = self.reader.read(thisscan)
            for i in range(block.field['CRVAL2'].shape[0]):
                for j in range(block.field['CRVAL2'].shape[1]):
                    cir = plt.Circle((block.field['CRVAL2'][i,j], 
                                      block.field['CRVAL3'][i,j]),
                                      radius=beamfwhp/2., 
                                      fc='b', 
                                      ec="none",
                                      alpha=0.003)
                    ax.add_patch(cir)
            ax.autoscale_view()
            #plt.scatter(block.field['CRVAL2'][:,:], block.field['CRVAL3'][:,:],
            #            c='w', alpha=0.3)
        plt.xlabel('Azimuth [deg]')
        plt.ylabel('Elevation [deg]')
        plt.tick_params(length=6, width=1.)
        plt.tick_params(which='minor', length=3, width=1.)
        plt.savefig('./png/parkes_test_elaz.png', format='png')

    def plotelaz_one(self):
        scan_inds = self.reader.scan_set
        plt.figure(figsize=(6,6))
        ax = plt.gca()
        beamfwhp = 14.0/60.
        for thisscan in scan_inds[0:1]:
            block = self.reader.read(thisscan)
            print block.field['CRVAL2'].shape
            print block.field['CRVAL3'].shape
            #plt.scatter(block.field['CRVAL2'][0,:], block.field['CRVAL3'][0,:],
            #            c='w', alpha=0.3)
            for i in range(13):
                cir = plt.Circle((block.field['CRVAL2'][0,i],
                                  block.field['CRVAL3'][0,i]),
                                  radius=beamfwhp/2., 
                                  fc='none', 
                                  ec='k')
                ax.add_patch(cir)
            ax.autoscale_view()
        plt.xlabel('Azimuth [deg]')
        plt.ylabel('Elevation [deg]')
        plt.tick_params(length=6, width=1.)
        plt.tick_params(which='minor', length=3, width=1.)
        plt.savefig('./png/parkes_test_elaz_one.png', format='png')
        
if __name__=="__main__":
    checkfits = CheckFitsFile(rawdatapath[0])

    #checkfits.printhead()
    #checkfits.printlabel()
    #checkfits.plotfreq()
    checkfits.plotradec()
    checkfits.plotelaz()
    checkfits.plotradec_one()
    checkfits.plotelaz_one()

