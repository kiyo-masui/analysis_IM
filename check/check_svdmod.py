#! /usr/bin/env python
import cPickle
import numpy as np
import matplotlib.pyplot as plt


class SVD(object):
    def __init__(self, path=None):
        self.svdlist = []
        self.labellist = []
        if not path==None:
            if isinstance(path, str):
                self.open(path)
                self.svdlist.append([self.svd_value, self.svd_lvect, self.svd_rvect])
            elif isinstance(path, list):
                for onepath in path:
                    self.open(onepath)
                    self.svdlist.append(
                        [self.svd_value, self.svd_lvect, self.svd_rvect])
                    self.labellist.append(onepath.split('.')[0].split('/')[-1])

    def open(self, path):
        self.svdmodpkl = cPickle.load(open(path))
        self.svd_value = self.svdmodpkl[0]
        #print self.svd_value
        self.svd_lvect = np.array(self.svdmodpkl[1])
        self.svd_rvect = np.array(self.svdmodpkl[2])
        print self.svd_lvect.shape, self.svd_lvect.dtype
        print self.svd_rvect.shape, self.svd_lvect.dtype
        self.tab_l = path.split('.')[0].split('_')[-3]
        self.tab_r = path.split('.')[0].split('_')[-1]

    def cal_eig(self, path):
        corr, count = cPickle.load(open(path))
        print corr.dtype
        corr_eigvalue, corr_eigvector = np.linalg.eig(corr[:,:,0])
        print corr_eigvector[:,0:3]
        print corr_eigvalue

    def plotvalue(self, savename, mod_start=0, mod_stop=50, bias=0, 
                  svd_start=0, svd_stop=-1, save=True, new=True):
        if new:
            plt.figure(figsize=(7,6))
        x = range(mod_start+bias, mod_stop+bias)
        #x = range(self.svd_lvect.shape[1])
        for (i, svd) in enumerate(self.svdlist[svd_start:svd_stop]):
            print i
            plt.plot(x, svd[0][mod_start:mod_stop], 'o',
                    label='eigenvalue %s'%self.labellist[i+svd_start],)
        plt.semilogy()
        plt.ylim(ymax=30, ymin=1.e-6)
        plt.xlim(xmin=-2)
        plt.xlabel('Mode Number')
        plt.ylabel('Eigenvalue')
        plt.legend(frameon=False, ncol=2)

        if save:
            plt.savefig('./png/'+savename+'_eigenvalue.png', format='png')

    def plot_real(self, savename, mod_start=0, mod_stop=6):
        svd_lvect_real = self.svd_lvect.real
        svd_rvect_real = self.svd_rvect.real

        self.plot(savename + '_real', mod_start, mod_stop, 
                  svd_lvect_real, svd_rvect_real)

    def plot_imag(self, savename, mod_start=0, mod_stop=6):
        svd_lvect_imag = self.svd_lvect.imag
        svd_rvect_imag = self.svd_rvect.imag

        self.plot(savename + '_imag', mod_start, mod_stop,
                  svd_lvect_imag, svd_rvect_imag)

    def plot(self, savename, mod_start=0, mod_stop=6,\
             svd_lvect_raw=None, svd_rvect_raw=None, combined=1):
        if svd_lvect_raw is None:
            svd_lvect_raw = self.svd_lvect
        if svd_rvect_raw is None:
            svd_rvect_raw = self.svd_rvect

        nfreq = 256 
        cutlist = [6, 7, 8, 15, 16, 18, 19, 20, 21, 22, 37, 80, 103, 104, 105, 106, \
                   107, 108, 130, 131, 132, 133, 134, 171, 175, 177, 179, 182, 183, \
                   187, 189, 192, 193, 194, 195, 196, 197, 198, 201, 204, 208, 209, \
                   212, 213, 218, 219, 229, 233, 237, 244, 254, 255]
        freq_list = tuple([ind for ind in range(nfreq) if ind not in cutlist])
        for i in range(1, combined):
            freq_list += tuple([ind + i*nfreq\
                                for ind in range(nfreq) if ind not in cutlist])
        nfreq = combined*nfreq

        xl = np.arange(nfreq)
        xr = np.arange(nfreq)
        svd_lvect = np.zeros(nfreq)
        svd_rvect = np.zeros(nfreq)

        f, ax = plt.subplots( mod_stop-mod_start, 1, sharex=True, figsize=(7,9))
        plt.subplots_adjust(hspace=0)
        for i in range(mod_start, mod_stop):
            np.put(svd_lvect, freq_list, svd_lvect_raw[i])
            svd_lvect[cutlist] = np.ma.masked
            ax[i-mod_start].plot(xl, svd_lvect, 
                                 label='%s svd mode %d'%(self.tab_l, i))

            np.put(svd_rvect, freq_list, svd_rvect_raw[i])
            svd_rvect[cutlist] = np.ma.masked
            ax[i-mod_start].plot(xr, svd_rvect, 
                                 label='%s svd mode %d'%(self.tab_r, i))

            ax[i-mod_start].set_ylim(ymin=-0.3, ymax=0.3)
            ax[i-mod_start].set_xlim(xmin=0, xmax=nfreq)
            ax[i-mod_start].legend(frameon=False, ncol=2)


        #plt.figure(figsize=(7,6))
        #f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7,6))
        #plt.subplots_adjust(hspace=0)

        #for i in range(mod_start, mod_stop):
        #    np.put(svd_lvect, freq_list, svd_lvect_raw[i])
        #    svd_lvect[cutlist] = np.ma.masked
        #    ax1.plot(xl, svd_lvect, label='%s svd mode %d'%(self.tab_l, i))
        #ax1.set_ylim(ymin=-0.3, ymax=0.3)
        #ax1.legend(frameon=False, ncol=2)

        #for i in range(mod_start, mod_stop):
        #    np.put(svd_rvect, freq_list, svd_rvect_raw[i])
        #    svd_rvect[cutlist] = np.ma.masked
        #    ax2.plot(xr, svd_rvect, label='%s svd mode %d'%(self.tab_r, i))
        #ax2.set_ylim(ymin=-0.3, ymax=0.3)
        #ax2.set_xlim(xmin=xr[0], xmax=xr[-1])
        #ax2.legend(frameon=False, ncol=2)


        plt.savefig('./png/'+savename+'_%d-%dmod.png'%(mod_start, mod_stop),
                    format='png')
        


if __name__=="__main__":
    #path = '/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/IxI5svd/'
    #path = '/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/IQUmap_clean_themselves/'
    #path = '/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/IQUmap_complex_clean_themselves/'
    #path = '/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_1gwj/IxI5svd/'
    #path = '/mnt/raid-project/gmrt/ycli/foreground_cleand/Parkes_sept12-14_west_legendre_modes_0gwj/mapmode_map/'
    #path = '/mnt/raid-project/gmrt/ycli/foreground_cleand/GBT_15hr_map_oldcal_legendre_modes_20gwj/mapmode_map/'
    #path = '/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj_conv/IQUmap_clean_themselves/'
    #path = '/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/IxI5svd/'
    path = '/mnt/raid-project/gmrt/ycli/foreground_cleand/IQUV_legendre_modes_0gwj/Emap_clean_themselves/'

    svd = SVD()

    #svdfile = 'SVD_pair_I_with_P.pkl'
    #svd.open(path+svdfile)
    #svd.plot_real(svdfile.split('.')[0], 0, 3)
    #svd.plot_imag(svdfile.split('.')[0], 0, 3)

    #svdfile = 'SVD_pair_I_with_I.pkl'
    #svd.open(path+svdfile)
    #svd.plot(svdfile.split('.')[0], 0, 4)

    #svdfile = 'SVD_pair_A_with_B.pkl'
    #corfile = 'foreground_corr_pair_A_with_B.pkl'
    #svd.cal_eig(path+corfile)
    #svd.open(path+svdfile)
    #svd.plot(svdfile.split('.')[0], 0, 4)
    #svd.plot(svdfile.split('.')[0], 4, 8)

    svdfile = 'SVD_pair_E_with_E.pkl'
    svd.open(path+svdfile)
    svd.plot(svdfile.split('.')[0], 0, 4, combined=4)
    #svd.plot(svdfile.split('.')[0], 4, 8)
    #svd.plotvalue(svdfile.split('.')[0])

    #svdfile = 'SVD_pair_I_with_U.pkl'
    #svd.open(path+svdfile)
    #svd.plot(svdfile.split('.')[0], 0, 4)
    #svd.plot(svdfile.split('.')[0], 4, 8)
    #svd.plotvalue(svdfile.split('.')[0])

    #path = '/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/'
    #svd = SVD([ path+'IQUmap_clean_themselves/SVD_pair_I_with_Q.pkl', 
    #            path+'IQUmap_clean_themselves/SVD_pair_I_with_U.pkl',
    #            path+'IQUmap_complex_clean_themselves/SVD_pair_I_with_P.pkl',
    #            path+'IxI5svd/SVD_pair_I_with_I.pkl',])
    #            #path+'IxI5svd_no_factorizable_noise/SVD_pair_I_with_I.pkl'])
    #svd.plotvalue('SVD_pair_IQ_IU_complex_eigenvalue', bias=3, svd_stop=3, save=False)
    #svd.plotvalue('SVD_pair_IQ_IU_complex_eigenvalue', svd_start=3, svd_stop=4, new=False)

    #svd = SVD([path+'SVD_pair_I_with_I.pkl',])
    #svd.plotvalue('parkes_SVD_pair_I_with_I', svd_stop=1)

