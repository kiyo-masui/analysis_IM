#! /usr/bin/env python
import cPickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from core import algebra


class SVD(object):
    def __init__(self, path=None):
        self.svdlist = []
        self.labellist = []
        if not path==None:
            if isinstance(path, str):
                self.open(path)
                self.svdlist.append([self.svd_value, self.svd_lvect, self.svd_rvect])
            elif isinstance(path, list):
                path = np.array(path)
                self.biaslist = path[:,0]
                self.lprelist = path[:,2]
                for onepath in path[:,1]:
                    self.open(onepath)
                    self.svdlist.append(
                        [self.svd_value, self.svd_lvect, self.svd_rvect])
                    self.labellist.append(onepath.split('.')[0].split('/')[-1])

    def open(self, path):
        self.path = path.replace(path.split('/')[-1], '')
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

    def plotvalue(self, savename, mod_start=0, mod_stop=50):
        plt.figure(figsize=(7,7))
        for i in range(len(self.svdlist)):
            bias = int(self.biaslist[i])
            x = range( mod_start+bias, mod_stop+bias )
            plt.plot(x, self.svdlist[i][0][mod_start:mod_stop], '-', marker='.',
                    label=self.lprelist[i]+'eigenvalue %s'%self.labellist[i],)
        plt.semilogy()
        plt.ylim(ymax=20, ymin=1.e-7)
        plt.xlim(xmin=-2)
        plt.xlabel('Mode Number')
        plt.ylabel('Eigenvalue')
        #plt.legend(frameon=False, ncol=2)
        plt.legend(frameon=False)

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

    def plot_svdmap(self, savename, mod_number, mod_space=5, right=False):
        if mod_number <1:
            print "Error: mod_number must >= 1"
            exit()
        mod_n = int((mod_number-1)%mod_space)
        mod_f = int(mod_number-1-mod_n+mod_space)
        mod_map = self.path
        if right:
            mod_map+= 'sec_%s_modes_clean_map_I_with_%s_%dmodes.npy'\
                     %(self.tab_r, self.tab_l, mod_f)
            mod = self.svd_rvect[mod_number-1]
        else:
            mod_map+= 'sec_%s_modes_clean_map_I_with_%s_%dmodes.npy'\
                     %(self.tab_l, self.tab_r, mod_f)
            mod = self.svd_lvect[mod_number-1]
        print "Mode Amp from: [%s]" % mod_map
        print 

        mod_amp = algebra.make_vect(algebra.load(mod_map))
        ra = mod_amp.get_axis('ra')
        dec= mod_amp.get_axis('dec')
        mod_amp = mod_amp[mod_n]

        map = mod[:,None,None] * mod_amp[None,:,:]

        map = np.ma.array(map)
        map[map==0] = np.ma.masked
        mod_amp = np.ma.array(mod_amp)
        mod_amp[mod_amp==0] = np.ma.masked

        halfindx = map.shape[0]//2
        vmax = 0.0002
        vmin = -0.0002

        dra = ra[1]-ra[0]
        ddec = dec[1]-dec[0]
        ra = np.resize(ra-0.5*dra,  ra.shape[0]+1)
        dec= np.resize(dec-0.5*ddec,dec.shape[0]+1)
        ra[-1] = ra[-2]+dra
        dec[-1]= dec[-2]+ddec
        f = plt.figure(figsize=(15, 21))
        ax = ImageGrid(f, 111,
                       nrows_ncols = (3, 1),
                       direction = "column",
                       axes_pad = 0.4,
                       add_all = True,
                       label_mode = "L",
                       share_all = True,
                       cbar_location = "right",
                       cbar_mode = "each",
                       cbar_size = "5%",
                       cbar_pad = 0.05,
                       )
        #im0 = ax[0].pcolormesh(ra, dec, np.ma.mean(map[:halfindx], 0))
        im0 = ax[0].pcolormesh(ra, dec, np.ma.mean(map[:halfindx], 0).swapaxes(0,1))
        im0.set_clim(vmin, vmax)
        ax[0].set_xlim(ra.min(), ra.max())
        ax[0].set_ylim(dec.min(), dec.max())
        ax[0].set_xlabel('RA [deg]')
        ax[0].set_ylabel('DEC[deg]')
        ax[0].set_title('average over [0:%d]'%halfindx)
        ax[0].cax.colorbar(im0)

        im1 = ax[1].pcolormesh(ra, dec, np.ma.mean(map[halfindx:], 0).swapaxes(0,1))
        im1.set_clim(vmin, vmax)
        ax[1].set_xlim(ra.min(), ra.max())
        ax[1].set_ylim(dec.min(), dec.max())
        ax[1].set_xlabel('RA [deg]')
        ax[1].set_ylabel('DEC[deg]')
        ax[1].set_title('average over [%d:%d]'%(halfindx, map.shape[0]))
        ax[1].cax.colorbar(im1)

        im2 = ax[2].pcolormesh(ra, dec, mod_amp.swapaxes(0,1))
        #im2.set_clim(vmin, vmax)
        ax[2].set_xlim(ra.min(), ra.max())
        ax[2].set_ylim(dec.min(), dec.max())
        ax[2].set_xlabel('RA [deg]')
        ax[2].set_ylabel('DEC[deg]')
        ax[2].set_title('Amp of mode %d (start with 0)'%(mod_number-1))
        ax[2].cax.colorbar(im2)

        plt.tick_params(length=6, width=1.)
        plt.tick_params(which='minor', length=3, width=1.)

        plt.savefig('./png/'+savename+'_%dmodmap.png'%(mod_number), format='png')

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
        svd_lvect = np.ones(nfreq)*1.e20
        svd_rvect = np.ones(nfreq)*1.e20

        f, ax = plt.subplots( mod_stop-mod_start, 1, sharex=True, figsize=(7,9))
        plt.subplots_adjust(hspace=0)
        for i in range(mod_start, mod_stop):
            np.put(svd_lvect, freq_list, svd_lvect_raw[i])
            svd_lvect = np.ma.masked_values(svd_lvect, 1.e20)
            ax[i-mod_start].plot(xl, svd_lvect, 
                                 label='%s svd mode %d'%(self.tab_l, i))

            np.put(svd_rvect, freq_list, svd_rvect_raw[i])
            svd_rvect = np.ma.masked_values(svd_rvect, 1.e20)
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

    #from mkpower import functions

    #map = algebra.make_vect(algebra.load('/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_AQUV_extend_legendre_modes_0gwj_conv/Emap_clean_themselves/sec_A_cleaned_clean_map_I_with_B_20modes.npy'))

    #cutlist = [ 171, 175, 177, 179, 182, 183, 187, 189, 192, 193, 194, 195, 196, 197, 198, 201, 204, 208, 209, 212, 213, 218, 219, 229, 233, 237, 244, 254, 255]
    #for i in range(len(cutlist)):
    #    cutlist[i] -= 143
    #print cutlist

    #print map.shape
    #print map.get_axis('freq')[0]
    #print map.get_axis('freq')[-1]
    ##map = functions.getmap_halfz(map, 'upper')
    #map = functions.getmap_halfz(map, 'lower')
    #print map.shape
    #print map.get_axis('freq')[0]
    #print map.get_axis('freq')[-1]
    #exit()



    do_plot_svd = False
    do_calc_eig = False
    do_plot_val = True
    do_plot_map = False

    #fg_root = '/mnt/raid-project/gmrt/ycli/foreground_cleand/'
    fg_root = '/mnt/scratch-gl/ycli/cln_result/'
    combined = 1

    '''>>>  parameters for plot svd eigenvector  <<<'''

    r'''I_AxI_{BCD}QUV'''
    #fg_file = '1hr_AQUV_extend_legendre_modes_0gwj_2conv/Emap_clean_themselves/'
    #fg_file = '15hr_AQU_extend_legendre_modes_0gwj_14conv_new/Emap_clean_themselves/'
    fg_file = '15hr_ABCD_legendre_modes_0gwj_14conv_new/Emap_clean_themselves/'
    svdfile = 'SVD_pair_A_with_B.pkl'
    #svdfile = 'SVD_pair_A_with_C.pkl'
    #svdfile = 'SVD_pair_A_with_D.pkl'
    combined= 3

    r'''I_BxI_{ACD}QUV'''
    #fg_file = '1hr_BQUV_extend_legendre_modes_0gwj_conv/Emap_clean_themselves/'
    #svdfile = 'SVD_pair_B_with_A.pkl'
    #svdfile = 'SVD_pair_B_with_C.pkl'
    #svdfile = 'SVD_pair_B_with_D.pkl'
    #combined= 4

    r'''I_CxI_{ABD}QUV'''
    #fg_file = '1hr_CQUV_extend_legendre_modes_0gwj_conv/Emap_clean_themselves/'
    #svdfile = 'SVD_pair_C_with_A.pkl'
    #svdfile = 'SVD_pair_C_with_B.pkl'
    #svdfile = 'SVD_pair_C_with_D.pkl'
    #combined= 4

    r'''I_DxI_{ABC}QUV'''
    #fg_file = '1hr_DQUV_extend_legendre_modes_0gwj_conv/Emap_clean_themselves/'
    #svdfile = 'SVD_pair_D_with_A.pkl'
    #svdfile = 'SVD_pair_D_with_B.pkl'
    #svdfile = 'SVD_pair_D_with_C.pkl'
    #combined= 4

    '''>>>  parameters for plot svd eigenvalue  <<<'''
    #svd_list =[
    #    [3, fg_root+'IQUmap_clean_themselves/SVD_pair_I_with_Q.pkl', ''],
    #    [3, fg_root+'IQUmap_clean_themselves/SVD_pair_I_with_U.pkl', ''],
    #    [3, fg_root+'IQUmap_complex_clean_themselves/SVD_pair_I_with_P.pkl', ''],
    #    [0, fg_root+'IxI5svd/SVD_pair_I_with_I.pkl', ''],]
    #svd_savename = 'SVD_pair_IQ_IU_complex_eigenvalue'

    svd_list = [
        #[0, fg_root+'15hr_AQU_extend_legendre_modes_0gwj_14conv_new/Emap_clean_themselves/SVD_pair_A_with_B.pkl', '15hr 1.4conv new '],
        #[0, fg_root+'15hr_AQU_extend_legendre_modes_0gwj_14conv_new/Emap_clean_themselves/SVD_pair_A_with_C.pkl', '15hr 1.4conv new '],
        #[0, fg_root+'15hr_BQU_extend_legendre_modes_0gwj_14conv_new/Emap_clean_themselves/SVD_pair_B_with_A.pkl', '15hr 1.4conv new '],
        #[0, fg_root+'15hr_AQU_extend_legendre_modes_0gwj_14conv/Emap_clean_themselves/SVD_pair_A_with_B.pkl', '15hr 1.4conv '],
        [0, fg_root+'15hr_ABCD_legendre_modes_0gwj_14conv_new/Emap_clean_themselves/SVD_pair_A_with_B.pkl', '15hr 1.4conv new IxI clip weight '],
        [0, fg_root+'15hr_ABCD_legendre_modes_0gwj_14conv_new_noclip/Emap_clean_themselves/SVD_pair_A_with_B.pkl', '15hr 1.4conv new IxI'],
        #[0, fg_root+'1hr_IQUV_extend_legendre_modes_0gwj_2conv/Emap_clean_themselves/SVD_pair_I_with_E.pkl', '1hr conv '],
        #[0, fg_root+'1hr_AQUV_extend_legendre_modes_0gwj_conv/Emap_clean_themselves/SVD_pair_A_with_B.pkl', '1hr conv '],
        #[0, fg_root+'1hr_AQUV_extend_legendre_modes_0gwj_2conv/Emap_clean_themselves/SVD_pair_A_with_B.pkl', '1hr 2conv '],
        #[0, fg_root+'1hr_AQUV_extend_legendre_modes_0gwj_2conv/Emap_clean_themselves/SVD_pair_A_with_C.pkl', '1hr 2conv '],
        #[0, fg_root+'1hr_AQUV_extend_legendre_modes_0gwj_2conv/Emap_clean_themselves/SVD_pair_A_with_D.pkl', '1hr 2conv '],
        ]
    svd_savename = 'SVD_pair_15hr_AxIQU_eigenvalue'


    if do_plot_svd:
        path = fg_root + fg_file
        svd = SVD()
        svd.open(path+svdfile)
        svd.plot(svdfile.split('.')[0], 0, 4, combined=combined)
        svd.plot(svdfile.split('.')[0], 4, 8, combined=combined)
    if do_calc_eig:
        path = fg_root + fg_file
        svd = SVD()
        corfile = 'foreground_corr_pair_A_with_B.pkl'
        svd.cal_eig(path+corfile)
        svd.open(path+svdfile)
        svd.plot(svdfile.split('.')[0], 0, 4, combined=combined)
        svd.plot(svdfile.split('.')[0], 4, 8, combined=combined)
    if do_plot_val:
        svd = SVD(svd_list)
        svd.plotvalue(svd_savename)
    if do_plot_map:
        path = fg_root + fg_file
        svd = SVD()
        svd.open(path+svdfile)
        svd.plot_svdmap(svdfile.split('.')[0], 20)

