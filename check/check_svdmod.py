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
        print self.svd_value
        self.svd_lvect = np.array(self.svdmodpkl[1])
        self.svd_rvect = np.array(self.svdmodpkl[2])
        print self.svd_lvect.shape
        print self.svd_rvect.shape
        self.tab_l = path.split('.')[0].split('_')[-3]
        self.tab_r = path.split('.')[0].split('_')[-1]

    def plotvalue(self, savename, mod_start=0, mod_stop=50):
        plt.figure(figsize=(8,4))
        x = range(mod_start, mod_stop)
        #x = range(self.svd_lvect.shape[1])
        for (i, svd) in enumerate(self.svdlist):
            plt.plot(x, svd[0][mod_start:mod_stop], 
                    label='eigenvalue %s'%self.labellist[i],)
        plt.semilogy()
        plt.ylim(ymax=1, ymin=1.e-6)
        plt.legend(frameon=False, ncol=2)

        plt.savefig('./png/'+savename+'_eigenvalue.png', format='png')

    def plot(self, savename, mod_start=0, mod_stop=6):
        xl = range(self.svd_lvect.shape[1])
        xr = range(self.svd_rvect.shape[1])
        plt.figure(figsize=(7,8))

        plt.subplot(211)
        for i in range(mod_start, mod_stop):
            plt.plot(xl, self.svd_lvect[i], label='%s svd mode %d'%(self.tab_l, i))
        plt.ylim(ymin=-0.3, ymax=0.3)
        plt.legend(frameon=False, ncol=2)

        plt.subplot(212)
        for i in range(mod_start, mod_stop):
            plt.plot(xr, self.svd_rvect[i], label='%s svd mode %d'%(self.tab_r, i))
        plt.ylim(ymin=-0.3, ymax=0.3)
        plt.legend(frameon=False, ncol=2)


        plt.savefig('./png/'+savename+'_%d-%dmod.png'%(mod_start, mod_stop),
                    format='png')
        


if __name__=="__main__":
    path = '/mnt/raid-project/gmrt/ycli/foreground_cleand/IQV_legendre_modes_0gwj/IxI5svd/'
    #path = '/mnt/raid-project/gmrt/ycli/foreground_cleand/IQV_legendre_modes_0gwj/IQUmap_clean_themselves/'

    svd = SVD()

    svdfile = 'SVD_pair_I_with_I.pkl'
    svd.open(path+svdfile)
    svd.plot(svdfile.split('.')[0], 3, 6)
    #svd.plotvalue(svdfile.split('.')[0])

    #svdfile = 'SVD_pair_I_with_Q.pkl'
    #svd.open(path+svdfile)
    #svd.plot(svdfile.split('.')[0], 0, 4)
    #svd.plot(svdfile.split('.')[0], 4, 8)
    #svd.plotvalue(svdfile.split('.')[0])

    #svdfile = 'SVD_pair_I_with_U.pkl'
    #svd.open(path+svdfile)
    #svd.plot(svdfile.split('.')[0], 0, 4)
    #svd.plot(svdfile.split('.')[0], 4, 8)
    #svd.plotvalue(svdfile.split('.')[0])

    #svd = SVD([ path+'SVD_pair_I_with_Q.pkl', path+'SVD_pair_I_with_U.pkl'])
    #svd.plotvalue('SVD_pair_IQ_IU_eigenvalue')
