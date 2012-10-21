#! /usr/bin/env python

import scipy as sp
import numpy as np
from numpy.fft import *
import scipy.linalg as linalg
import multiprocessing as mp

from core import algebra, hist
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
from scipy import integrate
from math import *
from sys import *
import matplotlib.pyplot as plt
import mkpower
from mpi4py import MPI


pi = 3.1415926
deg2rad = pi/180.

params_init = {
    'processes' : 1,
    'plot' : False,
    'resultf' : '',
    'input_root' : '../../../jkmap/map/',
    'output_root' : '../../../powerresult/',

    'imap_list' : [],
    'nmap_list' : [],
    'mmap_list' : [],

    'imap_pair' : [],
    'nmap_pair' : [],
    'mmap_pair' : [],

    #'hrlist' : (),
    #'ltlist' : (),
    #'hr' : ('15hr_40-41-43_','15hr_42_',),
    #'mid' : ('dirty_map_',),
    #'polarizations' : ('I',),
    #'last' : (),

    'cldir' : '',

    'boxshape' : (128,128,128),
    'boxunit' : 15., # in unit Mpc
    'discrete' : 3,
    #'Xrange' : (1400,),
    #'Yrange' : (-64*15,64*15),
    #'Zrange' : (0.,128*15),
    'kbinNum' : 20,
    'kmin' : -1,
    'kmax' : -1,

    'FKPweight' : False,
    'FKPpk' : 1.e-3,
    'jkerror' : False,
    'sme' : True,
}
prefix = 'pkc_'


class PowerSpectrumMaker(mkpower.PowerSpectrumMaker):
    """Calculate The Power Spectrum"""
    B = []
    B2= []
    Bk = []
    q = mp.JoinableQueue() # for 3d power
    qn = mp.JoinableQueue() # for 3d power
    q2 = mp.JoinableQueue() # for 3d power
    qn2= mp.JoinableQueue() # for 2d power

    def __init__(self, parameter_file_or_dict=None, feedback=1):
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, 
            params_init, prefix=prefix, feedback=feedback)

        self.feedback=feedback

    def mpiexecute(self, nprocesses=1):
        
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        params = self.params

        n_map = len(params['imap_list'])

        comm.barrier()

        if rank<n_map:
            self.params['imap_list'] = self.params['imap_list'][rank::size]
            self.params['nmap_list'] = self.params['nmap_list'][rank::size]
            if len(self.params['mmap_list'])!=0:
                self.params['mmap_list'] = self.params['mmap_list'][rank::size]
            PK, kn, PK2, kn2 = self.execute(rank=rank)
            if rank != 0:
                print "RANK %d: Send Data"%rank
                comm.ssend(PK,  dest=0, tag=11)
                comm.ssend(PK2, dest=0, tag=12)
                comm.ssend(kn,  dest=0, tag=13)
                comm.ssend(kn2, dest=0, tag=14)
                print "RANK %d: Send Data Succeed"%rank
            else:
                if size > n_map:
                    n_rank = n_map
                else:
                    n_rank = size
                for i in range(1,n_rank):
                    print "RANK %d: Receive Data from rank %d"%(rank, i)
                    PK = np.append(PK,  comm.recv(source=i, tag=11), axis=0)
                    PK2= np.append(PK2, comm.recv(source=i, tag=12), axis=0)
                    kn = np.append(kn,  comm.recv(source=i, tag=13), axis=0)
                    kn2= np.append(kn2, comm.recv(source=i, tag=14), axis=0)
                    print "RANK %d: Receive Data from rank %d Succeed"%(rank, i)
                for i in range(n_rank, size):
                    sig = comm.recv(source=i, tag=11)
                    if sig==1:
                        print 'RANK %d: rank %d is free'%(rank, i)
                    else:
                        print 'RANK %d: rank %d has problem'%(rank, i)

                self.save(PK, kn, PK2, kn2)
        else:
            if rank != 0:
                comm.ssend(1, dest=0, tag=11)

        comm.barrier()


    
    def execute(self, nprocesses=1, rank=-1):

        params = self.params
        n_processes = params['processes']

        if params['cldir']!='':
            # Read in the bias calibration data
            #self.B = sp.load(params['cldir']+'b_each_bias.npy')
            self.B  = sp.load(params['cldir']+'b_bias.npy')
            self.Bk = sp.load(params['cldir']+'k_bias.npy')
            self.B2 = sp.load(params['cldir']+'b2_bias.npy')
            #print self.B


        #### Process ####
        n_new = n_processes - 1
        n_map = len(params['imap_list'])
        if rank!=-1:
            print "RANK %d: "%rank,
        print "%d power need to be calculated"%n_map

        kbn = params['kbinNum']
        kmin = params['kmin']
        kmax = params['kmax']

        PK = np.zeros(shape=(n_map, kbn))
        kn = np.zeros(shape=(n_map, kbn))

        PK2= np.zeros(shape=(n_map, kbn, kbn))
        kn2= np.zeros(shape=(n_map, kbn, kbn))

        if n_new <= 0:
            raise ValueError('Process Should >= 2')
            for ii in range(n_map):
                self.process_map(ii, PK[ii])
        elif n_new > 32:
            raise ValueError('Process limit is 32')
        else:
            process_list = range(n_new)
            for ii in xrange(n_map + n_new):
                if ii >= n_new:
                    PK[ii-n_new] = self.q.get()
                    kn[ii-n_new] = self.qn.get()
                    PK2[ii-n_new]= self.q2.get()
                    kn2[ii-n_new] = self.qn2.get()
                    process_list[ii%n_new].join()
                    if process_list[ii%n_new].exitcode != 0:
                        raise RuntimeError("A thred faild with exit code"
                            + str(process_list[ii%n_new].exitcode))
                if ii < n_map:
                    process_list[ii%n_new] = mp.Process(
                        target=self.process_map, 
                        args=(ii, ii%n_new))
                    process_list[ii%n_new].start()

        if rank==-1:
            self.save(PK, kn, PK2, kn2)
        else:
            print "RANK %d: Power Calculation Finished"%rank
            return PK, kn, PK2, kn2

    def save(self, PK, kn, PK2, kn2):
        params = self.params
        resultf = params['resultf']
        FKPweight = params['FKPweight']
        kbn = params['kbinNum']
        kmin = params['kmin']
        kmax = params['kmax']
        if FKPweight:
            if params['jkerror']:
                sp.save(params['output_root']+resultf+'_p_each_jk_fkp', PK)
                sp.save(params['output_root']+resultf+'_n_each_jk_fkp', kn)
                sp.save(params['output_root']+resultf+'_p2_each_jk_fkp', PK2)
                sp.save(params['output_root']+resultf+'_n2_each_jk_fkp', kn2)
            else:
                sp.save(params['output_root']+resultf+'_p_each_fkp', PK)
                sp.save(params['output_root']+resultf+'_n_each_fkp', kn)
                sp.save(params['output_root']+resultf+'_p2_each_fkp', PK2)
                sp.save(params['output_root']+resultf+'_n2_each_fkp', kn2)
        else:
            if params['jkerror']:
                sp.save(params['output_root']+resultf+'_p_each_jk', PK)
                sp.save(params['output_root']+resultf+'_n_each_jk', kn)
                sp.save(params['output_root']+resultf+'_p2_each_jk', PK2)
                sp.save(params['output_root']+resultf+'_n2_each_jk', kn2)
            else:
                sp.save(params['output_root']+resultf+'_p_each', PK)
                sp.save(params['output_root']+resultf+'_n_each', kn)
                sp.save(params['output_root']+resultf+'_p2_each', PK2)
                sp.save(params['output_root']+resultf+'_n2_each', kn2)

        PKmean = PK.mean(axis=0)
        PKvar = PK.std(axis=0)
        knmean = kn.mean(axis=0)
        knvar = kn.std(axis=0)

        PK2mean = PK2.mean(axis=0)
        PK2var = PK2.std(axis=0)
        kn2mean = kn2.mean(axis=0)
        kn2var = kn2.std(axis=0)

        if params['jkerror']:
            print '\tJackKnife Error:'
            PKvar = PKvar*sqrt((PK.shape[0]-1))
            PK2var = PK2var*sqrt((PK.shape[0]-1))
        elif params['sme']:
            print '\tUsing the Standard Deviation of Sample Mean as error'
            PKvar = PKvar/sqrt(PK.shape[0])
            PK2var = PK2var/sqrt(PK.shape[0])

        print 
        print '===power before cal==='
        print PKmean
        print PK2mean
        print PKvar
        print PK2var
            

        #print 
        #print '===power before cal==='
        #print PKmean
        #print PKvar
        if params['cldir']!='':
            #print 
            #print '=== cal==='
            #print self.B
            PKmean = PKmean*self.B
            PKvar  = PKvar*self.B
            PK2mean = PK2mean*self.B2
            PK2var  = PK2var*self.B2
        #print 
        #print '===power after cal==='
        #print PKmean
        #print PKvar
        print 
        print '===power after cal==='
        print PKmean
        print PK2mean
        print PKvar
        print PK2var

        kunit = 2.*pi/(params['boxunit'])
        if (kmin==-1) or (kmax==-1):
            k = np.logspace(
                log10(1./params['boxshape'][0]), log10(sqrt(3)), num=kbn+1)
        else:
            kmin = kmin
            kmax = kmax
            k = np.logspace(log10(kmin), log10(kmax), num=kbn+1)
        #k = k*2.*pi/params['boxunit']
        k = k[:-1]
        k2= np.zeros(shape=(2, kbn))
        k2[0] = k
        k2[1] = k
        sp.save(params['output_root']+resultf+'_k_combined', k)
        sp.save(params['output_root']+resultf+'_k2_combined', k2)
        if FKPweight:
            if params['jkerror']:
                sp.save(params['output_root']+resultf+'_p_var_jk_fkp', PKvar)
                sp.save(params['output_root']+resultf+'_p2_var_jk_fkp', PK2var)
            else:
                sp.save(params['output_root']+resultf+'_p_var_combined_fkp', PKvar)
                sp.save(params['output_root']+resultf+'_p2_var_combined_fkp', PK2var)

                sp.save(params['output_root']+resultf+'_p_combined_fkp', PKmean)
                sp.save(params['output_root']+resultf+'_p2_combined_fkp', PK2mean)

                sp.save(params['output_root']+resultf+'_n_var_combined_fkp', knvar)
                sp.save(params['output_root']+resultf+'_n2_var_combined_fkp', kn2var)

                sp.save(params['output_root']+resultf+'_n_combined_fkp', knmean)
                sp.save(params['output_root']+resultf+'_n2_combined_fkp', kn2mean)
        else:
            if params['jkerror']:
                sp.save(params['output_root']+resultf+'_p_var_jk', PKvar)
                sp.save(params['output_root']+resultf+'_p2_var_jk', PK2var)
            else:
                sp.save(params['output_root']+resultf+'_p_var_combined', PKvar)
                sp.save(params['output_root']+resultf+'_p2_var_combined', PK2var)

                sp.save(params['output_root']+resultf+'_p_combined', PKmean)
                sp.save(params['output_root']+resultf+'_p2_combined', PK2mean)
    
                sp.save(params['output_root']+resultf+'_n_var_combined', knvar)
                sp.save(params['output_root']+resultf+'_n2_var_combined', kn2var)

                sp.save(params['output_root']+resultf+'_n_combined', knmean)
                sp.save(params['output_root']+resultf+'_n2_combined', kn2mean)
    
#   def findbias(self, hr, last, cllist, B):
#
#       def findidx(aa):
#           idx=('A', 'B', 'C', 'D')
#           for s in idx:
#               if aa.find(s)!=-1:
#                   return s
#           return 'NO'
#
#       bias = np.ones(B.shape[1])
#
#       for i in range(len(hr)):
#           idx = findidx(hr[i]) + findidx(last[i])
#           print idx
#           if idx=='NONO':
#               print 'no bias found'
#               continue
#           clidx = cllist.index(idx)
#           bias = bias*np.sqrt(B[clidx])
#           #print bias
#       return bias


    def process_map(self, mapnum, rank):
        params = self.params
        #params['hr'] = (params['hrlist'][mapnum][0],params['hrlist'][mapnum][1])
        #params['last'] =(params['ltlist'][mapnum][0],params['ltlist'][mapnum][1])
        params['imap_pair']=(params['imap_list'][mapnum][0],
                                    params['imap_list'][mapnum][1])

        params['nmap_pair']=(params['nmap_list'][mapnum][0],
                                    params['nmap_list'][mapnum][1])
        if (len(params['mmap_list'])!=0):
            params['mmap_pair']=(params['mmap_list'][mapnum][0],
                                        params['mmap_list'][mapnum][1])
            
        PK, kn, k, PK2, kn2, k2 = self.GetPower()

        print "PK:"
        print PK
        print "PK2:"
        print PK2
        print

        self.q.put_nowait(PK)
        self.qn.put_nowait(kn)
        self.q2.put_nowait(PK2)
        self.qn2.put_nowait(kn2)

#       jkbin = int(params['jknumber']/size)
#       jk = np.array(range(jkbin))
#       jk = (rank*jkbin)+jk
#
#
#       kbn = params['kbinNum']
#       PK = np.zeros(shape=(jk.shape[0],kbn))
#       num = 0
#       for i in jk:
#           params['mid'] = ('jk'+str(i)+mid[0], 'jk'+str(i)+mid[1])
#           #params['mid'] = ('jk'+str(i)+mid[0], mid[1])
#           #print params['mid'][0]
#           kiyopy.utils.mkparents(params['output_root'])
#           inifile = params['output_root']+ 'rank' + str(rank) +'params.ini'
#           parse_ini.write_params(params, inifile ,prefix='pk_')
#           PK[num] = mkpower.PowerSpectrumMaker(
#               inifile, feedback=self.feedback).execute()
#           num = num + 1
#       
#       if rank !=0 :
#           comm.send(PK, dest=0, tag=11)
#       if rank ==0:
#           print 'Calculate the error!'
#           for i in range(1,size):
#               print 'Receive ' + str(i)
#               PK = np.append(PK,comm.recv(source=i, tag=11),axis=0)
#
#           if FKPweight:
#               sp.save(params['output_root']+\
#                   'PKjk_fkp_' + resultf, PK)
#           else:
#               sp.save(params['output_root']+'PKjk_' + resultf, PK)
#           #PKmean = sp.load(params['input_root'] + 'PK.npy')
#           PKmean = PK.mean(axis=0)
#           PK[:] = (PK[:]-PKmean)**2
#           PKvar = np.sum(PK, axis=0)
#           PKvar = PKvar*(params['jknumber']-1)/params['jknumber']
#           PKvar = np.sqrt(PKvar)
#           print PKvar
#           if FKPweight:
#               sp.save(params['output_root']+\
#                   'PKvar_fkp_' + resultf, PKvar)
#           else:
#               sp.save(params['output_root']+\
#                   'PKvar_' + resultf, PKvar)

if __name__ == '__main__':
    import sys
    if len(sys.argv)==2 :
        PowerSpectrumMaker(str(sys.argv[1])).execute()
    elif len(sys.argv)>2 :
        print 'Maximun one argument, a parameter file name.'
    else :
        PowerSpectrumMaker().execute()
