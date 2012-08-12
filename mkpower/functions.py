"""
    Public Function are defined here


"""
# public module
import numpy as np
import scipy as sp
from scipy import integrate
from math import *

# kiyo module
from core import algebra
from utils import distance

# li module
import MakePower

pi = 3.1415926
deg2rad = pi/180.
#----------------------------------------------------------------------#
def E(z, omegam, omegal):
    return (1.+z)**2*(1.+omegam*z)-z*(2.+z)*omegal

def funcdl(z, omegam, omegal):
    func = lambda z, omegam, omegal: (E(z, omegam, omegal))**(-0.5)
    dl, dlerr = integrate.quad(func, 0, z, args=(omegam, omegal))

    if omegam+omegal>1. :
        k = (omegam+omegal-1.)**(0.5)
        return sin(k*dl)/k
    elif omegam+omegal<1.:
        k = (1.-omegam-omegal)**(0.5)
        return sinh(k*dl)/k
    elif omegam+omegal==1.:
        return dl

def xyzv(ra, de, r, ra0=0.):
    x = r*sin(0.5*pi-de)*cos(ra-ra0)
    y = r*sin(0.5*pi-de)*sin(ra-ra0)
    z = r*cos(0.5*pi-de)
    v = r**2*sin(0.5*pi-de)
    return x, y, z, v

def fq2r(freq, freq0=1.42e9 , c_H0 = 2.99e3, Omegam=0.27, Omegal=0.73):
    """change the freq to distence"""
    zz =  freq0/freq - 1.
    cosmo = {'omega_M_0' : Omegam, 'omega_lambda_0' : Omegal, 'h' : 0.72}
    cosmo = distance.set_omega_k_0(cosmo)
    zz = 0.72*distance.comoving_distance_transverse(zz, **cosmo)
    #for i in range(0, zz.shape[0]):
    #   zz[i] = c_H0*funcdl(zz[i], Omegam, Omegal)
    return zz

def discrete(params, array):
    """discrete the data pixel into small size"""
    newarray = sp.zeros(params['discrete']*(array.shape[0]-1)+1)
    for i in range(0, array.shape[0]-1):
        delta = (array[i+1]-array[i])/float(params['discrete'])
        for j in range(0, params['discrete']):
            newarray[i*params['discrete']+j] = array[i] + j*delta
    newarray[-1] = array[-1]
    return newarray
#----------------------------------------------------------------------#

def getedge(map):
    ra   = map.get_axis('ra')*deg2rad
    ra   = ra - ra[round(ra.shape[0]/2.)]
    dec  = map.get_axis('dec')*deg2rad
    freq = map.get_axis('freq')
    r    = fq2r(freq)

    r.shape   = r.shape + (1, 1)
    ra.shape  = (1,) + ra.shape + (1,)
    dec.shape = (1, 1) + dec.shape
    radec_map = np.zeros((3,)+map.shape)
    radec_map[0,...] = r
    radec_map[1,...] = ra
    radec_map[2,...] = dec
    xyz_map = np.zeros((3,)+map.shape)
    xyz_map[0,...]=radec_map[0,...]*np.cos(radec_map[2,...])*np.cos(radec_map[1,...])
    xyz_map[1,...]=radec_map[0,...]*np.cos(radec_map[2,...])*np.sin(radec_map[1,...])
    xyz_map[2,...]=radec_map[0,...]*np.sin(radec_map[2,...])

    x = xyz_map[0].flatten()
    y = xyz_map[1].flatten()
    z = xyz_map[2].flatten()

    xrange = (floor(x.min()), floor(x.max())+1.)
    yrange = (floor(y.min()), floor(y.max())+1.)
    zrange = (floor(z.min()), floor(z.max())+1.)

    return xrange, yrange, zrange


#----------------------------------------------------------------------#
def fill(params, imap, nmap, mmap=None):
    """
    Function that used to fill the fftbox with the intensity map
    
    params : the params dict for each module
    imap : the direction to the intensity maps
    nmap : the direction to the noise maps
    mmap : the direction to the mock maps

    It will return the fft box and nbox.
    box is for the intensity maps, while
    nbox is for the noise intensity maps.
    If the mmap!=None, it also return mbox
    which for the mock maps
    """
    #params = self.params
    
    mapshape = np.array(imap.shape)

    r  = fq2r(imap.get_axis('freq'))
    ra = imap.get_axis('ra')*deg2rad
    de = imap.get_axis('dec')*deg2rad
    ra0= ra[int(ra.shape[0]/2)]
    ra = ra - ra0
    dra= ra.ptp()/ra.shape[0]
    dde= de.ptp()/de.shape[0]


    #ra_far = 0
    #de_far = 0
    #if ra.min()*ra.max()>0:
    #   if fabs(ra.min())<fabs(ra.max()):
    #       ra_far = ra.min()
    #   else:
    #       ra_far = ra.max()
    #if de.min()*de.max()>0:
    #   if fabs(de.min())<fabs(de.max()):
    #       de_far = de.min()
    #   else:
    #       de_far = de.max()
    #point = []
    #for i in range(3):
    #   point.append([xyzv(ra.min(), de.min(), r.min())[i], 
   #                 xyzv(ra.max(), de.min(), r.min())[i],
   #                 xyzv(ra.min(), de.max(), r.min())[i],
   #                 xyzv(ra.max(), de.max(), r.min())[i],
   #                 xyzv(ra.min(), de.min(), r.max())[i],
   #                 xyzv(ra.max(), de.min(), r.max())[i],
   #                 xyzv(ra.min(), de.max(), r.max())[i],
   #                 xyzv(ra.max(), de.max(), r.max())[i],
   #                 xyzv(ra_far,   de_far,   r.max())[i],
    #                     ])
    #   point[i].sort()
    #print point
    #print params['boxshape']
    #print (point[0][-1]-point[0][0])/params['boxshape'][0]
    #print (point[1][-1]-point[1][0])/params['boxshape'][1]
    #print (point[2][-1]-point[2][0])/params['boxshape'][2]
    #return 0

    #print r.min(), r.max()
    #print xyz(ra.min(), de.min(), r.min())
    #print xyz(ra.max(), de.min(), r.min())
    #print xyz(ra.min(), de.max(), r.min())
    #print xyz(ra.max(), de.max(), r.min())
    #print xyz(ra.min(), de.min(), r.max())
    #print xyz(ra.max(), de.min(), r.max())
    #print xyz(ra.min(), de.max(), r.max())
    #print xyz(ra.max(), de.max(), r.max())

    ###return 0

    mapinf = [ra.min(), dra, de.min(), dde]
    mapinf = np.array(mapinf)

    box = algebra.info_array(sp.zeros(params['boxshape']))
    box.axes = ('x','y','z')
    box = algebra.make_vect(box)
    
    box_xrange = params['Xrange']
    box_yrange = params['Yrange']
    box_zrange = params['Zrange']
    box_unit = params['boxunit']
    box_disc = params['discrete']

    print box_xrange, box_yrange, box_zrange, box_unit, 

    box_x = np.arange(box_xrange[0], box_xrange[1], box_unit/box_disc)
    box_y = np.arange(box_yrange[0], box_yrange[1], box_unit/box_disc)
    box_z = np.arange(box_zrange[0], box_zrange[1], box_unit/box_disc)

    #print box_x.shape
    #print box_y.shape
    #print box_z.shape

    boxshape = np.array(box.shape)*int(box_disc)

    boxinf0 = [0, 0, 0]
    boxinf0 = np.array(boxinf0)
    boxinf1 = [boxshape[0], boxshape[1], boxshape[2]]
    boxinf1 = np.array(boxinf1)

    print "MapPrepare: Filling the FFT BOX"
    MakePower.Filling(
        imap, r, mapinf, box, boxinf0, boxinf1, box_x, box_y, box_z)


    nbox = algebra.info_array(sp.zeros(params['boxshape']))
    nbox.axes = ('x','y','z')
    nbox = algebra.make_vect(nbox)

    #nbox = algebra.info_array(sp.ones(params['boxshape']))
    #nbox.axes = ('x','y','z')
    #nbox = algebra.make_vect(nbox)

    #print boxinf1
    MakePower.Filling(
        nmap, r, mapinf, nbox, boxinf0, boxinf1, box_x, box_y, box_z)


    if mmap != None:
        mbox = algebra.info_array(sp.zeros(params['boxshape']))
        mbox.axes = ('x','y','z')
        mbox = algebra.make_vect(mbox)
        MakePower.Filling(
            mmap, r, mapinf, mbox, boxinf0, boxinf1, box_x, box_y, box_z)
        return box, nbox, mbox
    else:
        return box, nbox


#----------------------------------------------------------------------#
def getresultf(params):
    """
    get the resultf 
    """
    resultf = params['resultf']
    if resultf == '':                               
        resultf = 'testresult'
        #resultf = params['hr'][0]                  
        #if len(params['last']) != 0:
        #   resultf = resultf + params['last'][0]
        #resultf = resultf + '-' + params['hr'][1]
        #if len(params['last']) != 0:
        #   resultf = resultf + params['last'][1]
    return resultf

#----------------------------------------------------------------------#
def getmap(imap_fname, nmap_fname, mmap_fname=None):
    """
    get the matrix of intensity map and noise map
    """
    #in_root = params['input_root']

    imap = algebra.load(imap_fname)
    imap = algebra.make_vect(imap)
    #print "--The neam value for imap is:",imap.flatten().mean(),"--"
    #imap = imap - imap.flatten().mean()
    if imap.axes != ('freq', 'ra', 'dec') :
        raise ce.DataError('AXES ERROR!')

    try:
        nmap = algebra.load(nmap_fname)
        nmap = algebra.make_vect(nmap)

        bad = nmap<1.e-5*nmap.flatten().max()
        nmap[bad] = 0.
        non0 = nmap.nonzero()
        #imap[non0] = imap[non0]/nmap[non0]
    except IOError:
        print 'NO Noise File :: Set Noise to One'
        nmap = algebra.info_array(sp.ones(imap.shape))
        nmap.axes = imap.axes
        nmap = algebra.make_vect(nmap)
    nmap.info = imap.info
    if nmap.axes != ('freq', 'ra', 'dec') :
        raise ce.DataError('AXES ERROR!')

    if mmap_fname != None:
        try:
            mmap = algebra.load(mmap_fname)
            mmap = algebra.make_vect(mmap)
        except IOError:
            print 'NO Mock File :: Make it!'
            mmap = algebra.info_array(
                2.*np.random.rand(imap.shape[0],imap.shape[1], imap.shape[2])-0.5)
            mmap.axes = imap.axes
            mmap = algebra.make_vect(mmap)
        
        return imap, nmap, mmap
    else:
        return imap, nmap
    
#----------------------------------------------------------------------#
def grothfunction(z, omegam, omegal):
    func = lambda z, omegam, omegal: (1+z)*(E(z, omegam, omegal))**(-3)
    D, Derr = integrate.quad(func, z, np.inf, args=(omegam, omegal))

    D = D*(5/2)*omegam*E(z, omegam, omegal)

    return D
#----------------------------------------------------------------------#
def getpower_th(params):
    """
        Get the theoretial power spectrum result given by camb
        Return : Tb*Tb * G * P_th
    """
    PKcamb = sp.load(params['input_root']+'pk_camb.npy')
    #PKnonl = sp.load(params['input_root']+'nonlPK_'+resultf+'.npy')
    #PKnonl_k = sp.load(params['input_root']+'k_nonlPK_'+resultf+'.npy')
    OmegaHI= params['OmegaHI']
    Omegam = params['Omegam']
    OmegaL = params['OmegaL']
    z = params['z']

    def get_distortion(b):
        f = (Omegam*(1+z)**3)**0.55
        t = (1.+(2./3.)*(f/b)+(1./5.)*(f/b)**2)
        return t

    #b_opt = 1.2
    #t_opt = get_distortion(b_opt)

    b_opt = 1.2 
    t_opt = 1.
        
    if params['optical']:
        print '\tb=%f'%b_opt, '\tt=%f'%t_opt
        PKcamb[1] = PKcamb[1]*b_opt**2*t_opt**2
        PKcamb[1] = PKcamb[1]*(PKcamb[0]**3)/2./3.1415926/3.1415926
        return PKcamb

    #b_gbt = 1.7
    #t_gbt = get_distortion(b_gbt)

    b_gbt = 1.35
    t_gbt = 1

    # Get the Tb
    a3 = (1+z)**(-3)
    Tb = 0.3*(OmegaHI)*((Omegam + a3*OmegaL)/0.29)**(-0.5)*((1.+z)/2.5)**0.5
    if params['PKunit']=='mK':
        Tb = Tb/1.e-3

    # Get G
    #def getG(z):
    #   xx = ((1.0/Omegam)-1.0)/(1.0+z)**3
    #   num = 1.0 + 1.175*xx + 0.3046*xx**2 + 0.005335*xx**3
    #   den = 1.0 + 1.875*xx + 1.021 *xx**2 + 0.1530  *xx**3

    #   G = (1.0 + xx)**0.5/(1.0+z)*num/den
    #   return G
    #G = getG(z)
    
    ##print G**2*(1+z)**(-2)
    #D  = grothfunction(z, Omegam, OmegaL)
    #D0 = grothfunction(0, Omegam, OmegaL)
    #print (D/D0)**2

    #PKcamb[1] = PKcamb[1]*(D/D0)**2
    #PKcamb[1] = PKcamb[1]*(G**2*(1.+z)**(-2)*6**2)

    if params['cross']: 
        PKcamb[1] = PKcamb[1]*Tb*b_gbt*b_opt*sqrt(t_opt*t_gbt)
    else: 
        print '\tb=%f'%b_gbt, '\tt=%f'%t_gbt
        PKcamb[1] = PKcamb[1]*(Tb**2)*b_gbt**2*t_gbt**2
        #PKcamb[1] = PKcamb[1]*(Tb**2)

    PKcamb[1] = PKcamb[1]*(PKcamb[0]**3)/2./3.1415926/3.1415926

    return PKcamb
#----------------------------------------------------------------------#

