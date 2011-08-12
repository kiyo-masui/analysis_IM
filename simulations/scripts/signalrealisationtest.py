import numpy as np

from core import algebra
from map import beam
from simulations import corr, corr21cm, ps_estimation

from utils import cubicspline as cs
from utils import units




def bincov(cov, bins = None):
    l = cov.shape[0]
    xr = np.arange(l)

    bins = xr if bins == None else bins

    rows, cols = np.meshgrid(xr, xr)
    dist = rows - cols

    df1 = dist.flat[np.where(dist.flat >= 0)]
    cf1 = cov.flat[np.where(dist.flat >= 0)]

    ds1 = df1[np.argsort(df1)]
    cs1 = cf1[np.argsort(df1)]

    dg1 = np.digitize(ds1, bins)

    #pdb.set_trace()

    rind = np.insert(np.where(dg1[1:] - dg1[:-1]), 0, -1) + 1

    nr = rind[1:] - rind[:-1]

    csim = np.insert(np.cumsum(cs1), 0, 0.0)

    tbin = csim[rind[1:]] - csim[rind[:-1]]

    return tbin / nr

thetax = 32.0
thetay = 32.0
nx = 128
ny = 128

f1 = 850
f2 = 650
nf = 256

z1 = units.nu21 / f1
z2 = units.nu21 / f2

c1 = cs.LogInterpolater.fromfile("data/ps_z1.5.dat")
kstar = 2.0
ps = lambda k: np.exp(-0.5 * k**2 / kstar**2) * c1(k)



#cr = corr.RedshiftCorrelation(ps_vv = ps, redshift = 1.5)
cr = corr21cm.Corr21cm(ps_vv = ps, redshift = 1.5)


#rf = cr.realisation(1.0, 1.0, 1.95, 2.0, 256, 256, 256)

#cb = cr.realisation_dv([1024.0, 1024.0, 1024.0], [128, 128, 128])

#df = cb[0]
#vf = cb[1]

#rfv = cb[2]

#del cb

#tf = df + vf

#tp, kpar, kperp = ps_estimation.ps_azimuth(df, width=[64.0, 64.0, 64.0], kmodes = True)
#kvec = np.rollaxis(np.array(np.meshgrid(kpar, kperp)), 0, 3)
#pst = ps((kvec**2).sum(axis=2)**0.5)
#mi = np.fft.irfftn(rfv._kweight* 2**0.5)
#tpm, kparm, kperpm = ps_estimation.ps_azimuth(mi, width=[64.0, 64.0, 64.0], kmodes = True)


rf = cr.realisation(z1, z2, thetax, thetay, nf, nx, ny)[::-1,...]

psnw, kpar, kperp = ps_estimation.ps_azimuth(rf, window = False)
psww, kpar, kperp = ps_estimation.ps_azimuth(rf, window = True)

lag0 = np.cov(rf.reshape((256, 128*128)))

bc = bincov(lag0)

a = algebra.make_vect(rf, axis_names = ('freq', 'ra', 'dec'))

a.set_axis_info('freq', (f1+f2)/2.0, (f1-f2)/nf)
a.set_axis_info('ra', 0.0, thetax / nx)
a.set_axis_info('dec', 0.0, thetay / ny)

b = beam.GaussianBeam(width = [0.25, 0.25*f2/f1], freq = [f2, f1])

ab = b.apply(a)













