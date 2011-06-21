import numpy as np

import corr

import cubicspline as cs
import ps_estimation

import pdb


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

c1 = cs.LogInterpolater.fromfile("data/ps.dat")
kstar = 2.0
ps = lambda k: np.exp(-0.5 * k**2 / kstar**2) * c1(k)

c2 = cs.LogInterpolater.fromfile("data/ps_z1.5.dat")
kstar = 0.5
ps2 = lambda k: np.exp(-0.5 * k**2 / kstar**2) * c2(k)



cr = corr.RedshiftCorrelation(ps_vv = ps2)


#rf = cr.realisation(1.0, 1.0, 1.95, 2.0, 256, 256, 256)

cb = cr.realisation_dv([1024.0, 1024.0, 1024.0], [128, 128, 128])

df = cb[0]
#vf = cb[1]

#rfv = cb[2]

#del cb

#tf = df + vf

tp, kpar, kperp = ps_estimation.ps_azimuth(df, width=[64.0, 64.0, 64.0], kmodes = True)
#kvec = np.rollaxis(np.array(np.meshgrid(kpar, kperp)), 0, 3)
#pst = ps((kvec**2).sum(axis=2)**0.5)
#mi = np.fft.irfftn(rfv._kweight* 2**0.5)
#tpm, kparm, kperpm = ps_estimation.ps_azimuth(mi, width=[64.0, 64.0, 64.0], kmodes = True)


rf = cr.realisation(32.0, 32.0, 0.5, 1.0, 128, 128, 256)

psnw, kpar, kperp = ps_estimation.ps_azimuth(rf, window = False)
psww, kpar, kperp = ps_estimation.ps_azimuth(rf, window = True)

lag0 = np.cov(rf.reshape((256, 128*128)))

bc = bincov(lag0)







