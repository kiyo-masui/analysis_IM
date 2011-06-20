import numpy as np

import corr

import cubicspline as cs
import ps_estimation

ps = cs.LogInterpolater.fromfile("data/ps.dat")

cr = corr.RedshiftCorrelation(ps_vv = ps)


#rf = cr.realisation(1.0, 1.0, 1.95, 2.0, 256, 256, 256)

#cb = cr.realisation_dv([256.0, 64.0, 64.0], [128, 128, 128])

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


rf, rcube = cr.realisation(3.0, 3.0, 0.6, 1.1, 128, 128, 128)







