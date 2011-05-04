


import numpy as np
from gaussianfield import RandomField

class _LofarGSDE_3D(RandomField):

    delta = -4.0

    def powerspectrum(self, karray):
        r"""Power law power spectrum"""

        ps = (karray**2).sum(axis=3)**(self.delta / 2.0)
        ps[0,0,0] = 0.0

        return ps

class LofarGSDE(object):

    widthx = 5.0
    widthy = 5.0
    
    numx = 128
    numy = 128

    numf = 128

    nu_0 = 325.0

    nu_l = 120.0
    nu_h = 325.0
    
    correlated = False

    A_amp = 20
    A_std = A_amp * 0.02
    
    beta_mean = -2.55
    beta_std = 0.1

    alpha = -2.7
    
    def getfield(self):
        r"""Lofar synchrotron."""

        numz = int((self.numx + self.numy)/2)

        # Set up 3D field generator
        npix = [self.numx, self.numy, numz]
        wsize = [5.0 / self.widthx, 5.0 / self.widthy, 1.0]
        lf = _LofarGSDE_3D(npix = npix, wsize = wsize)
        #lf.delta = self.alpha - 1.0
        lf.delta = self.alpha
        
        A = lf.getfield()
        beta = A if self.correlated else lf.getfield()

        A = ((1.0*self.A_amp) / numz) + A * (self.A_std / A.sum(axis = 2).std())
        beta = self.beta_mean + beta * (self.beta_std / beta.std())

        freq = np.linspace(self.nu_l, self.nu_h, self.numf) / self.nu_0

        Tb = np.zeros((self.numf, self.numx, self.numy))

        for i in range(self.numf):
            Tb[i,:,:] = (A * freq[i]**beta).sum(axis=2)

        return Tb


