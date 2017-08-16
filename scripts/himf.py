import numpy as np
import pylab as pl
from scipy import integrate
from sympy import integrate as inte

# parameters currently from ALFALFA (Martin et al 2010)
mass = (10.**(9.91))
param = {'theta': 0.006, 'alpha': -1.25, 'm': mass}

x = np.logspace((6), (25), 1000000)
print x[0:10]
mf = param['theta']*((x) / param['m'])**(1+(param['alpha'])) * 
             np.exp(-x / param['m']) / (-1-param['alpha'])

mbins = np.logspace(6.2, 11, 25)

print integrate.simps(mf, x) / 10**9.5
int = inte(mf, x) / x

pl.plot(x, (mf))
pl.plot(x, int)
pl.xlim([2.512e6, 1.e11])
pl.ylim([1.58e-6, 3.16e3])
pl.xscale('log')
pl.yscale('log')
pl.show()
