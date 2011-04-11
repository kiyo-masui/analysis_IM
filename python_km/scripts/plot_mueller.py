import pylab
from numpy import *
import scipy as sp

mueller_params = loadtxt('mueller_params_calc.txt')
#print mueller_params[0,:] #is a set of values for highest frequency
freq = mueller_params[:,0]
deltaG = mueller_params[:,1]
alpha = mueller_params[:,2]
psi = mueller_params[:,3]
phi = mueller_params[:,4]
epsilon = mueller_params[:,5]
Q = sp.zeros(260)
for ii in range(0,260):
    Q[ii] = abs(mueller_params[ii,6])
U = sp.zeros(260)
for ii in range(0,260):
    U[ii] = abs(mueller_params[ii,7])

pylab.plot(freq, deltaG)
pylab.plot(freq, epsilon)
pylab.axis('auto')
pylab.legend(('deltaG', 'epsilon'))
pylab.xlabel('frequency')
pylab.savefig('mueller_params_deltaG_epsilon_39-40.png')
pylab.clf()

pylab.plot(freq, alpha)
pylab.axis('auto')
pylab.xlabel('frequency')
pylab.ylabel('alpha (degrees)')
pylab.savefig('mueller_params_alpha_39-40.png')
pylab.clf()

pylab.plot(freq, psi)
pylab.axis('auto')
pylab.xlabel('frequency')
pylab.ylabel('psi (degrees)')
pylab.savefig('mueller_params_psi_39-40.png')
pylab.clf()

pylab.plot(freq, phi)
pylab.axis('auto')
pylab.xlabel('freqency')
pylab.ylabel('phi (degrees)')
pylab.savefig('mueller_params_phi_39-40.png')
pylab.clf()

pylab.plot(freq, Q)
pylab.plot(freq, U)
pylab.axis('auto')
pylab.legend(('Q_src', 'U_src'))
pylab.xlabel('frequency')
pylab.ylabel('fractional polarization')
pylab.savefig('mueller_params_Q_U_39-40.png')
pylab.clf()

ratio = sp.zeros(260)
for ii in range(0,260):
    ratio[ii] = abs(mueller_params[ii,6]/mueller_params[ii,7])

pylab.plot(freq, ratio)
pylab.axis([650,950,0,5])
pylab.xlabel('frequency')
pylab.ylabel('ratio of Q/U')
pylab.savefig('mueller_params_ratio_39-40.png')
pylab.clf()

