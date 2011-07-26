import pylab
from numpy import *
import scipy as sp

prefix = '66-68'

mueller_params = loadtxt('mueller_params_calc.txt')
m_err = loadtxt('mueller_params_error.txt')
size = len(mueller_params[:,0])
#print mueller_params[0,:] #is a set of values for highest frequency
freq = mueller_params[:,0]
deltaG = mueller_params[:,1]
alpha = mueller_params[:,2]
psi = mueller_params[:,3]
phi = mueller_params[:,4]
epsilon = mueller_params[:,5]
chi = mueller_params[:,6]
flux = mueller_params[:,7]
#Q = mueller_params[:,6]
#U = mueller_params[:,7]
#chi = mueller_params[:,8]
#Q = sp.zeros(size)
#for ii in range(0,size):
#    Q[ii] = abs(mueller_params[ii,6])
#U = sp.zeros(size)
#for ii in range(0,size):
#    U[ii] = abs(mueller_params[ii,7])

heiles_params = sp.zeros((8,8))
heiles_error = sp.zeros((8,8))

heiles_params[:,1] = [0.676, 0.545, -0.013, -0.398, -0.299, -0.217, -0.116, -0.182]
heiles_error[:,1] = [0.026, 0.025, 0.030, 0.025, 0.021, 0.015, 0.015, 0.037]
heiles_params[:,2] = [-4.9, -1.6, -3.2, 10.0, -2.5, 0.2, 1.7, -1.3]
heiles_error[:,2] = [17.9, 13.9, 18.0, 9.0, 9.0, 8.8, 10.3, 20.9]
heiles_params[:,3] = [134.1, 170.1, 157.8, 174.0, -177.7, -166.7, 167.2, 19.0]
heiles_error[:,3] = [36.4, 27.8, 36.1, 19.3, 18.1, 17.5, 20.6, 41.9]
heiles_params[:,4] = [42.9, 6.2, 18.7, 2.9, -0.9, -14.8, 4.7, 143.4]
heiles_error[:,4] = [42.5, 35.5, 45.1, 24.5, 23.8, 22.7, 26.2, 75.3]
heiles_params[:,5] = [0.017, 0.016, 0.016, 0.024, 0.020, 0.015, 0.014, 0.008]
heiles_error[:,5] = [0.007, 0.006, 0.008, 0.006, 0.005, 0.004, 0.004, 0.009]
heiles_params[:,6] = [-0.002, -0.004, -0.002, 0.014, 0.010, 0.002, 0.002, 0.018]
heiles_error[:,6] = [0.009, 0.009, 0.011, 0.010, 0.008, 0.005, 0.005, 0.013]
heiles_params[:,7] = [0.021, 0.026, 0.024, 0.037, 0.032, 0.025, 0.022, 0.018]
heiles_error[:,7] = [0.009, 0.009, 0.011, 0.009, 0.007, 0.005, 0.005, 0.013]
heiles_params[:,0] = [904, 874, 844, 814, 784, 754, 724, 694]

#deltaG plot
pylab.plot(freq, deltaG,label='deltaG')
#pylab.errorbar(freq,deltaG,m_err[:,1],label = 'deltaG')
#pylab.errorbar(heiles_params[:,0],heiles_params[:,1],heiles_error[:,1],label='h_deltaG')
pylab.scatter(freq,deltaG,label='deltaG values')
deltaG_ave = mean(deltaG)
#pylab.axhline(y=deltaG_ave,color='y',label='mean_deltaG')
pylab.ylim(-0.5,0.5)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_params_deltaG_'+prefix+'.png')
pylab.clf()

#epsilon plot
pylab.plot(freq, epsilon,label='epsilon')
#pylab.errorbar(freq,epsilon,m_err[:,5],label = 'epsilon')
#pylab.errorbar(heiles_params[:,0],heiles_params[:,5],heiles_error[:,5],label='h_epsilon')
pylab.scatter(freq, epsilon, label='epsilon values')
epsilon_ave = mean(epsilon)
#pylab.axhline(y=epsilon_ave,color='m',label='mean_epsilon')
pylab.ylim(-0.05,0.05)
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.legend()
pylab.savefig('mueller_params_epsilon_'+prefix+'.png')
pylab.clf()


#Alpha Plots
pylab.plot(freq, alpha,label='generated')
#pylab.errorbar(freq,alpha,m_err[:,2],label='generated')
#pylab.errorbar(heiles_params[:,0],heiles_params[:,2],heiles_error[:,2],label='heiles')
pylab.scatter(freq, alpha, label='data points')
alpha_ave = mean(alpha)
pylab.axhline(y=90,color='m',label='90 degree rotated linear feed')
pylab.xlim(freq[-1],freq[0])
pylab.ylim(0,180)
pylab.xlabel('frequency')
pylab.ylabel('alpha (degrees)')
pylab.legend()
pylab.savefig('mueller_params_alpha_'+prefix+'.png')
pylab.clf()

#Psi Plots
pylab.plot(freq, psi,label='generated')
#pylab.errorbar(freq, psi, m_err[:,3], label='generated')
#pylab.errorbar(heiles_params[:,0],heiles_params[:,3],heiles_error[:,3],label='heiles')
pylab.scatter(freq, psi, label='data points')
psi_ave = mean(psi)
#pylab.axhline(y=psi_ave,color='m',label='mean_gen')
pylab.xlim(freq[-1],freq[0])
pylab.ylim(-180,180)
pylab.xlabel('frequency')
pylab.ylabel('psi (degrees)')
pylab.legend()
pylab.savefig('mueller_params_psi_'+prefix+'.png')
pylab.clf()

#Phi Plots
pylab.plot(freq, phi,label='generated')
#pylab.errorbar(freq, phi, label = 'generated')
#pylab.errorbar(heiles_params[:,0],heiles_params[:,4],heiles_error[:,4],label='heiles')
pylab.scatter(freq, phi, label='data points')
phi_ave = mean(phi)
#pylab.axhline(y=phi_ave,color='m',label='mean_gen')
pylab.xlim(freq[-1],freq[0])
pylab.ylim(-180,180) 
pylab.xlabel('freqency')
pylab.ylabel('phi (degrees)')
pylab.legend()
pylab.savefig('mueller_params_phi_'+prefix+'.png')
pylab.clf()

#Chi Plots
pylab.plot(freq, chi,label='generated')
#pylab.errorbar(freq, chi, label = 'generated')
#pylab.errorbar(heiles_params[:,0],heiles_params[:,4],heiles_error[:,4],label='heiles')
pylab.scatter(freq, chi, label='data points')
chi_ave = mean(chi) 
#pylab.axhline(y=chi_ave,color='m',label='mean_gen')
pylab.xlim(freq[-1],freq[0])
pylab.ylim(-180,180)
pylab.xlabel('freqency')
pylab.ylabel('chi (degrees)')
pylab.legend()
pylab.savefig('mueller_params_chi_'+prefix+'.png')
pylab.clf()

#Flux plot
pylab.plot(freq,flux,label='generated')
pylab.scatter(freq,flux,label='data points')
pylab.xlim(freq[-1],freq[0])
pylab.xlabel('frequency')
pylab.ylabel('flux')
pylab.legend()
pylab.savefig('mueller_params_flux_'+prefix+'.png')
pylab.clf()

# Q, U plots
#pylab.plot(freq, Q,label='Q_generated')
#pylab.errorbar(freq,Q,m_err[:,6],label='Q_generated')
#pylab.scatter(freq,Q,color='b',label='Q_gen data points')
#pylab.plot(freq, U,label='U_generated')
#pylab.scatter(freq,U,color='g',label='U_gen data points')
#pylab.errorbar(freq,U,m_err[:,7],label='U_generated')
#pylab.errorbar(heiles_params[:,0],heiles_params[:,6],heiles_error[:,6],label='Q_heiles')
#pylab.errorbar(heiles_params[:,0],heiles_params[:,7],heiles_error[:,7],label='U_heiles')
#Q_ave = mean(Q)
#pylab.axhline(y=Q_ave,color='y',label='mean_Q')
#U_ave = mean(U)
#pylab.axhline(y=U_ave,color='m',label='mean_U')
#pylab.xlim(freq[-1],freq[0])
#pylab.ylim(0,0.1) 
#pylab.legend()
#pylab.xlabel('frequency')
#pylab.ylabel('fractional polarization')
#pylab.savefig('mueller_params_Q_U_'+prefix+'.png')
#pylab.clf()

#ratio = sp.zeros(size)
#for ii in range(0,size):
#    ratio[ii] = abs(mueller_params[ii,6]/mueller_params[ii,7])

#Plot Ratio of Q/U (should be a constant)
#pylab.plot(freq, ratio)
#pylab.ylim([0,5])
#pylab.xlim(freq[-1],freq[0])
#pylab.xlabel('frequency')
#pylab.ylabel('ratio of Q/U')
#pylab.savefig('mueller_params_ratio_'+prefix+'.png')
#pylab.clf()

