import numpy as np
from scipy import interpolate

pwrspec_data = np.genfromtxt("simulations/data/wigglez_halofit_z0.8.dat")
(log_k, log_pk) = (np.log(pwrspec_data[:,0]), \
                   np.log(pwrspec_data[:,1]))

pk_test = interpolate.interp1d(pwrspec_data[:,0],
                                pwrspec_data[:,1])

logpk_interp = interpolate.interp1d(log_k, log_pk)
pk_interp = lambda k: np.exp(logpk_interp(np.log(k)))

ktest = pwrspec_data[:,0]

np.min(ktest), np.max(ktest)

rangek = np.array([0.00005, 0.0001, 1, 11.584, 12])
print rangek

rangek[rangek < np.min(ktest)] = np.min(ktest)
rangek[rangek > np.max(ktest)] = np.max(ktest)

print rangek

#for kval in pwrspec_data[:,0]:
#    print kval, pk_interp(kval), pk_interp(kval*1.01), pk_test(kval)
