import shelve
import os

import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
from scipy import random, optimize

plt.ion()

path = os.getenv('RAID_PRO')
path += 'eswitzer/GBT/pwrspec/GBT_15hr_map_oldcal_noise_simulations_xspec.shelve'

db = shelve.open(path, 'r')

C_sim = db["stats"]["0modes"]["pk_1d_from_2d_stat"]["cov"]
P_sim = db["stats"]["0modes"]["pk_1d_from_2d_stat"]["mean"]

k = db['k_1d_from_2d']['center']

mask = np.isfinite(P_sim)
P_sim = P_sim[mask]
C_sim = C_sim[mask][:,mask]
k = k[mask]
n = len(k)

# Generate model for Covariance.
C_mod, f = get_C_model(C_sim, 10)

# Make a fake power spectrum.
P_fake = 1.3 * P_sim + 3.1 * np.mean(P_sim) * (k / np.mean(k))**1.4
e, v = linalg.eigh(C_sim)
noise = np.sqrt(e) * random.randn(n)
noise = np.dot(v, noise)
P_fake += noise

# Choose PS and covariance to process.
P = P_fake
C = C_sim
C_inv = linalg.inv(C)

# Process and fit.
fg_residuals = P - P_sim

def PL_model(pars):
    A = pars[0]
    beta = pars[1]
    return A * (k/k[0])**beta
    

def chi2(pars):
    resid = fg_residuals - PL_model(pars)
    return np.dot(resid, np.dot(C_inv, resid))


out = optimize.fmin(chi2, [np.mean(fg_residuals), 2])
print chi2(out) / n

plt.figure()
plt.xscale('log')
plt.yscale('log')
plt.plot(k, fg_residuals, '.')
plt.plot(k, PL_model(out))








A = norm_mat(linalg.inv(C_mod))
plt.figure()
plt.pcolor(k, k, A)
plt.xscale('log')
plt.yscale('log')
plt.colorbar()







def get_C_model(C, m):

    n = C.shape[0]
    if C.shape != (n, n):
        raise ValueError("Must pass square matrix")
    # First normalize.
    norm = np.sqrt(np.diag(C))
    C_norm = C / (norm[None,:] * norm[:,None])
    # Subtract diagonal.
    C_norm.flat[::n + 1] = 0
    # XXX
    #plt.figure()
    #plt.imshow(C_norm, vmin=-.2, vmax=0.5)
    #plt.colorbar()
    
    # Resynthesize the diagonal.
    C_replaced = C_norm.copy()
    for q in range(n + 1):
        for ii in range(5):
            new_diag = truncated_mat(C_replaced, q).flat[::n + 1]
            C_replaced.flat[::n + 1] = new_diag

    # Now take the m modes and put them in the matrix.
    C_out = truncated_mat(C_replaced, m)
    # Reset the diagonal to the right level.
    C_out.flat[::n + 1] = 1
    # XXX
    #plt.figure()
    #plt.imshow(C_out, vmin=-.2, vmax=0.5)
    #plt.colorbar()
    

    fom = mat_resid(C_out, C_norm)
    C_out *= (norm[None,:] * norm[:,None])

    return C_out, fom


def norm_mat(C):
    norm = np.sqrt(np.diag(C))
    return C / (norm[None,:] * norm[:,None])

def truncated_mat(C, m):
    n = C.shape[0]
    e, v = linalg.eigh(C)
    if m > 0:
        out = np.sum(e[None,-m:,None] * v[:,-m:,None]
                       * v.T[None,-m:,:], 1)
    else:
        out = np.zeros((n, n))
    return out


def mat_resid(A, B):

    if A.shape != B.shape:
        raise ValueError("Arrays must be the same shape")

    return np.sum((A - B)**2) / A.size



