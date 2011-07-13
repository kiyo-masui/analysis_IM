

import numpy as np

import cubicspline as cs

dmfunc = np.vectorize(cs.Interpolater.fromfile('data/dmcorr.dat').value)

def c_syn(l, nu_1, nu_2):

    A = 700.0
    beta = 2.4
    alpha = 2.80
    zeta = 1.0

    nu_f = 130.0

    c = A * (1e3 / l)**beta * (nu_f**2 / (nu_1 * nu_2))**alpha *  np.exp(-0.5 * np.log(nu_1 / nu_2)**2 / zeta**2)

    return c

def cf_syn(nu_1, nu_2 = None):

    if(nu_2 == None):
        l = np.arange(1, 5000, dtype=np.float64)[:,np.newaxis]
        n1 = nu_1[np.newaxis,:] 
        n2 = nu_1[np.newaxis,:] 
    else:
        l = np.arange(1, 5000, dtype=np.float64)[:,np.newaxis,np.newaxis]
        n1 = nu_1[np.newaxis,:,np.newaxis] 
        n2 = nu_2[np.newaxis,np.newaxis,:] 
    
    c = c_syn(l, n1, n2)

    cf = (c * (2.0*l + 1.0)).sum(axis=0) / (4*np.pi)

    return cf

fgs = [ ('syn', 700.0, 2.80, 2.4, 4.0),
        ('ps',   57.0, 2.07, 1.1, 1.0),
        ('eff', 0.014, 2.10, 1.0, 35.0),
        ('gff', 0.088, 2.15, 3.0, 35.0) ]


def c_all(l, nu_1, nu_2):

    c = None
    nu_f = 130.0

    for name, A, alpha, beta, zeta in fgs:
        c1 = A * (1e3 / l)**beta * (nu_f**2 / (nu_1 * nu_2))**alpha *  np.exp(-0.5 * np.log(nu_1 / nu_2)**2 / zeta**2)

        c = c1 if c == None else c+c1

    return c

def cf_all(nu_1, nu_2 = None):

    if(nu_2 == None):
        l = np.arange(1, 5000, dtype=np.float64)[:,np.newaxis]
        n1 = nu_1[np.newaxis,:] 
        n2 = nu_1[np.newaxis,:] 
    else:
        l = np.arange(1, 5000, dtype=np.float64)[:,np.newaxis,np.newaxis]
        n1 = nu_1[np.newaxis,:,np.newaxis] 
        n2 = nu_2[np.newaxis,np.newaxis,:] 
    
    c = c_all(l, n1, n2)

    cf = (c * (2.0*l + 1.0)).sum(axis=0) / (4*np.pi)

    return cf


def cf_all2(nu_1, nu_2 = None):

    def ct(l):
        n1 = nu_1[np.newaxis,:,np.newaxis] 
        n2 = nu_2[np.newaxis,np.newaxis,:]
        la = l[:,np.newaxis,np.newaxis]
        return (c_all(la, n1, n2) * (2.0*la + 1.0)).sum(axis=0) / (4*np.pi)

    cv = np.array([ct(ls) for ls in np.array_split(np.arange(1, 5000, dtype=np.float64), 20)])
    #print cv

    return cv.sum(axis=0)



def dmcorr(nu_1, nu_2 = None):

    n1 = nu_1[:,np.newaxis]
    if(nu_2 == None):
        nu_2 = nu_1
    n2 = nu_2[np.newaxis,:]

    piv = 5.95 * np.abs(n1 - n2)

    return 0.3**2 * dmfunc(piv) / (2.5**2)

    
