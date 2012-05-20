#! /usr/bin/env python

import numpy as np
from scipy.special import legendre as leg
from scipy.integrate import romberg
from scipy.integrate import simps
from scipy.integrate import quadrature
import matplotlib.pyplot as plt
from math import *

def p(n, fn=512):
    x = np.linspace(-1, 1, fn)
    return np.polyval(leg(n), x)

def fdot(u, v):
    #return romberg(lambda x: u(x)*v(x), -1, 1)
    fn = len(u)
    x = np.linspace(-1, 1, fn)
    return simps(u*v, x)
    #return quadrature(lambda x: u(x)*v(x), -1, 1)[0]

def proj(u, v):
    r'''
    projection operator: 
    project u to v
    u, v should be array
    '''

    return fdot(u, v)/fdot(v, v)*v

def gs(v):
    r'''
    gram-schmidt operator:
    it will generate a new base using legendre 
    v : 2d array, old bases

    return: the new base
    '''

    vn = v.shape[0]
    x = np.linspace(-1, 1, v.shape[1])
    u = leg(vn)(x)
    for i in range(vn):
        u = u - proj(u, v[i])
    return u

def gs_array(v):
    r'''
    same as gs(v), but the return is the whole base array
    '''

    u = gs(v)
    v = np.resize(v, (v.shape[0]+1, v.shape[1]))
    v[-1] = u
    return v
    


if __name__=='__main__':
    
    x = np.linspace(-1, 1, 500)
    num = 30
    v = np.zeros((num,len(x)))
    for i in range(num):
        v[i] = leg(i)(x)

    plt.figure(figsize=(8,7))
    #plt.plot(x, gs(v), label='gram-schmidt')

    for i in range(31):
        v = gs_array(v)

    plt.subplot(311)
    plt.plot(x, v[30], label='gs legendre P(30)')

    plt.plot(x, leg(num)(x), label='scipy legendre P(30)')
    plt.ylim(ymax=1, ymin=-1)
    plt.legend()
    
    plt.subplot(312)
    plt.plot(x, v[35], label='gs legendre P(35)')

    plt.plot(x, leg(num+5)(x), label='scipy legendre P(35)')
    plt.ylim(ymax=1, ymin=-1)
    plt.legend()

    plt.subplot(313)
    plt.plot(x, v[60], label='gs legendre P(60)')

    plt.plot(x, leg(num+30)(x), label='scipy legendre P(60)')
    plt.ylim(ymax=1, ymin=-1)
    plt.legend()

    plt.savefig('./png/comp.png', format='png')
    plt.show()


