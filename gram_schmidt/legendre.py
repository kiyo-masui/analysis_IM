#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre as leg

def p(n):
    x = np.linspace(-1, 1, 5000)
    return leg(n)(x)

plt.figure(figsize=(8,7))
x = np.linspace(-1, 1, 5000)
plt.plot(x, p(0), label='0')
plt.plot(x, p(1), label='1')
plt.plot(x, p(2), label='2')
plt.plot(x, p(3), label='3')
plt.plot(x, p(4), label='4')
plt.plot(x, p(10), label='10')
plt.plot(x, p(20), label='20')
plt.plot(x, p(30), label='30')
plt.plot(x, p(200), label='200')
plt.legend(frameon=False)

plt.savefig('./png/lengdre.png')
plt.show()

