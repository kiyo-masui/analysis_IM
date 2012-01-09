from numpy import *
from scipy.optimize import leastsq
import pylab

def peval(x,p):
   return( (x/p[0])**p[1])

def residuals(p,y,x,errors):
   err = (y-peval(x,p))/errors
   return err

#wavelength from vla calibrator list
#w_cm = array([90.0,20.0,6.0,3.7,2.0,1.3,0.7])
#w_cm = array([90.0,20.0,6.0])
## convert to MHz
#x = 299.792458/(w_cm/100.0)
# from ned list for 3c48
x = array([8300.00,
5000.00,
5000.00,
5000.00,
5000.00,
5000.00,
4850.00,
4850.00,
4850.00,
4830.00,
2700.00,
2700.00,
2700.00,
2700.00,
2700.00,
2700.00,
2250.00,
1420.00,
1400.00,
1400.00,
1400.00,
1400.00,
1400.00,
1400.00,
1400.00,
750.00,
750.00,
750.00,
750.00,
750.00,
408.00,
408.00,
365.00,
318.00,
178.00,
178.00,
178.00,
178.00,
160.00,
160.00
])
print x
#flux density from vla calibrator list
#flux = array([42.0, 16.5, 5.48,3.25,1.78,1.13,0.64])
#flux = array([42.0, 16.5, 5.48])
#flux density from ned
flux = array([3.33E+00,
5.33E+00,
5.33E+00,
5.48E+00,
5.37E+00,
5.37E+00,
5.73E+00,
5.75E+00,
5.48E+00,
6.08E+00,
9.17E+00,
9.32E+00,
9.10E+00,
8.97E+00,
9.07E+00,
9.07E+00,
1.14E+01,
1.59E+01,
1.59E+01,
1.58E+01,
1.60E+01,
1.56E+01,
1.53E+01,
1.57E+01,
1.57E+01,
2.55E+01,
2.55E+01,
2.44E+01,
2.55E+01,
2.70E+01,
3.87E+01,
3.75E+01,
4.23E+01,
3.74E+01,
5.10E+01,
5.50E+01,
5.99E+01,
5.59E+01,
6.98E+01,
6.77E+01
])
print flux.shape, x.shape
# go to temperature GBT reported as 2.0K/Jy rather flat for low freq, dont care about high will uses gbtidl model
ep = 390e-6
c = 299792458.0
na = 0.71*exp(-(4*pi*x*1e6/c)**2)
to = 0.008 + exp(sqrt(x/1e3))/8000.0
nl = 0.99
KpJy = 2.85 * na * nl *exp(-to)
# dropped the KpJy, just get flux, do convert to temp elsewhere
y = flux
#error = array([1.0,1.0,1.0,1.0,1.0,1.0,1.0])
error = array([1.0,
2.70E-01,
7.00E-02,
5.00E-02,
7.00E-02,
2.70E-01,
8.59E-01,
7.37E-01,
1.00E+00,
1.00E+00,
4.60E-01,
1.00E+00,
7.30E-01,
4.50E-01,
1.00E-01,
1.01E-01,
1.00E+00,
1.00E+00,
5.30E-01,
1.00E+00,
4.81E-01,
2.50E-01,
7.70E-01,
8.00E-01,
3.60E-01,
5.20E-01,
3.20E-01,
1.22E+00,
1.30E+00,
3.40E-01,
3.02E+00,
1.88E+00,
3.45E-01,
1.40E+00,
4.08E+00,
5.50E+00,
3.00E+00,
2.90E+00,
9.10E+00,
1.00E+01])

p0 = array([2000,-0.7])
## leastsq takes function to minimize (residuals), first argument for function( p0 ), additional arguments are args ( y,x, errors)
plsq = leastsq(residuals,p0, args=(y,x,error),full_output=1, maxfev=2000)
fitp = plsq[0]
errp = plsq[1]

dev = -1.0 * residuals(fitp, y, x, error)
chisq = sum(dev**2)
dof = len(y) - len(p0)
# multiply by the residual standard deviation to get a better error estimate of covariance matrix
# but should only do it if chisq/dof > 1
chisq_dof = chisq/dof
if ( chisq_dof > 1):
   perrs = errp*chisq/dof
else:
   perrs = errp
   print "chisq/dof < 1, too small! not sure about error estimates"


print "(r/p0)**p1"
print fitp
print "error matrix in paramters:"
print perrs
print "chisq/dof:"
print (chisq_dof)

xfit = arange(100,8300,50)

yfit = peval(xfit, fitp)

pylab.errorbar(x,y, yerr=error)
pylab.plot(xfit, yfit, 'r-')
pylab.show()
