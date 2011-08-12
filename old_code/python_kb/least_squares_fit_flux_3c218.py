from numpy import *
from scipy.optimize import leastsq
import pylab

def peval(x,p):
   return( (x/p[0])**p[1])

def residuals(p,y,x,errors):
   err = (y-peval(x,p))/errors
   return err

#x is frequency in megahertz
#x = array([160.0,160.0, 365.0,408.0,408.0,468.0,635.0,750.0,750.0,960.0,1400.0,1400.0, 1400.0,1410.0,1410.0,2650.0,2700.0,2700.0,2700.0,4850.0,5000.0,5010.0,5010.0,8000.0])
#x = array([160.0,160.0, 408.0,408.0,468.0,635.0,750.0,750.0,960.0,1400.0,1400.0, 1400.0,1410.0,1410.0,2650.0,2700.0,2700.0,2700.0,4850.0,5000.0,5010.0,5010.0,8000.0])
x = array([8000.00,
5010.00,
5010.00,
5000.00,
4850.00,
2700.00,
2700.00,
2700.00,
2650.00,
1410.00,
1410.00,
1400.00,
1400.00,
1400.00,
960.00,
750.00,
750.00,
635.00,
468.00,
408.00,
408.00,
160.00,
160.00])
#print x
#flux density took out the 365MHz measurement
#flux = array([243.0,245.8,72.939,132.0,132.0,114.99,97.09,78.92,83.58,65.44,40.8,43.4,45.05,44.55,43.5,24.38,24,23.5,23.5,13.982,13.1,13.91,14.04,8.7])
#flux = array([243.0,245.8,132.0,132.0,114.99,97.09,78.92,83.58,65.44,40.8,43.4,45.05,44.55,43.5,24.38,24,23.5,23.5,13.982,13.1,13.91,14.04,8.7])
flux = array([8.7,
14.04,
13.91,
13.1,
13.982,
23.5,
23.5,
24,
24.38,
43.5,
44.55,
45.05,
43.4,
40.8,
65.44,
83.58,
78.92,
97.09,
114.99,
132,
132,
245.8,
243])
# go to temperature GBT reported as 2.0K/Jy rather flat for low freq, dont care about high (will use gbtidl model)
ep = 390e-6
c = 299792458.0
na = 0.71*exp(-(4*pi*x*1e6/c)**2)
to = 0.008 + exp(sqrt(x/1e3))/8000.0
nl = 0.99
KpJy = 2.85 * na * nl *exp(-to)

y = flux

#mask = (x > 400.0) & (x < 1500.0)
#x = x[mask]
#y = y[mask]
#error = array([1.0,1.0,1.0,1.0,1.0,1.0,1.0])
#error = ones(len(x))
error = array([4.35E-01,
4.30E-01,
4.20E-01,
1.00E+00,
9.90E-02,
7.10E-01,
2.00E+00,
1.20E+00,
1.10E-01,
1.00E+00,
3.30E-01,
4.20E-01,
4.00E-01,
1.30E+00,
3.30E-01,
5.70E-01,
5.40E-01,
6.90E-01,
5.50E-01,
13,
1.30E+01,
3.19E+01,
6.00E+01])

print x.shape, y.shape, error.shape

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

xfit = arange(100,10000,100)

yfit = peval(xfit, fitp)

pylab.errorbar(x,y, yerr=error)
pylab.plot(xfit, yfit, 'r-')
pylab.show()
