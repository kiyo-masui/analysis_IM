from simulations import corr21cm

from utils import cubicspline as cs

## Modify with path to wigglez power spectrum
## Should be a file with columns of k-vals and ps-values.

## By default this will use the first two as k and ps, if not, try
## using the colspec keyword argument
ps_wigglez = cs.LogInterpolater.fromfile("<wigglez_ps_file>")


## Put an exponential cutoff into the powerspectrum at kstar in order
## to regularise small scales for integrals (eventually I'll fold this
## into the Corr21cm class).
kstar = 10.0  # units: h Mpc^{-1}
ps = lambda k: np.exp(-0.5 * (k / kstar)**2) * ps_wigglez(k)


## Set the redshift that the powerspectrum is normalised about.
z_ps = <wigglez ps effective redshift>


## Create the 21cm correlations object
cr = corr21cm.Corr21cm(ps = ps, redshift = z_ps)
