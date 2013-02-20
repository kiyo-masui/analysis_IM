from utils import cosmology
from simulations import corr
from utils.cosmology import Cosmology

redshift = 1.

cosmo = Cosmology()
proper = cosmo.proper_distance(redshift)
comoving = cosmo.comoving_distance(redshift)
comoving_inv = cosmology.inverse_approx(cosmo.comoving_distance, 0.6, 1.)

print proper, comoving, comoving_inv(2355.35909781)

