
import numpy as np

from simulations import foregroundsck, lofar

tf = foregroundsck.Synchrotron()

#f1 = tf.getfield()

freq = tf.nu_pixels

cv = tf.frequency_covariance(*np.meshgrid(freq, freq))
#nm, aff = tf.getfield()

f1 = tf.getfield()

fdt = f1 * ((freq/tf.nu_0)**2.80)[:,np.newaxis,np.newaxis]

lf = lofar.LofarGDSE()

f2 = lf.getfield()
