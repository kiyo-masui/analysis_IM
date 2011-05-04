
import foregroundmap
import numpy as np

import lofar

tf = foregroundmap.ForegroundSCK()

#f1 = tf.getfield()

freq = np.linspace(tf.nu_l, tf.nu_h, tf.numf)

cv = tf.frequency_covariance(*np.meshgrid(freq, freq))
#nm, aff = tf.getfield()

f1 = tf.getfield()

fdt = f1 * ((freq/tf.nu_0)**2.80)[:,np.newaxis,np.newaxis]

lf = lofar.LofarGSDE()

f2 = lf.getfield()
