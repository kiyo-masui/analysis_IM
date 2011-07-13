import numpy as np
import corr21cm

c = corr21cm.Corr21cm.from_file_matterps(fname="data/corr0.dat")

rv = np.linspace(0.0, 200.0, 1000)

xv = c.correlation_flatsky(rv, 0.0, z1 = 2.0)

