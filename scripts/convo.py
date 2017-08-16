import map.beam as beam
import numpy as np
import scipy as sp

streaming_dispersion = 500.*0.72 / np.sqrt(2.)
freq_data = sp.array([1250, 1275, 1300, 1325, 1350, 1430], dtype=float)
# beam_data in unit of dgree
beam_data = sp.array([14.4, 14.4, 14.4, 14.4, 14.4, 14.4])/60.
beam_data = beam_data*1420/freq_data
freq_data *= 1.0e6
beam_data = sp.sqrt(beam_data**2)
degrade = 1.4
beam_diff = np.sqrt(max(degrade * beam_data) ** 2 - (beam_data) ** 2)

def convolve(map):
    reso = beam.GaussianBeam(beam_diff, freq_data)
    degrade = reso.apply(map)

    return degrade
