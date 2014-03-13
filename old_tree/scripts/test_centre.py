import numpy as np
import math
from core import algebra

nsamp = 102
x_axis = np.linspace(0, 2. * math.pi, num = nsamp, endpoint=True)
y_axis = np.random.normal(size=(nsamp))*math.sqrt(nsamp/2./math.pi)

fftarr = np.fft.fftshift(np.fft.fftn(y_axis))
k_axes = ["x_axis"]
xspec_arr = algebra.make_vect(fftarr, axis_names=k_axes)

info = {'axes': k_axes, 'type': 'vect'}
delta_axis = abs(x_axis[1] - x_axis[0])

k_axis = np.fft.fftshift(np.fft.fftfreq(nsamp, d=delta_axis))
delta_k_axis = abs(k_axis[1] - k_axis[0])

info["x_axis_delta"] = delta_k_axis
info["x_axis_centre"] = 0.

xspec_arr.info = info

print k_axis-xspec_arr.get_axis("x_axis")
print xspec_arr.get_axis("x_axis")
