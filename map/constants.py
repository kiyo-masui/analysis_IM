"""Some constants used by the map makers.

These constants are defined here so other modules (namely
time_stream.measure_noise) don't need to import dirty map to have access to
them.
"""

# Absolute constants.

# Constant that represents a very high noise level.  Setting the noise of a
# mode to this number deweights that mode.
# 100 K**2 seems like a good number.  Bigger than any actual noise mode we would
# encounter, but not so big that it screws up numerical stability.
T_infinity = 100.0  # Kelvin**2
T_small = T_infinity * 1.e-10
# For deweighting by a lot, but lower risk of messing with numerical stability.
T_large = 1.0
T_medium = 0.01
T_sys = 625  # K**2


# Relative constants.

f_infinity = 1.e6
f_large = 1.e4
f_medium = 1.e2

