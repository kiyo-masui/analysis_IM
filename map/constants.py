"""Some constants used by the map makers.

This file defines a hierarchy of noise temperatures and are all based on
physical temperatures and noise levels.  This hierarchy is experiment dependant.

All temperatures in Kelvin.

These constants are defined here so other that modules (namely
time_stream.measure_noise) don't need to import dirty map to have access to
them.
"""

# Absolute constants.

# Constant that represents an infinite noise level.  Only use this when there
# is no risk of creating numerical instability.
T_infinity = 1000.0  # Kelvin

# Noise level much higher than any scale in the system but not so high that it
# guarantees numerical instability if used.
T_huge = 10

# The scale of the foregrounds in the map.  Still very high compared to the
# signal.
T_large = 1.0

T_medium = 0.01

T_small = 1.e-4

# The scale of the telescope system temperature. Approximate.
T_sys = 25.  # K


# Relative constants.

f_infinity = 1.e6
f_large = 1.e4
f_medium = 1.e2

