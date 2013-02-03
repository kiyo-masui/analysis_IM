"""
Important info about structure of the binary file

This file is a binary that has a header of one float(4) which you can discard,
then
it has as many lines as halos, and 28 columns:
1:3 are position of the peak
4:6 are position of the centre of mass (CM)
7:9 are velocity of CM
10:12 are angular momentum in CM frame
13:15 are velocity dispersion about the CM
16 is the radius
17:19 are different estimates of the halo mass. 17 is the most accurate.
20:22 are some kind of variance about the the position of the CM,
23:28 are the 6 independent elements of the Inertia matrix.

What you need to know is that the positions are in grid cells units (2048^3
cells total),
and the file I gave you is a sub-volume (1-64th) of a [505 Mpc/h]^3 simulation.
So the file probably will have posits running from 0.0 to 512.0.

The conversion is 2048 cells = 505 Mpc/h.

The mass is in simulation unit as well, and for this particular simulation, the 
conversion unit is:

1 mass unit = 1.3313e10 Solar mass.
"""
# run code : "!python this_file.py > data_file.dat" to produce data_file

import struct

# binary file we want to convert
file_name = "/cita/h/home-2/mufma/Secondary/0.696halo0.dat"

# set sequence in which binary file would be read
record_header = struct.Struct("f")
halo_struct = struct.Struct("<" + "f"*28)

# scaling factors to physical units
physical_scale = 505. * (1./ 2048.)
mass_scaling = 1.3313e10
# velocity_scale = 

# reading binary : opened_file.read(move_cursor.size) - it moves cursor along
# the "opened_file" on "move_cursor.size" steps
# the while loop works until the variable binary_info becomes False
with open(file_name,"rb") as halofile:
    header = halofile.read(record_header.size)
    binary_info = halofile.read(halo_struct.size)
    while binary_info:
        halo_info = halo_struct.unpack(binary_info)
        print halo_info[0] * physical_scale, \
              halo_info[1] * physical_scale, \
              halo_info[2] * physical_scale, \
              halo_info[5] * velocity_scale, \
              halo_info[6] * velocity_scale, \
              halo_info[7] * velocity_scale, \
              halo_info[16] * mass_scaling
        binary_info = halofile.read(halo_struct.size)
