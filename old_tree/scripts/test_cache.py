import filetools as ft
from core import algebra
import os
from utils import batch_handler as bh
from map import physical_gridding as pg

def get_cached_physical(filename, refinement=2, pad=5, order=1):
    basename = filename.split(".")[0]
    phys_cachename = basename + "_physical.npy"
    chksum_cachename = basename + ".md5"
    print phys_cachename, chksum_cachename

    curr_chksum = ft.hashfile(filename)
    # ALSO CHECK IF THE PARAMS CHANGED!

    # try to get an existing checksum
    try:
        chkfile = open(chksum_cachename, "r")
        old_chksum = chkfile.read()
        chkfile.close()
        if old_chksum == curr_chksum:
            chksum_not_changed = True
    except IOError as e:
        chksum_not_changed = False

    if os.path.isfile(phys_cachename) and chksum_not_changed:
        print "using the cached file: " + phys_cachename
        ret_data = algebra.make_vect(algebra.load(chksum_cachename))
    else:
        print "writing a physical cache for: " + filename
        # calculate the physical coordinate box
        obs_map = algebra.make_vect(algebra.load(filename))
        ret_data = bh.repackage_kiyo(pg.physical_grid(obs_map,
                                           refinement=refinement,
                                           pad=pad, order=order))
        algebra.save(chksum_cachename, ret_data)

        # save the new checksum
        chkfile = open(chksum_cachename, "w")
        chkfile.write(curr_chksum)
        chkfile.close()

    return ret_data


get_cached_physical("test.npy")
