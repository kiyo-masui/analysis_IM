"""Module contains class that represents a calibration source."""

import numpy as np
import ephem

from utils import misc

catalogue = {
    '3C147': {
         'RA': "05:42:36.155",
         'DEC': "49:51:07.28"
         },
    '3C295': {
        'RA': "14:11:20.6",
        'DEC': "52:12:21"
        },
    '3C286': {
        'RA': "13:31:08.288",
        'DEC': "30:30:32.960"
        },
    '3C48': {
        'RA': "01:37:41.3",
        'DEC': "33:09:35."
        },
    }


class Source(object):
    """Represents a calibration source.

    Parameters
    ----------
    source_name: string
        Name of source catalogue source.
    """

    def __init__(self, source_name):
        # Check if the source is in the catalogue.
        if not source_name in catalogue.keys():
            raise ValueError("Source not in catalogue see "
                             "(`source.catalogue`).")
        self.source_name = source_name
        cat_info = catalogue[source_name]
        # Convert souce location to degrees.
        Loc = ephem.Equatorial(cat_info['RA'], cat_info['DEC'])
        RA, DEC = Loc.get()
        self.RA = RA * 180 / np.pi
        self.DEC = DEC * 180 / np.pi

    def azelGBT(self, UT):
        if isinstance(UT, np.ndarray):
            n = len(UT)
            az = np.empty(n, dtype=np.float64)
            el = np.empty(n, dtype=np.float64)
            for ii in range(n):
                az[ii], el[ii] = misc.radec2azelGBT(self.RA, self.DEC, UT[ii])
        else:
            az, el = misc.radec2azelGBT(self.RA, self.DEC, UT)
        return az, el


