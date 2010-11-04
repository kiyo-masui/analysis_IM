#!/usr/bin/python
"""Puts time stream data in units of frequency.
"""

import base_single

class Calibrate(base_single.BaseSingle) :
    """Pipeline module converts data from units of cal temperture to Kelvins.

    This module reads the calibrator temperture from a fits file (as a function
    of polarization and frequency) and multiplies it by the time stream data.
    If the time stream data was in units of calibrator temperture, it will end
    up in units of acctual temperature (K).
    """
    
    prefix = 'cl_'
    params_init = {
                   #TODO: Have a real default (like a test file).
                   'cal_temperture_files' : ('some_file_name.fits',)
                   }



# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Calibrate(str(sys.argv[1])).execute()
