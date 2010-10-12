"""Some general utilities that are useful."""

import time

import scipy as sp
import ephem

def elaz2radec_lst(el, az, lst, lat = 38.43312) :
    """DO NOT USE THIS ROUTINE FOR ANTHING THAT NEEDS TO BE RIGHT.  IT DOES NOT
    CORRECT FOR PRECESSION.
    
    Calculates the Ra and Dec from elavation, aximuth, LST and Latitude.
    
    This function is vectorized with numpy so should be fast.  Standart numpy
    broadcasting should also work.

    All angles in degrees, lst in seconds. Latitude defaults to GBT.
    """
    
    # Convert everything to radians.
    el = sp.radians(el)
    az = sp.radians(az)
    lst = sp.array(lst, dtype = float)*2*sp.pi/86400
    lat = sp.radians(lat)
    # Calculate dec.
    dec = sp.arcsin(sp.sin(el)*sp.sin(lat) + 
                    sp.cos(el)*sp.cos(lat)*sp.cos(az))
    # Calculate the hour angle
    ha = sp.arccos((sp.sin(el) - sp.sin(lat)*sp.sin(dec)) / 
                   (sp.cos(lat)*sp.cos(dec)))
    ra = sp.degrees(lst - ha) % 360

    return ra, sp.degrees(dec)

def elaz2radecGBT(el, az, UT) :
    """Calculates the Ra and Dec from the elevation, azimuth and UT for an
    observer at GBT.
    
    All input should be formated to correspond to the data in a GBT fits file.
    El, and Az in degrees and UT a string like in the GBT DATE-OBS field.

    Largely copied from Kevin's code.
    """

    GBT = ephem.Observer()
    GBT.long = '-79:50:23.4'
    GBT.lat = '38:25:59.23'
    GBT.pressure = 0 # no refraction correction.
    GBT.temp = 0

    UT_wholesec, partial_sec = UT.split('.', 1)
    time_obj = time.strptime(UT_wholesec, "%Y-%m-%dT%H:%M:%S")
    UT_reformated = time.strftime("%Y/%m/%d %H:%M:%S", time_obj)
    GBT.date = UT_reformated + "." + partial_sec
    
    el_r = el*sp.pi/180.0
    az_r = az*sp.pi/180.0
    ra, dec = GBT.radec_of(az_r,el_r)

    return ra*180.0/sp.pi, dec*180.0/sp.pi


def mk_map_grid(centre, shape, spacing) :
    """Make a grid of coordinates in Ra and Dec.

    This function accepts a field centre (tuple 2 floats), map-shape (tuple 2
    ints) and pixel spacing (float) for a map.  It returns two arrays with
    the specified map shape, and meshed like sp.meshgrid.  However, the spacing
    in Ra is converted to real degrees (devided by cos(dec) at field centre).

    All units in degrees, values are pixel centres.
    """

    dec = centre[1] + spacing*sp.arange(-(shape[1]-1.)/2., shape[1]/2.)
    ra = centre[0] + (spacing/sp.cos(centre[1]*sp.pi/180.) * 
                      sp.arange(-(shape[0]-1.)/2., shape[0]/2.))

    grid_ra, grid_dec = sp.meshgrid(ra, dec)

    return grid_ra, grid_dec



