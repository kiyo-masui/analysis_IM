"""Functions usefull for all map makers.
"""

import scipy as sp
import numpy.ma as ma

from core import data_map

# Fields that should be copied from times stream data to maps.
fields_to_copy = (
                  'BANDWID',
                  'OBJECT'
                  )

def calc_inds(pointing, centre, shape, spacing=1.) :
    """Calculates a 1D map index corresponding to a location.

    Given a position, this funciton returns the 1D index on map of that
    position.  The map is specified by its centre, number of pixels across and
    the pixel width.

    Function is vectorized for input pointing.
    
    This fuction might be sufficiently general that it could be moved to
    core.utils.  We will see.
    """
    
    if type(pointing) is list or type(pointing) is tuple :
        pointing = sp.array(pointing, dtype=float)
    shape = int(shape)

    inds = sp.array((pointing - centre) / spacing + shape/2.0, dtype=int)
    
    return inds

def calc_bins(centre, shape, spacing=1., edge='lower') :
    """Calculates bin edges.

    Calculates a 1D array of regularly spaced bins.  User specifies whethar
    lower, middle, or upper side bins are to be calculated.
    """
    
    if edge == 'left' :
        offset = 0.
    elif edge == 'middle' :
        offset = spacing/2.0
    elif edge == 'right' :
        offset = float(spacing)
    else :
        raise ValueError("Argument edge must be either 'left', 'middle' or "
                         "'right'")
    shape = int(shape)
    
    bins = (spacing * sp.arange(-shape/2.0, shape/2.0, dtype=float)
            + centre + offset)
    
    return bins


def get_map_params(Map) :
    """Calculates the centre, shape, and pixel spacing of a map."""
    
    Map.calc_axes()
    shape = (len(Map.long), len(Map.lat), len(Map.freq))
    centre = ((Map.long[-1] + Map.long[0])/2.0, (Map.lat[-1] + Map.lat[0])/2.0,
              (Map.freq[-1] + Map.freq[0])/2.0)
    spacing = (Map.field['CDELT1'].item(), Map.field['CDELT2'].item(), 
               Map.field['CDELT3'].item())

    return centre, shape, spacing



def set_up_map(Data, centre, shape, spacing) :
    """Sets up a data_map.DataMap object.

    Initializes a data map by coping some information from a time stream data,
    which the angular information passed explicitly.

    Args:
        Data: A data_block.DataBlock object.  Time stream data from which to
        copy frequency information.
        centre, shape, spaceing:  All tuples with 2 elements each. Map centre
        (RA, DEC), map dimensions (pixels on the side) and pizel spacing
        (degrees).
    """
    
    if not (hasattr(centre, '__iter__') and hasattr(shape, '__iter__')
            and hasattr(spacing, '__iter__')) :
        raise ValueError('Arguments centre, shape, and spacing must all be'
                         ' length 2 iterables.')
    elif not (len(centre) == 2 and len(shape) == 2 and len(spacing) == 2) :
        raise ValueError('Arguments centre, shape, and spacing must all be'
                         ' length 2 iterables.')
    
    nf = Data.dims[-1]
    Map = data_map.DataMap(ma.zeros(tuple(shape) + (nf,), dtype=float))
    
    # Copy some data over that all the data_blocks should have
    # in common.
    # Do not copy RA and DEC to the map, since they cannot be
    # stored in the fits file.
    #Map.set_field('RA', ra_bins, ('long',), '1E') 
    #Map.set_field('DEC', dec_bins, ('lat',), '1E')
    Map.history = Data.history
    for key in fields_to_copy :
        Map.set_field(key, Data.field[key], 
           Data.field_axes[key], Data.field_formats[key])
    Map.set_field('CTYPE3', 'FREQ--HZ', (), '32A')
    Map.set_field('CTYPE1', 'RA---DEG', (), '32A')
    Map.set_field('CTYPE2', 'DEC--DEG', (), '32A')
    # Copy frequency axis (now the third axis not the first).
    Map.set_field('CRVAL3', Data.field['CRVAL1'], (), 'D')
    Map.set_field('CRPIX3', Data.field['CRPIX1'], (), 'D')
    Map.set_field('CDELT3', Data.field['CDELT1'], (), 'D')
    # Set the other two axes.
    ra_bins = calc_bins(centre[0], shape[0], spacing[0], 'middle')
    dec_bins = calc_bins(centre[1], shape[1], spacing[1], 'middle')
    Map.set_field('CRPIX1', shape[0]//2 + 1, (), 'D')
    Map.set_field('CRVAL1', ra_bins[shape[0]//2], (), 'D')
    Map.set_field('CDELT1', spacing[0], (), 'D')
    Map.set_field('CRPIX2', shape[1]//2 + 1, (), 'D')
    Map.set_field('CRVAL2', dec_bins[shape[1]//2], (), 'D')
    Map.set_field('CDELT2', spacing[1], (), 'D')
    
    Map.verify()

    return Map

def calc_time_var_file(Blocks, pol_ind=0, cal_ind=0) :
    """Calculates this time variance over several data blocks.

    Given a tuple of DataBlock objects (all with compatible dimensions), the
    time varience is calculated at each frequency.  pol_ind and cal_ind specify
    the polarization and cal.
    """

    # These all become arrays on first iteration.
    var = 0.0
    mean = 0.0
    counts = 0
    for Data in Blocks :
        var += ma.sum(Data.data[:,pol_ind,cal_ind,:]**2,0)
        mean += ma.sum(Data.data[:,pol_ind,cal_ind,:],0)
        counts += (Data.dims[0] 
                   - ma.count_masked(Data.data[:,pol_ind,cal_ind,:], 0))
    var = var/counts - (mean/counts)**2
    var[counts <= 1] = ma.masked
    var[var <= 0.0] = ma.masked

    return var
    
