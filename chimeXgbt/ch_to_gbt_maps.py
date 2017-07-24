from core import algebra as al

def CHImap_to_GBTvect(map, axes, info):
    ''' Takes contents of a CHIME map container and converts it to a 
        vector from Kiyo's algebra module, for use with the GBT/Parkes
        pipeline. Not currently set up to handle axes without 
        associated 'cent' and 'delt' values.

        map:    numpy array containing the map
        axes:   list of axis names. len(axes) must equal map.ndim
                eg: ['freq', 'ra']
        info:   dictionary with metadata for each axis. len(info) must
                equal 2*len(axes)
                eg: {'freq_cent': 10., 'freq_delt': 0.5, 'ra_cent':...}

        Returns a GBT vector object.
    '''

    # make sure info and axes are compatible
    if 2*len(axes) != len(info):
        print 'Error: each axis should have associated centre and delta info'
        return -1

    # create algebra vector
    vect = al.make_vect(map, axis_names=axes)

    # set metadata
    for axis in axes:
        vect.set_axis_info(axis, info[axis + '_cent'], info[axis + '_delt'])

    return vect
