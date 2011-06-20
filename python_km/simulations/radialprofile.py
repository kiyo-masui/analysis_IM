import numpy as np

import pdb

def azimuthalAverage(image, center=None, bw = 3):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fractional pixels).
    
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    # Get sorted radii
    maxr = np.array([center[0], image.shape[0] - center[0], center[1], image.shape[1] - center[1]]).min()

    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    maxind = np.searchsorted(r_sorted, maxr+0.00001)
    r_sorted = r_sorted[:maxind]
    i_sorted = i_sorted[:maxind]

    numbins = int(maxr / bw)
    maxr = numbins * bw

    bn = np.linspace(0.0, maxr, numbins+1)

    bc = 0.5*(bn[1:] + bn[:-1])
    bl = bn[:-1]

    # Get the integer part of the radii (bin size = 1)
    #r_int = r_sorted.astype(int)
    r_int = np.digitize(r_sorted, bn)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    #pdb.set_trace()
    rind = np.where(deltar)[0]       # location of changed radius
    rind = np.insert(rind, 0, -1)
    nr = rind[1:] - rind[:-1]        # number of radius bin


    
    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    csim = np.insert(csim, 0, 0.0)
    tbin = csim[rind[1:]+1] - csim[rind[:-1]+1]

    radial_prof = tbin / nr

    return bl, radial_prof
