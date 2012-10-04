import numpy as np
import h5py

def rotate_pol_moon(Data, rotationfile):
    r"""Rotation the polarizations of a TOD in such a way as to make the
    observation of the moon pure-I
    `Data` is the TOD
    `rotationfile` is the hd5 file with the 4x4 rotation per frequency
    The rotation is in {XX, XY, YX, YY}
    """
    if tuple(Data.field['CRVAL4']) != (-5, -7, -8, -6) :
        raise ValueError

    rotation = h5py.File(rotationfile, "r")
    rotation = rotation['moon_scan1'].value
    print rotation

    # could be faster with tensordot?
    nfreq = Data.dims[3]
    # data are time_index,polindex,cal_index,freq
    Data.data = np.rollaxis(Data.data, 1, 4)

    for freq in range(0, nfreq):
        Data.data[:, :, freq, :] = np.dot(Data.data[:, :, freq, :],
                                          rotation[:, :, freq])

    Data.data = np.rollaxis(Data.data, 3, 1)

