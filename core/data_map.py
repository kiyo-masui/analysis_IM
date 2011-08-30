"""Data container for maps."""

import scipy as sp

import base_data
import kiyopy.custom_exceptions as ce

class DataMap(base_data.BaseData) :
    """This class holds data from a map."""

    axes = ('long', 'lat', 'freq')
    
    # Add some additional checks to verify().
    def verify(self) :
        """Mostly the same as BaseData.verify, but with added checks.
        """

        for field_axes in self.field_axes.itervalues() :
            # The fields in the map must be stored in the header (not the data)
            # and as such, they must be 0 dimensional (scalars).
            if len(field_axes) != 0 :
                raise ce.DataError("Maps can only have scalar fields.")
        # Now call verify from the base class.
        base_data.BaseData.verify(self)

    def calc_axes(self) :
        """Calculates that frequency, lat and long axes.

        These are not stored as fields as all fields have to be writable to the
        fits file.
        """
        if (self.field.has_key('CRPIX3') and self.field.has_key('CDELT3') and 
            self.field.has_key('CRVAL3')) :
            self.freq = ((sp.arange(self.dims[2], dtype=float) + 1.0 - 
                         self.field['CRPIX3'])*self.field['CDELT3'] + 
                         self.field['CRVAL3'])
        if (self.field.has_key('CRPIX2') and self.field.has_key('CDELT2') and 
            self.field.has_key('CRVAL2')) :
            self.lat = ((sp.arange(self.dims[1], dtype=float) + 1.0 - 
                         self.field['CRPIX2'])*self.field['CDELT2'] + 
                         self.field['CRVAL2'])
        if (self.field.has_key('CRPIX1') and self.field.has_key('CDELT1') and 
            self.field.has_key('CRVAL1')) :
            self.long = ((sp.arange(self.dims[0], dtype=float) + 1.0 - 
                         self.field['CRPIX1'])*self.field['CDELT1'] + 
                         self.field['CRVAL1'])

# Clone some extra functions:
from hist import merge_histories
