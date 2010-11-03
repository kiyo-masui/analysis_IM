"""Data container for maps."""

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


# Clone some extra functions:

merge_histories = base_data.merge_histories

