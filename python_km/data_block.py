"""This module contains the class that holds an IF and scan of GBT data"""

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce

class DataBlock() :
    
    # These are the valid axes that a data field can vary over.
    axes = ('time', 'pol', 'cal', 'freq')

    def __init__(self, data=None) :
        """Can either be initialized with a raw data array (4D) or with None"""
        
        # Dictionary that holds all data other than .data.  This is safe to 
        # be accessed and updated by the user.
        self.field = {}
        # Dictionary with the same keys as field but holds the axes over which
        # a parameter varies.  For instance, the LST variable varies over the
        # 'time' axis.  axes['LST'] should thus be ('time') and
        # shape(field['LST']) should be (ntimes, ).
        self.field_axes = {}

        if data is None :
            self.data = ma.zeros((0,0,0,0), float)
            self.data_set = False
        else :
            self.set_data(data)

    def set_data(self, data) :
        """Set the data to passed array."""
        self.data = ma.array(data)
        self.data_set = True
        self.dims = sp.shape(data)

    def set_field(self, field_name, field_data, axis_names=()) :
        """
        Set field data to be stored.

        Note that these operation can also be done by accessing the 'field' and
        'field_axes' dictionaries directly, but using this function combines a
        few operations that go together.  It also does some sanity checks.
        Using this function is safer.
        """

        if type(axis_names) is str :
            a_names = (axis_names,)
        else :
            a_names = axis_names
        
        self._verify_single_axis_names(a_names)
        self.field[field_name] = field_data
        self.field_axes[field_name] = tuple(a_names)

    def _verify_single_axis_names(self, axis_names) :
        for name in axis_names :
            if not name in self.axes :
                raise ValueError("Field axes must contain only entries from: ",
                                 str(self.axes))
        # XXX: If someone decides they want to implement multi dimensional
        # fields, it shouldn't be to bad.  Make sure you update the fits
        # and writer though.
        if len(axis_names) > 1 :
            raise NotImplementedError("There is no reason we couldn't handle "
                                      "multi dimensional fields.")

    def verify(self) :
        """Verifies that all the data is consistant.

        This method should be run every time you muck around with the data
        and field entries.  It simply checks that all the data is consistant
        (axes, lengths etc.).

        Note that even if you know that your DataBlock will pass the verify,
        you still need to verify as this tells the DataBlock that you are done
        messing with the data.  It then sets some internal variables.
        """
        
        if self.data_set :
            self.dims = sp.shape(self.data)

        # Will delete these keys if they are found in 'field', then see if any
        # are left over.
        axes_keys = self.field_axes.keys()
        for field_name in self.field.iterkeys() :
            # Check for keys in fields and not in field_axes, the oposite is
            # done outside this loop.
            if not self.field_axes.has_key(field_name) :
                raise ce.DataError("Dictionaries 'field' and 'field_axes'"
                                   " must have the same keys.")
            axes_keys.remove(field_name)
            # Check all the axes
            axes = self.field_axes[field_name] # for saving keystrokes only
            self._verify_single_axis_names(axes)
            # Check the shape.
            field_data_shape = sp.shape(self.field[field_name])
            for ii in range(len(axes)) :
                axis_ind = list(self.axes).index(axes[ii])
                if field_data_shape[ii] != self.dims[axis_ind] :
                    raise ce.DataError("The shape of the data in one of the "
                                       "fields is incompatible with the shape "
                                       "of the main data. field: "+field_name)
        # The opposite of the first check in the loop.
        if len(axes_keys) :
            raise ce.DataError("Dictionaries 'field' and 'field_axes'"
                               " must have the same keys.")


        

