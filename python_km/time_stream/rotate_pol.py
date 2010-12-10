#!/usr/bin/python
"""Module that converts between representations of polarization."""

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
import base_single


class RotatePol(base_single.BaseSingle) :
    """Pipeline module that converts correlation power to stoke parameters.
    
    See the main function of this module: rotate_pol.rotate for a detailed
    doc string.
    """

    prefix = 'rp_'
    params_init = {
                   # The polarizations that should be included in the output
                   # data. 1, 2, 3, 4 are I, Q, U ,V respectively (SD fits
                   # convension).
                   'new_pols' : (1,),
                   'average_cals' : False
                   }

    def action(self, Data):
        rotate(Data, self.params['new_pols'])
        Data.add_history('Rotated polarizations parameters.', ('Rotated to:' +
                                                str(self.params['new_pols']),))
        return Data


def rotate(Data, new_pols=(1,), average_cals=False) :
    """Changes the basis of the polarization axis.

    Passed a data_block.DataBlock object and the new polarization axis.
    Polarizations follow the SDfits convensions: 1=I, 2=Q, 3=U, 4=V, -5=XX,
    -6=YY, -7=XY, -8=YX.

    Right now this can only convert XX, etc to I, but eventually it should be
    expanded to go from any complete basis to any other set of polarizations.

    It also optionally takes the average of cal_on and cal_off.
    """
    
    # Here we check the polarizations indicies
    xx_ind = 0
    yy_ind = 3
    xy_inds = [1,2]
    if (Data.field['CRVAL4'][xx_ind] != -5 or
        Data.field['CRVAL4'][yy_ind] != -6 or
        Data.field['CRVAL4'][xy_inds[0]] != -7 or
        Data.field['CRVAL4'][xy_inds[1]] != -8) :
            raise ce.DataError('Polarization types not as expected.')
    on_ind = 0
    off_ind = 1
    if (Data.field['CAL'][on_ind] != 'T' or
        Data.field['CAL'][off_ind] != 'F') :
            raise ce.DataError('Cal states not in expected order.')
    if len(new_pols) == 1 and new_pols[0] == 1 :
        I = (Data.data[:,[xx_ind],:,:] + Data.data[:,[yy_ind],:,:])/2.0
        if average_cals :
            I = (I[:,:,[0],:] + I[:,:,[1],:])/2.0
        Data.set_data(I)
        Data.field['CRVAL4'] = sp.array([1])
    elif tuple(new_pols) == (1,2,3,4) :
        new_data = ma.empty(Data.dims)
        new_data[:,[0],:,:] = (Data.data[:,[xx_ind],:,:] + 
                             Data.data[:,[yy_ind],:,:])/2.0
        new_data[:,[1],:,:] = (Data.data[:,[xx_ind],:,:] - 
                             Data.data[:,[yy_ind],:,:])/2.0
        new_data[:,[2],:,:] = Data.data[:,[xy_inds[0]],:,:] 
        new_data[:,[3],:,:] = Data.data[:,[xy_inds[1]],:,:] 
        if average_cals :
            new_data = (new_data[:,:,[0],:] + new_data[:,:,[1],:])/2.0
        Data.set_data(new_data)
        Data.field['CRVAL4'] = sp.array([1,2,3,4])
    else :
        raise NotImplementedError('Right now we can only calculate I'
                                  ' or (I, Q, U V)')
    

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    RotatePol(str(sys.argv[1])).execute()

