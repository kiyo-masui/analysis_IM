#!/usr/bin/python
"""Module that converts between representations of polarization."""

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
import base_single
#import combine_cal


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
        rotate(Data, self.params['new_pols'],
               average_cals=self.params['average_cals'])
        Data.add_history('Rotated polarizations parameters.', ('Rotated to:' +
                                                str(self.params['new_pols']),))
        if self.params['average_cals'] :
            Data.add_history('Averaged cal states.')
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
    
    # Two supported cases for input polarizations.
    if (tuple(Data.field['CRVAL4']) == (1, 2)) :
        if tuple(new_pols) == (1,):
            new_data = Data.data[:,:,[0],:,:] + Data.data[:,:,[1],:,:]
        else:
            msg = ("Converstion to " + str(tuple(new_pols)) + " from " 
                   + str(tuple(Data.field['CRVAL4'])) + " is not supported.")
            raise NotImplementedError(msg)
    elif (tuple(Data.field['CRVAL4']) == (1, 2, 3, 4)) :
        I_ind = 0
        Q_ind = 1
        U_ind = 2
        V_ind = 3
        # Two supported cases for output polarizations.
        if tuple(new_pols) == (1,) :
            new_data = Data.data[:, [0], :, :]
        elif tuple(new_pols) == (-5, -7, -8, -6) :
            new_data = ma.empty(Data.dims)
            new_data[:,[0],:,:] = (Data.data[:,[I_ind],:,:] 
                                   - Data.data[:,[Q_ind],:,:])
            new_data[:,[1],:,:] = Data.data[:,[U_ind],:,:] 
            new_data[:,[2],:,:] = Data.data[:,[V_ind],:,:] 
            new_data[:,[3],:,:] = (Data.data[:,[I_ind],:,:] 
                                   + Data.data[:,[Q_ind],:,:])
        else :
            msg = ("Converstion to " + str(tuple(new_pols)) + " from " 
                   + str(tuple(Data.field['CRVAL4'])) + " is not supported.")
            raise NotImplementedError(msg)
    elif (tuple(Data.field['CRVAL4']) == (-5, -7, -8, -6)):
        xx_ind = 0
        yy_ind = 3
        xy_ind = 1
        yx_ind = 2
        # Two supported cases for output polarizations.
        if tuple(new_pols) == (1,) :
            new_data = (Data.data[:,[xx_ind],:,:] + 
                        Data.data[:,[yy_ind],:,:])/2.0
        elif tuple(new_pols) == (1,2,3,4) :
            new_data = ma.empty(Data.dims)
            new_data[:,[0],:,:] = (Data.data[:,[xx_ind],:,:] + 
                                 Data.data[:,[yy_ind],:,:])/2.0
            new_data[:,[1],:,:] = (-Data.data[:,[xx_ind],:,:] + 
                                 Data.data[:,[yy_ind],:,:])/2.0
            new_data[:,[2],:,:] = Data.data[:,[xy_ind],:,:] 
            new_data[:,[3],:,:] = Data.data[:,[yx_ind],:,:] 
        else :
            msg = ("Converstion to " + str(tuple(new_pols)) + " from " 
                   + str(tuple(Data.field['CRVAL4'])) + " is not supported.")
            raise NotImplementedError(msg)
    else :
        raise NotImplementedError('For now polarizations must be (I, Q, U, V)'
                                  ' or (XX, XY, YX, YY) in those orders')
    print new_data.shape
    Data.set_data(new_data)
    Data.field['CRVAL4'] = sp.array(new_pols)
    Data.field['TSYS'] = Data.field['TSYS'][:,:,[0],:]+Data.field['TSYS'][:,:,[1],:]

    ## Now deal with the cal if desired.
    #if average_cals :
    #    combine_cal.combine(Data, weights=(0.5, 0.5), sub_mean=False,
    #                        average_cals=True)


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    RotatePol(str(sys.argv[1])).execute()

