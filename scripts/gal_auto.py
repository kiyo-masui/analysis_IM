import h5py 
import numpy as np
from core import algebra as al
from mkpower import ps_analysis as psa
from mkpower import ps_summary as pss
import copy 

file_start = 486
fields = 4
type = '2df'

powers_2d = []
kmns = []

for ii in range(fields):

    file = h5py.File('full%s_%s_ps_10mode/ps_result.hd5'%(np.str(file_start+ii),type))
    power_2df, kmn, kpedge, kvedge = pss.load_power_spectrum('%s_ps_10mode_2dpow'%type,\
                                     file)

    powers_2d.append(power_2df)
    kmns.append(kmn)

std = (np.std(np.array(powers_2d), axis=0))
std[std==0] = np.inf
wt = std**(-2)

fullpow, err, k = pss.convert_2dps_to_1dps_each(copy.deepcopy(powers_2d), 
    copy.deepcopy(kmns), np.ones_like(powers_2d[0]), None)

gal_1d_list = []

for ii in range(len(powers_2d)):
    gal_1d, gal_1d_std, k_cent = pss.convert_2dps_to_1dps_each(
        powers_2d[ii], kmns[ii], np.ones_like(powers_2d[ii]), None)

    gal_1d_list.append(gal_1d)
    gal_1d[fullpow==0] = 0

    gal_pos, gal_neg = psa.seperate_positive_and_negative_power(
        gal_1d, err, k_cent)

    fname = 'full%s_%s_ps_10mode/comp_varweight_cros_ps_10mode_1dpow_%s.txt'

    np.savetxt(fname%(np.str(file_start+ii),type,'positive'), gal_pos.T, fmt='%15.12e')
    np.savetxt(fname%(np.str(file_start+ii),type,'negative'), gal_neg.T, fmt='%15.12e')

psa.plot_2d_power_spectrum(powers_2d[0], None, kpedge, kvedge, 
                           filename='%s_ps_2dpow'%type, label_list='', 
                           output='full%s_%s_ps_10mode'%(np.str(file_start),type))

std[std==np.inf] = 0
psa.plot_2d_power_spectrum(std, None, kpedge, kvedge, 
                           filename='std_ps_2dpow', label_list='', 
                           output='full%s_%s_ps_10mode'%(np.str(file_start),type))

fullpos, fullneg = psa.seperate_positive_and_negative_power(
    np.mean(gal_1d_list, axis=0), err, k_cent)

np.savetxt('redpos', fullpos)
np.savetxt('redneg', fullneg)





