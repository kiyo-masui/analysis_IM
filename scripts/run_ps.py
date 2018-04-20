import numpy as np

from core import algebra

import os

from mkpower import ps_estimator2 as ps

from mkpower import ps_summary2 as psum

pobject = ps.PowerSpectrumEstimator()

map_root = os.getenv('MAPROOT')

field = os.getenv('FIELD')

path21 = map_root + field + os.getenv('MAP')
#path21 = '/scratch2/p/pen/bh5/analysis_IM-mkpower/21cm_array.npy'
path22 = path21

#path22 = '/scratch2/p/pen/bh5/analysis_IM-mkpower/21cm_array.npy'

noise_path = map_root + field + os.getenv('NOISE')

#type = os.getenv('TYPE')

#output = '/scratch2/p/pen/andersoc/2df_data/power_spectra/'
output = os.getenv('OUTPUT') + field

#nfreq = 64

kbin_min = 0.02
kbin_max = 10
#kbin_max = 0.2
kbin_num = 30

pobject.params['kbin_min'] = 0.05

pobject.params['kbin_max'] = 2

pobject.params['kbin_num'] = 6

ps_est = pobject.ps_estimate(pobject.params, [[path21],[path22]], [[noise_path],[noise_path]], path21)

params = pobject.params

k_bin = pobject.get_kbin_centr(pobject.params)

k_axes_2d = ("k_p", "k_v")
info_2d = {'axes': k_axes_2d, 'type': 'vect'}
info_2d['k_p_delta']  = k_bin[1]/k_bin[0]
info_2d['k_p_centre'] = k_bin[params['kbin_num']//2]
info_2d['k_v_delta']  = k_bin[1]/k_bin[0]
info_2d['k_v_centre'] = k_bin[params['kbin_num']//2]

ps_2d = ps_est[0]
kn_2d = ps_est[1]

ps_2d  = algebra.make_vect(ps_2d, axis_names=k_axes_2d)
kn_2d  = algebra.make_vect(kn_2d, axis_names=k_axes_2d)
ps_2d.info = info_2d
kn_2d.info = info_2d

onedpwr = psum.convert_2dps_to_1dps_each([ps_2d],[kn_2d], None)

#print type(onedpwr[0])

np.ma.dump(onedpwr[0], output + '1d_power')

np.ma.dump(onedpwr[2], output + '1d_k')


