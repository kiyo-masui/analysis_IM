# save_2d.py
# This script saves the average 2d power and corresponding standard
# deviation across the six cross powers. This data is to be used in 
# the 2d Bayesian analysis. The data is saved in avg_result.hd5, 
# in the directory of the first cross power.

import h5py 
import core.algebra as al
import numpy as np
from mkpower import ps_summary as pss
from mkpower import ps_analysis as psa

# number of cross powers
map_num = 4
# run number of first cross power
first_cross = 446
#foreground modes
mode = 10
#fields
fields = ['ra165', 'ra199', 'ran18', 'ra33']

def cross_pow(map_num, first_run, mode):
    root = '/scratch2/p/pen/nluciw/parkes/analysis_IM/full%s_cros_ps_%smode/'
    pows = []
    stds = []

    for i in np.arange(first_run, first_run+map_num):
        print 'Attempting to open', root%(i,mode)
        ps_file = h5py.File(root%(i,mode) + 'ps_result.hd5', 'r')

        power_spectrum, k_num, k_p_edges, k_v_edges =\
            pss.load_power_spectrum('cros_ps_%dmode_2dpow'%mode, ps_file)
        transfer_function = \
            pss.load_transfer_function('cros_rf_%dmode_2dpow'%mode, 
                'cros_tr_%dmode_2dpow'%mode, ps_file, False)[0]
        sim_pow = ps_file['cros_si_%dmode_2draw'%mode].value *\
                       transfer_function[None,...]
    
        sim_list = []
        for j in range(sim_pow.shape[0]):
            each = al.make_vect(sim_pow[j], 
                   axis_names=power_spectrum.info['axes'])
            sim_list.append(each)
        std = np.std(np.array(sim_list), axis=0)

        stds.append(std)
        pows.append(power_spectrum)

    avg_file = h5py.File(root%(first_run,mode) + 'avg_result.hd5', 'a')

    stds = np.array(stds)
    pows = np.array(pows)
    stds[stds == 0] = np.inf

    pow = np.sum(pows*stds**(-2), axis=0)
    pow = div0(pow, np.sum(stds**(-2), axis=0))
    pow = al.make_vect(pow, axis_names=power_spectrum.info['axes'])
    pow.info = power_spectrum.info
    al.save_h5(avg_file, 'cros_ps_%dmode_2davg'%mode, pow)

    std = np.sqrt(div0(1., np.sum(np.array(stds)**(-2), axis=0)))
    std = al.make_vect(std, axis_names=power_spectrum.info['axes'])
    std.info = power_spectrum.info
    al.save_h5(avg_file, 'cros_std_%smode_2davg'%mode, std)    

def auto_pow(fields):
    root = '/scratch2/p/pen/nluciw/parkes/analysis_IM/field_%s_transfer_select/'
    pows = []
    stds = []

    for field in fields:
        ps_file = h5py.File(root%field + 'ps_result.hd5', 'r')

        power_spec, knum, k_p_edges, k_v_edges =\
            pss.load_power_spectrum('ps_2d', ps_file)
        std, knum, k_p_edges, k_v_edges =\
            pss.load_power_spectrum('std_2d', ps_file)

        pows.append(power_spec)
        stds.append(std)

    pows = np.array(pows)
    stds = np.array(stds)

    avg_file = h5py.File(root%fields[0] + 'avg_result.hd5', 'a')

    stds[stds == 0] = np.inf

    avg_pow = np.sum(pows * stds**(-2), axis=0)
    avg_pow = div0(avg_pow,  np.sum(stds**(-2), axis=0))
    avg_pow = al.make_vect(avg_pow, axis_names=power_spec.info['axes'])
    avg_pow.info = power_spec.info
    al.save_h5(avg_file, 'ps_2d', avg_pow)

    avg_std = np.sqrt(div0(1., np.sum(stds**(-2), axis=0)))
    avg_std = al.make_vect(avg_std, axis_names=power_spec.info['axes'])
    avg_std.info = power_spec.info
    al.save_h5(avg_file, 'std_2d', avg_std)

    psa.plot_2d_power_spectrum(avg_pow, None, k_p_edges, k_v_edges,
            filename='all_fields_2dpow', label_list='', output=root%fields[0])
    psa.plot_2d_power_spectrum(avg_std, None, k_p_edges, k_v_edges, 
            filename='all_fields_2dstd', label_list='', output=root%fields[0])

def div0(a, b):
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        c[~np.isfinite(c)] = 0
    return c

if __name__=='__main__':

    cross_pow(map_num, first_cross, mode)
#    auto_pow(fields)
