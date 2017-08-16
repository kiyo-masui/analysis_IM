import h5py
import core.algebra as al
import numpy as np
from mkpower import ps_summary as pss
import pylab as pl
from mkpower import ps_analysis as psa
import os 

select_type = 'sn'

modes = [0]
runStart = 228

runs = [x for x in np.arange(runStart, runStart+len(modes))]

inputroot = '/scratch2/p/pen/nluciw/parkes/analysis_IM/'
outputroot = '/scratch2/p/pen/nluciw/parkes/analysis_IM/%s_thru_%s_%s_select/'\
             %(runs[0],runs[-1],select_type)
if not os.path.exists(outputroot):
    os.makedirs(outputroot)

rf_root = 'cros_rf_%dmode_2dpow'
tr_root = 'cros_tr_%dmode_2dpow'
si_root = 'cros_si_%dmode_2draw'

sims = []
stds = []
si_2df = []
pow_2df = []
pow = []
transfer = []
sn = []

for i in range(len(modes)):
    ps_root = inputroot + 'full%s_cros_ps_%02dmode/'%(runs[i],modes[i]) + 'ps_result.hd5'
    ps_result = h5py.File(ps_root, 'r')

    power_spectrum_sim, k_mode_number_sim, sim_p_edges, sim_v_edges =\
        pss.load_power_spectrum('cros_si_%dmode_2dpow'%modes[i], ps_result)

    transfer_function = \
            pss.load_transfer_function(rf_root%modes[i], tr_root%modes[i], ps_result, \
                    False)[0]
    transfer.append(transfer_function)

    pow_2d_raw_sim = ps_result[si_root%modes[i]].value
    power_2d_raw_kmn_sim = ps_result['cros_si_%dmode_2draw_kmn'%modes[i]].value

    power_2d_raw_sim = pow_2d_raw_sim * transfer_function[None,...]

    sim_list = []
    sim_kmn_list = []
    for j in range(power_2d_raw_sim.shape[0]):
        sim_each = al.make_vect(power_2d_raw_sim[j],
                axis_names=power_spectrum_sim.info['axes'])
        sim_each.info = power_spectrum_sim.info
        kmn_each = al.make_vect(power_2d_raw_kmn_sim[j],
                axis_names=power_spectrum_sim.info['axes'])
        kmn_each.info = power_spectrum_sim.info
        sim_list.append(sim_each)
        sim_kmn_list.append(kmn_each)
    std = np.std(np.array(sim_list), axis=0)
    
    sims.append(sim_list)
    stds.append(std)

    sim_spectrum_2df, k_mode_number_si2df, k_p_edges, k_v_edges =\
            pss.load_power_spectrum('2df_si_%dmode_2dpow'%modes[i], ps_result)
    si_2df.append(sim_spectrum_2df)    

    #short_noise, short_noise_kmn, k_p_edges, k_v_edges =\
    #        pss.load_power_spectrum('2df_sn_%dmode_2dpow'%modes[i], ps_result)
    power_spectrum_2df, k_mode_number_2df, k_p_edges, k_v_edges =\
            pss.load_power_spectrum('2df_ps_%dmode_2dpow'%modes[i], ps_result)
    #power_spectrum_2df -= short_noise
    pow_2df.append(power_spectrum_2df)

    cros_short_noise, short_noise_kmn, k_p_edges, k_v_edges =\
            pss.load_power_spectrum('cros_sn_%dmode_2dpow'%modes[i], ps_result)
    power_spectrum, k_mode_number, k_p_edges, k_v_edges =\
            pss.load_power_spectrum('cros_ps_%dmode_2dpow'%modes[i], ps_result)
    power_spectrum -= cros_short_noise
    power_spectrum *= transfer_function
    pow.append(power_spectrum)

    if modes[i] == 10:
        psa.plot_2d_power_spectrum(power_spectrum, None, k_p_edges, k_v_edges,
                               filename='comp_cros_ps_%dmode'%(modes[i]),
                               label_list='', output=outputroot)
                
sims = np.array(sims)
stds = np.array(stds)
si_2df = np.array(si_2df)
pow_2df = np.array(pow_2df)
pow = np.array(pow)
transfer = np.array(transfer)

sim_select = np.zeros_like(sims[0])
std_select = np.zeros_like(stds[0])
si_2df_select = np.zeros_like(si_2df[0])
pow_2df_select = np.zeros_like(pow_2df[0])
pow_select = np.zeros_like(pow[0])
transfer_select = np.zeros_like(transfer[0])
mode_map = np.zeros_like(transfer[0])

sn = np.zeros_like(pow)
stds[stds == 0] = np.inf
for i in range(len(modes)):
    sn[i] = np.mean(sims[i], axis=0) / stds[i]
stds[stds == np.inf] = 0

if select_type == 'sn':
    for i in range(std_select.shape[0]):
        for j in range(std_select.shape[1]):
            print 'Options:', sn[:,i,j]
            std_select[i,j] = stds[np.where(sn[:,i,j] == sn[:,i,j].max())[0][0],i,j]
            sim_select[:,i,j] = sims[np.where(sn[:,i,j] == sn[:,i,j].max())[0][0],:,i,j]
            #si_2df_select[i,j] = si_2df[np.where(stds[:,i,j] == stds[:,i,j].min())[0][0],i,j]
            #pow_2df_select[i,j] = pow_2df[np.where(stds[:,i,j] == stds[:,i,j].min())[0][0],i,j]
            pow_select[i,j] = pow[np.where(sn[:,i,j] == sn[:,i,j].max())[0][0],i,j]
            transfer_select[i,j] = transfer[np.where(sn[:,i,j] == sn[:,i,j].max())[0][0],i,j]
            mode_map[i,j] = modes[np.where(stds[:,i,j] == stds[:,i,j].max())[0][0]]

if select_type == 'variance':
    for i in range(std_select.shape[0]):
        for j in range(std_select.shape[1]):
            std_select[i,j] = stds[:,i,j].min()
            sim_select[:,i,j] = sims[np.where(stds[:,i,j] == stds[:,i,j].min())[0][0],:,i,j]
            si_2df_select[i,j] = si_2df[np.where(stds[:,i,j] == stds[:,i,j].min())[0][0],i,j]
            pow_2df_select[i,j] = pow_2df[np.where(stds[:,i,j] == stds[:,i,j].min())[0][0],i,j]
            pow_select[i,j] = pow[np.where(stds[:,i,j] == stds[:,i,j].min())[0][0],i,j]
            transfer_select[i,j] = transfer[np.where(stds[:,i,j] == stds[:,i,j].min())[0][0],i,j]
            mode_map[i,j] = np.where(stds[:,i,j] == stds[:,i,j].min())[0][0] * 5 + 5 

transfer[transfer == 0] = np.inf
transfer = 1. / transfer
if select_type == 'transfer':
    for i in range(transfer_select.shape[0]):
        for j in range(transfer_select.shape[1]):
            transfer_select[i,j] = transfer[np.argmin\
                (np.abs(0.8 - transfer[:,i,j])),i,j]
            std_select[i,j] = stds[np.argmin\
                (np.abs(0.8 - transfer[:,i,j])),i,j]
            sim_select[:,i,j] = sims[np.argmin\
                (np.abs(0.8 - transfer[:,i,j])),:,i,j]
            si_2df_select[i,j] = si_2df[np.argmin\
                (np.abs(0.8 - transfer[:,i,j])),i,j]
            pow_2df_select[i,j] = pow_2df[np.argmin\
                (np.abs(0.8 - transfer[:,i,j])),i,j]
            pow_select[i,j] = pow[np.argmin\
                (np.abs(0.8 - transfer[:,i,j])),i,j]
            mode_map[i][j] = np.argmin(np.abs(0.8 - transfer[:,i,j])) * 5 + 5

print mode_map
mode_map[pow_select == 0] = 0
psa.plot_2d_power_spectrum(mode_map, None, k_p_edges, k_v_edges, 
                       filename='mode_map', label_list='', output=outputroot,
                       logscale=False, cmax=20., cmin=5.)

psa.plot_2d_power_spectrum(np.mean(sim_select, axis=0), None, sim_p_edges, sim_v_edges,
                       filename='comp_cros_si_%d-%dmode'%(modes[0],modes[-1]),
                       label_list='', output=outputroot)

psa.plot_2d_power_spectrum(std_select, None, k_p_edges, k_v_edges,
                           filename='comp_cros_std_%d-%dmode'%(modes[0],modes[-1]),
                           label_list='', output=outputroot)

psa.plot_2d_power_spectrum(pow_select, None, k_p_edges, k_v_edges,
                           filename='comp_cros_ps_%d-%dmode'%(modes[0],modes[-1]),
                           label_list='', output=outputroot)

psa.plot_2d_power_spectrum(transfer_select, None, k_p_edges, k_v_edges,
                           filename='transferf_%d-%dmode'%(modes[0],modes[-1]),
                           label_list='', output=outputroot, logscale=False,
                           cmax=1., cmin=0.3)

var = std_select**2
var[var == 0] = np.inf
var = var**(-1)
print power_spectrum_sim.info['axes']
print sim_select.shape
sim_list = []
for arr in sim_select:
    sim = al.make_vect(arr, axis_names=power_spectrum_sim.info['axes'])
    sim.info = power_spectrum_sim.info
    sim_list.append(sim)
si_1d_mean, si_1d_std, si_k_centre = pss.convert_2dps_to_1dps_each(
        sim_list, sim_kmn_list, var, None)
si_positive, si_negative = psa.seperate_positive_and_negative_power(
        si_1d_mean, si_1d_std, si_k_centre)
filename = '%scros_si_%d-%dmode_1dpow_positive.txt'%('', modes[0], modes[-1])
np.savetxt('%s/%s'%(outputroot, filename), si_positive.T, fmt='%15.12e')
filename = '%scros_si_%d-%dmode_1dpow_negative.txt'%('', modes[0], modes[-1])
np.savetxt('%s/%s'%(outputroot, filename), si_negative.T, fmt='%15.12e')

ps_1d_mean, ps_1d_std, k_centre = pss.convert_2dps_to_1dps_each(
        power_spectrum, k_mode_number, var, None)

positive, negative = psa.seperate_positive_and_negative_power(ps_1d_mean,
                                                              ps_1d_std,
                                                              k_centre)
filename = '%scros_ps_%d-%dmode_1dpow.txt'%('comp_varweight_', modes[0], modes[-1])
np.savetxt('%s/%s'%(outputroot, filename), np.concatenate((k_centre[None,:],\
        ps_1d_mean[None,:], si_1d_std[None,:]),axis=0).T, fmt='%15.12e')
filename = '%scros_ps_%d-%dmode_1dpow_positive.txt'%('comp_varweight_', modes[0], modes[-1])
np.savetxt('%s/%s'%(outputroot, filename), positive.T, fmt='%15.12e')
filename = '%scros_ps_%d-%dmode_1dpow_negative.txt'%('comp_varweight_', modes[0], modes[-1])
np.savetxt('%s/%s'%(outputroot, filename), negative.T, fmt='%15.12e')

