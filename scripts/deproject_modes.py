#Script for deprojecting modes from maps or vectors.

import numpy as np
import core.algebra as al
import h5py

def left_deproject(mat, vec):
    #First dimension of mat must equal dim of vec.
    #Normalize vec.
    vec /= np.sqrt(np.dot(vec, vec))
    proj = np.transpose(np.transpose(np.ones(mat.shape)*np.tensordot(mat,vec, axes=[0,0]))*vec)
    ans = mat - proj
    return ans

if __name__ == '__main__':
    modes_dir = '/scratch2/p/pen/andersoc/GBT_cleaning_results/1hr/1hr_80-68_azrm_short_svd_flagged5/chris_cuts_131chan/I/'
    orig_dir = '/scratch2/p/pen/andersoc/GBT_maps/1hr/map_differences_1hr_80-68_azrm_short_svd_flagged5/' 
    new_dir = orig_dir + 'chris_cuts_deproject2modes/'
    weight_dir = '/scratch2/p/pen/andersoc/GBT_maps/1hr/'
    freqs = np.arange(256)
    freqs = np.delete(freqs, range(0,63) + [103,104,105,106,107,108,130,131,132,133,134] + range(179,218) + range(246,256))
    noise_list = ['A_minus_B', 'A_minus_C',  'A_minus_D', 'B_minus_C', 'B_minus_D', 'C_minus_D']
    key_list = ['beam1_with_beam2','beam1_with_beam3','beam1_with_beam4','beam2_with_beam3','beam2_with_beam4','beam3_with_beam4']
    
    noise_mode_dir = '/scratch2/p/pen/andersoc/GBT_cleaning_results/1hr/map_diff_autos_1hr_80-68_azrm_short_svd_flagged5/chris_cuts_deproj_mapmodes/I/'
    orig_map_dir = '/scratch2/p/pen/andersoc/GBT_maps/1hr/'

    weights = ['A', 'B', 'C', 'D']
    '''    for let in weights:
        weight = al.load(weight_dir + 'sec' + let + '_1hr_80-68_azrm_short_svd_flagged5_noise_inv_diag_I_800.npy')
        #weight = weight[freqs,:,:]
        al.save(new_dir + 'sec' + let + '_1hr_80-68_azrm_short_svd_flagged5_noise_inv_diag_I_800.npy', weight)
    '''    

    for pair in noise_list[1:]:
        #map1 = al.load(orig_map_dir + 'sec' + pair[0] + '_1hr_80-68_azrm_short_svd_flagged5_clean_map_I_800.npy')
        #map2 = al.load(orig_map_dir + 'sec' + pair[-1] + '_1hr_80-68_azrm_short_svd_flagged5_clean_map_I_800.npy')
        map1 = al.load(orig_map_dir + 'secA' + '_1hr_80-68_azrm_short_svd_flagged5_clean_map_I_800.npy')
        map2 = al.load(orig_map_dir + 'secB' + '_1hr_80-68_azrm_short_svd_flagged5_clean_map_I_800.npy')
        map3 = al.load(orig_map_dir + 'secC' + '_1hr_80-68_azrm_short_svd_flagged5_clean_map_I_800.npy')
        map4 = al.load(orig_map_dir + 'secD' + '_1hr_80-68_azrm_short_svd_flagged5_clean_map_I_800.npy')
        svd = h5py.File(noise_mode_dir + 'SVD.hd5')
        for ind in range(30):
            mode_short = svd['svd_modes1'][pair][ind]
            mode = np.zeros(256)
            mode[63:103] = mode_short[0:40]
            mode[109:130] = mode_short[40:61]
            mode[135:179] = mode_short[61:105]
            mode[218:246] = mode_short[105:133]
            map1 = left_deproject(map1, mode)
            map2 = left_deproject(map2, mode)
            map3 = left_deproject(map3, mode)
            map4 = left_deproject(map4, mode)
        #al.save(orig_map_dir + 'deproj_30noisemodes/' + pair[0] + pair[-1] + '/' + 'sec' + pair[0] + '_1hr_80-68_azrm_short_svd_flagged5_clean_map_I_800.npy', map1)
        #al.save(orig_map_dir + 'deproj_30noisemodes/' + pair[0] + pair[-1] + '/' + 'sec' + pair[-1] + '_1hr_80-68_azrm_short_svd_flagged5_clean_map_I_800.npy', map2)
        al.save(orig_map_dir + 'deproj_30noisemodes/' + pair[0] + pair[-1] + '/' + 'secA' + '_1hr_80-68_azrm_short_svd_flagged5_clean_map_I_800.npy', map1)
        al.save(orig_map_dir + 'deproj_30noisemodes/' + pair[0] + pair[-1] + '/' + 'secB' + '_1hr_80-68_azrm_short_svd_flagged5_clean_map_I_800.npy', map2)
        al.save(orig_map_dir + 'deproj_30noisemodes/' + pair[0] + pair[-1] + '/' + 'secC' + '_1hr_80-68_azrm_short_svd_flagged5_clean_map_I_800.npy', map3)
        al.save(orig_map_dir + 'deproj_30noisemodes/' + pair[0] + pair[-1] + '/' + 'secD' + '_1hr_80-68_azrm_short_svd_flagged5_clean_map_I_800.npy', map4)

    '''    for pair in zip(noise_list, key_list):
        name = pair[0]
        key = pair[1]
        noise = al.load(orig_dir + name)
        #noise = noise[freqs,:,:]
        svd_modes = h5py.File(modes_dir + 'SVD.hd5')
        mode1_short = svd_modes['svd_modes1'][key][0]
        mode1 = np.zeros(256)
        mode1[63:103] = mode1_short[0:40]
        mode1[109:130] = mode1_short[40:61]
        mode1[135:179] = mode1_short[61:105]
        mode1[218:246] = mode1_short[105:133]
        mode2_short = svd_modes['svd_modes1'][key][1]
        mode2 = np.zeros(256)
        mode2[63:103] = mode2_short[0:40]
        mode2[109:130] = mode2_short[40:61]
        mode2[135:179] = mode2_short[61:105]
        mode2[218:246] = mode2_short[105:133]
        noise = left_deproject(noise, mode1)
        noise = left_deproject(noise, mode2)
        mode1_short = svd_modes['svd_modes2'][key][0]
        mode1 = np.zeros(256)
        mode1[63:103] = mode1_short[0:40]
        mode1[109:130] = mode1_short[40:61]
        mode1[135:179] = mode1_short[61:105]
        mode1[218:246] = mode1_short[105:133]
        mode2_short = svd_modes['svd_modes2'][key][1]
        mode2 = np.zeros(256)
        mode2[63:103] = mode2_short[0:40]
        mode2[109:130] = mode2_short[40:61]
        mode2[135:179] = mode2_short[61:105]
        mode2[218:246] = mode2_short[105:133]
        noise = left_deproject(noise, mode1)
        noise = left_deproject(noise, mode2)
        al.save(new_dir + name, noise)
    ''' 
