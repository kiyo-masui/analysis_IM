#! /usr/bin/env python 

from parkes import plotmap

plot_opt = False
plot_raw = True
plot_gbt = False
plot_cln = False

if plot_gbt:

    prefix = 'secA'
    suffix = '800'
    field = '15hr_41-80_avg_fdgp_new'

    degrade_factor = 1.4
    degrade_map = True
    map_root = '/Users/ycli/DATA/gbt_raw_map/15hr/'
    mapdict = {}
    mapdict['imap'] = map_root + "%s_%s_clean_map_I_%s.npy"%(prefix, field, suffix)
    mapdict['nmap'] = map_root + "%s_%s_noise_inv_diag_I_%s.npy"%(prefix, field, suffix)
    mapdict['name'] = "%s_%s_clean_map"%(prefix, field)
    mapdict['noise_name'] = "%s_%s_noise_map"%(prefix, field)
    a = plotmap.PlotMap(mapdict, 
            freq_cut = [0,15,16,17,18,19,20,21,22,23,103,104,105,106,107,
                        108,130,131,132,133,134],
            degrade_factor = degrade_factor,
            degrade_map = degrade_map,
            plot_size=(8,6),
            nvss_catalog_path = "/Users/ycli/DATA/nvss/nvss_catlog_15hr.dat")
    a.plot_map(with_nvss=True)
    a.check_nvss(flux_limit=[0.3, 100], rescale=False)

if plot_opt:

    prefix = 'pks'
    field = 'p3500n3000_parkes_2010_10_24-28'
    degrade_factor = 1.1
    degrade_map = True
    map_root = '/home/ycli/data/2df/map/%s_%s/'%(prefix, field)
    mapdict = {}
    mapdict['imap'] = map_root + 'real_map_2df.npy'
    #mapdict['nmap'] = map_root + 'sele_map_2df.npy'
    mapdict['nmap'] = map_root + 'sele_map_2df_separable.npy'
    mapdict['name'] = '%s_%s_opt_real_map'%(prefix, field)
    mapdict['noise_name'] = '%s_%s_opt_sele_map'%(prefix, field)
    a = plotmap.PlotMap(mapdict, 
                        freq_cut = [0,1,2,3,4,5,59,60,61,62,63], 
                        degrade_factor=degrade_factor, 
                        degrade_map=degrade_map,
                        plot_size=(15, 6))
    a.plot_map(with_nvss=False)
    a.check_nvss(flux_limit=[0.5, 100], rescale=False)

if plot_raw:
    #prefix = 'test_allbeams'
    #suffix = '1315'
    #field = '27n30_10by7'
    #degrade_factor = 1.1
    #degrade_map = False
    #map_root = '/home/ycli/data/parkes/maps_by_chris/pks_%s/'%field
    
    #prefix = 'fir'
    #suffix = '1316'
    #field = 'p3500n3000_parkes_2010_10_24-28'
    #degrade_factor = 1.1
    #degrade_map = False
    #map_root = '/home/ycli/data/parkes/maps/pks_%s/'%field
    
    #prefix = 'fir'
    #suffix = '1316'
    ##field = 'testflagging_p3500n3000_parkes_2010_10_24'
    ##field = 'testflagging_old_p3500n3000_parkes_2010_10_24'
    #field = 'testflagging_p3500n3000_parkes_2010_10_24'
    #degrade_factor = 1.1
    #degrade_map = False
    ##map_root = '/home/ycli/data/parkes/maps/parkes_testflagging/'
    #map_root = '/home/ycli/data/parkes/maps/parkes_testflagging/threshold_3pt5/'
    ##map_root = '/home/ycli/data/parkes/maps/parkes_testpointing/threshold_3pt5/'

    #suffix = '800'
    ##field = '15hr_41-80_pointcorr'
    #field = '15hr_41-80_avg_fdgp_new'
    #degrade_factor = 1.4
    #degrade_map = True
    #map_root = '/home/ycli/data/gbt/gbt_%s/'%field

    prefix = 'fir'
    suffix = '1316'

    field = 'p3500n3000_parkes_2010_10_24-28'
    #field = 'p3500n3000_parkes_2010_10_24-28_beam_0'
    #field = 'p3500n3000_parkes_2010_10_24-28_beam_1'
    #field = 'p3500n3000_parkes_2010_10_24-28_beam_2'
    #field = 'p3500n3000_parkes_2010_10_24-28_beam_5'
    #field = 'p3500n3000_parkes_2010_10_24-28_beam_6'
    #field = 'p3500n3000_parkes_2010_10_24-28_beam_7'
    #field = 'p3500n3000_parkes_2010_10_24-28_beam_8'
    #field = 'p3500n3000_parkes_2010_10_24-28_beam_9'

    degrade_factor = 1.1
    degrade_map = False
    #map_root = '/scratch/p/pen/ycli/map_result/maps/parkes_parallel/'
    #map_root = '/Users/ycli/DATA/parkes/maps/pks_p3500n3000_parkes_2010_10_24-28/'
    #map_root = '/Users/ycli/DATA/parkes/maps/small_pks_p3500n3000_parkes_2010_10_24-28/'
    #map_root = '/Users/ycli/DATA/parkes/maps/pks_p3500n3000_parkes_2010_10_24-28/'
    map_root = '/Users/ycli/DATA/parkes/maps/pks_p3500n3000_parkes_2010_10_24-28_in_J/'
    
    mapdict = {}
    mapdict['imap'] = map_root + '%s_%s_clean_map_I_%s.npy'%(prefix, field, suffix)
    mapdict['nmap'] = map_root + '%s_%s_noise_inv_diag_I_%s.npy'%(prefix, field, suffix)
    mapdict['name'] = '%s_clean_map_I_in_J'%(field)
    mapdict['noise_name'] = '%s_noise_map_I_in_J'%(field)

    #mapdict['imap'] = map_root + 'secA_%s_clean_map_I_%s.npy'%(field, suffix)
    ##mapdict['nmap'] = map_root + 'secA_%s_noise_weight_I_%s.npy'%(field, suffix)
    #mapdict['nmap'] = map_root + 'secA_%s_noise_inv_diag_I_%s.npy'%(field, suffix)
    #mapdict['name'] = '%s_clean_map_I'%(field)
    #mapdict['noise_name'] = '%s_noise_map_I'%(field)

    #mapdict['imap_sec'] = [
    #        map_root + '%s_%s_A_clean_map_I_%s.npy'%(prefix, field, suffix),
    #        map_root + '%s_%s_B_clean_map_I_%s.npy'%(prefix, field, suffix), 
    #        map_root + '%s_%s_C_clean_map_I_%s.npy'%(prefix, field, suffix), 
    #        map_root + '%s_%s_D_clean_map_I_%s.npy'%(prefix, field, suffix), 
    #        map_root + '%s_%s_E_clean_map_I_%s.npy'%(prefix, field, suffix), 
    #                      ]
    #mapdict['nmap_sec'] = [
    #        map_root + '%s_%s_A_noise_weight_I_%s.npy'%(prefix, field, suffix),
    #        map_root + '%s_%s_B_noise_weight_I_%s.npy'%(prefix, field, suffix), 
    #        map_root + '%s_%s_C_noise_weight_I_%s.npy'%(prefix, field, suffix), 
    #        map_root + '%s_%s_D_noise_weight_I_%s.npy'%(prefix, field, suffix), 
    #        map_root + '%s_%s_E_noise_weight_I_%s.npy'%(prefix, field, suffix), 
    #                      ]
    #mapdict['name'] += '_1pt1_cov'
    #mapdict['noise_name'] += '_1pt1_cov'
    
    a = plotmap.PlotMap(mapdict, 
                        #freq_cut = [0,1,2,3,4,5,59,60,61,62,63], 
                        degrade_factor=degrade_factor, 
                        degrade_map=degrade_map,
                        #plot_size=(10, 8))
                        #plot_size=(12, 6))
                        plot_size=(12, 4))
    #a.mapping_coord(plot=True)
    a.plot_map(with_nvss=True)
    #a.plot_map(with_nvss=True, diff=True)
    #a.check_nvss(flux_limit=[0.5, 100], rescale=False)
    #a.check_nvss(flux_limit=[0.5, 100])
    #a.plot_map(with_nvss=True)
    #a.plot_map(with_nvss=True, diff=True)

if plot_cln:
    #mode = 2 
    #prefix = 'test_allbeams'
    #field = '27n30_10by7'
    #clean = 'ABCD_freq_cut_1pt1_cov_bgcal'
    #degrade_factor = 1.1
    #degrade_map = True
    #map_root = '/home/ycli/data/cln_result/PKS_%s_%s/'%(field, clean)
    #map_root = map_root + 'Emap_clean_themselves_noise_prior2En1/'

    mode = 5 
    prefix = 'fir'
    field = 'p3500n3000_parkes_2010_10_24-28'
    #clean = 'ABCD_freq_cut_1pt1_cov_cal'
    clean = 'ABCDE_freq_cut_1pt1_cov_bgcal_new'
    degrade_factor = 1.1
    degrade_map = True
    map_root = '/home/ycli/data/cln_result/PKS_%s_%s/'%(field, clean)
    map_root = map_root + 'Emap_clean_themselves/'

    #mode = 25
    #prefix = 'gbt'
    ##field = '15hr_ptcorr_finefreq'
    #field = '15hr_41-80_avg_fdgp_new'
    ###clean = 'ABCD_freq_cut_1pt1_cov_cal'
    ##clean = 'ABCDE_freq_cut_1pt1_cov_bgcal'
    #degrade_factor = 0
    #degrade_map = True
    ##map_root = '/home/ycli/data/cln_result/15hr_AA_fine_freq_11conv/'
    #map_root = '/home/ycli/data/cln_result/GBT_15hr_41-80_avg_fdgp_new_ABCD_1pt4_cov/'
    #map_root = map_root + 'Emap_clean_themselves/'

    mapdict = {}
    mapdict['imap'] = map_root + 'combined_clean_map_%dmodes.npy'%mode
    mapdict['nmap'] = map_root + 'combined_clean_weight_%dmodes.npy'%mode
    mapdict['name'] = '%s_%s_clean_map_I_%dmode_1pt1_cov_cal'%(prefix, field, mode)
    mapdict['noise_name'] = '%s_%s_noise_map_I_%dmode_1pt1_cov'%(prefix, field, mode)
    #mapdict['name'] = '%s_%s_clean_map_I_%dmode_1pt1_cov_cal_noise_prior2En1'%(prefix, field, mode)
    #mapdict['noise_name'] = '%s_%s_noise_map_I_%dmode_1pt1_cov_cal_noise_prior2En1'%(prefix, field, mode)
    mapdict['imap_sec'] = [
            #map_root + 'sec_A_cleaned_clean_map_I_with_A_%dmodes.npy'%mode,

            map_root + 'sec_A_cleaned_clean_map_I_with_B_%dmodes.npy'%mode,
            map_root + 'sec_A_cleaned_clean_map_I_with_C_%dmodes.npy'%mode, 
            map_root + 'sec_B_cleaned_clean_map_I_with_A_%dmodes.npy'%mode, 
            map_root + 'sec_B_cleaned_clean_map_I_with_C_%dmodes.npy'%mode, 
            map_root + 'sec_C_cleaned_clean_map_I_with_A_%dmodes.npy'%mode, 
            map_root + 'sec_C_cleaned_clean_map_I_with_B_%dmodes.npy'%mode, 

            map_root + 'sec_A_cleaned_clean_map_I_with_D_%dmodes.npy'%mode, 
            map_root + 'sec_B_cleaned_clean_map_I_with_D_%dmodes.npy'%mode, 
            map_root + 'sec_C_cleaned_clean_map_I_with_D_%dmodes.npy'%mode, 
            map_root + 'sec_D_cleaned_clean_map_I_with_A_%dmodes.npy'%mode, 
            map_root + 'sec_D_cleaned_clean_map_I_with_B_%dmodes.npy'%mode, 
            map_root + 'sec_D_cleaned_clean_map_I_with_C_%dmodes.npy'%mode, 

            map_root + 'sec_A_cleaned_clean_map_I_with_E_%dmodes.npy'%mode, 
            map_root + 'sec_B_cleaned_clean_map_I_with_E_%dmodes.npy'%mode, 
            map_root + 'sec_C_cleaned_clean_map_I_with_E_%dmodes.npy'%mode, 
            map_root + 'sec_D_cleaned_clean_map_I_with_E_%dmodes.npy'%mode, 
            map_root + 'sec_E_cleaned_clean_map_I_with_A_%dmodes.npy'%mode, 
            map_root + 'sec_E_cleaned_clean_map_I_with_B_%dmodes.npy'%mode, 
            map_root + 'sec_E_cleaned_clean_map_I_with_C_%dmodes.npy'%mode, 
            map_root + 'sec_E_cleaned_clean_map_I_with_D_%dmodes.npy'%mode, 
            ]
    mapdict['nmap_sec'] = [
            #map_root + 'sec_A_cleaned_noise_inv_I_with_A_%dmodes.npy'%mode,

            map_root + 'sec_A_cleaned_noise_inv_I_with_B_%dmodes.npy'%mode,
            map_root + 'sec_A_cleaned_noise_inv_I_with_C_%dmodes.npy'%mode, 
            map_root + 'sec_B_cleaned_noise_inv_I_with_A_%dmodes.npy'%mode, 
            map_root + 'sec_B_cleaned_noise_inv_I_with_C_%dmodes.npy'%mode, 
            map_root + 'sec_C_cleaned_noise_inv_I_with_A_%dmodes.npy'%mode, 
            map_root + 'sec_C_cleaned_noise_inv_I_with_B_%dmodes.npy'%mode, 

            map_root + 'sec_A_cleaned_noise_inv_I_with_D_%dmodes.npy'%mode, 
            map_root + 'sec_B_cleaned_noise_inv_I_with_D_%dmodes.npy'%mode, 
            map_root + 'sec_C_cleaned_noise_inv_I_with_D_%dmodes.npy'%mode, 
            map_root + 'sec_D_cleaned_noise_inv_I_with_A_%dmodes.npy'%mode, 
            map_root + 'sec_D_cleaned_noise_inv_I_with_B_%dmodes.npy'%mode, 
            map_root + 'sec_D_cleaned_noise_inv_I_with_C_%dmodes.npy'%mode, 

            map_root + 'sec_A_cleaned_noise_inv_I_with_E_%dmodes.npy'%mode, 
            map_root + 'sec_B_cleaned_noise_inv_I_with_E_%dmodes.npy'%mode, 
            map_root + 'sec_C_cleaned_noise_inv_I_with_E_%dmodes.npy'%mode, 
            map_root + 'sec_D_cleaned_noise_inv_I_with_E_%dmodes.npy'%mode, 
            map_root + 'sec_E_cleaned_noise_inv_I_with_A_%dmodes.npy'%mode, 
            map_root + 'sec_E_cleaned_noise_inv_I_with_B_%dmodes.npy'%mode, 
            map_root + 'sec_E_cleaned_noise_inv_I_with_C_%dmodes.npy'%mode, 
            map_root + 'sec_E_cleaned_noise_inv_I_with_D_%dmodes.npy'%mode, 
            ]
    a = plotmap.PlotMap(mapdict, 
                        freq_cut = [0,1,2,3,4,5,59,60,61,62,63], 
                        degrade_factor=degrade_factor, 
                        degrade_map=degrade_map,
                        plot_size=(12, 6))
                        #plot_size=(15, 6))
    #a.mapping_coord(plot=True)
    a.plot_map(with_nvss=False)
    #a.plot_map(with_nvss=True, diff=True)
    #a.check_nvss(flux_limit=[0.5, 100], rescale=False)


