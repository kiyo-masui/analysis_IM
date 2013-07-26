import numpy as np
import core.algebra as al
import random
import copy as cp
import scipy.interpolate as si
from utils import units
from utils import cosmology as cosmo


def radius_binning(map,bins,mask,info):
	l,m,n = map.shape
	z = (np.arange(l) - l/2)*info['z_delta'] + info['z_centre']
	x = (np.arange(m) - m/2)*info['x_delta'] + info['x_centre']
	y = (np.arange(n) - n/2)*info['y_delta'] + info['y_centre']
	rad = (z*z)[:,None,None] * np.ones((1,m,n))
	rad += (x*x)[None,:,None] * np.ones((l,1,n))
	rad += (y*y)[None,None,:] * np.ones((l,m,1))
	rad = np.sqrt(rad)
	output = np.histogram(rad,bins=bins,weights=map)[0]
	count = np.histogram(rad,bins=bins,weights=mask)[0]
	return output,count


def radius_weighting(weight_1d,bins,map,info):
	l,m,n = map.shape
	z = (np.arange(l) - l/2)*info['z_delta'] + info['z_centre']
	x = (np.arange(m) - m/2)*info['x_delta'] + info['x_centre']
	y = (np.arange(n) - n/2)*info['y_delta'] + info['y_centre']
	rad = (z*z)[:,None,None] * np.ones((1,m,n))
	rad += (x*x)[None,:,None] * np.ones((l,1,n))
	rad += (y*y)[None,None,:] * np.ones((l,m,1))
	rad = np.sqrt(rad)
	w = si.interp1d(bins,weight_1d)
	valid_interp = (rad>=bins[0])&(rad<=bins[-1])
	weight_3d = np.zeros_like(map)
	weight_3d[valid_interp] = w(rad[valid_interp])
	return map*weight_3d



if __name__ == 'main':
	path_cat = '/mnt/raid-project/gmrt/goblot/JD_catalogs/'
	path_temp = '/mnt/raid-project/gmrt/goblot/'
	path = '/mnt/scratch-gl/ycli/2df_catalog/map/' # works only on prawn! if on other machine, switch with comments for eq_coord_info, f and ratio_xyz
	
	cat_ycli_xyz = np.load(path_temp+'2df/cat_ycli_xyz.npy')
	sel_func_ycli_xyz = np.load(path_temp+'2df/sel_func_ycli_xyz.npy')
	map_ycli_xyz = al.load(path_temp+'2df/delta_ycli_xyz.npy')
	map_ycli_xyz = al.make_vect(map_ycli_xyz)
	info = cp.deepcopy(map_ycli_xyz.info)
	
	cat_ycli = al.load(path+'map_2929.5_full_selection_10000mock/real_map_2df.npy')
	eq_coord_info = cp.deepcopy(cat_ycli.info)
		
#	eq_coord_info = {
#		'axes': ('freq', 'ra', 'dec'),
# 		'dec_centre': -29.5,
# 		'dec_delta': 0.35,
# 		'freq_centre': 1221900000.0,
# 		'freq_delta': -1000000.0,
# 		'ra_centre': 7.0,
# 		'ra_delta': -0.35,
# 		'type': 'vect'
# 	}
 	
	f = cat_ycli.get_axis('freq')
#	l = cat_ycli.shape[0]
#	f = (np.arange(l)-l/2)*eq_coord_info['freq_delta'] + eq_coord_info['freq_centre']

	ratio_xyz = cat_ycli_xyz.sum()/cat_ycli.sum()
#	ratio_xyz = cat_ycli_xyz.sum()/58698.0
	
	cosmology = cosmo.Cosmology()
	z_21 = units.nu21*1e6/f - 1
	d_21 = cosmology.comoving_distance(z_21)
	d_21_shift = np.zeros_like(d_21)
	d_21_shift[:-1] = d_21[1:]
	d_21_shift[-1] = d_21[-1] + (d_21[-1]-d_21[-2])
	
	vol_pix = info['z_delta']*info['x_delta']*info['y_delta']
	vol_ang = -eq_coord_info['dec_delta']*np.pi/180.*eq_coord_info['ra_delta']*np.pi/180*(d_21**2)*(d_21_shift - d_21)*np.cos(eq_coord_info['dec_centre']*np.pi/180)

	cat_ycli_xyz_w = radius_weighting(vol_pix/vol_ang,d_21,cat_ycli_xyz,info)
	sel_func_ycli_xyz_w = radius_weighting(vol_pix/vol_ang,d_21,sel_func_ycli_xyz,info)

	sim_2df_like = al.load(path_cat+'gal_10_in_2df_box.npy')
	nul = sim_2df_like==0
	n_bar = sim_2df_like.sum()/(vol_pix*1384*261*177)

	proba_r = (sel_func_ycli_xyz/ratio_xyz)/(n_bar*vol_pix)
	mock_r = np.random.binomial(np.int_(sim_2df_like+nul),np.minimum(proba,1))*(1-nul)
	
	proba_w = sel_func_ycli_xyz_w/(n_bar*vol_pix)
	mock_w = np.random.binomial(np.int_(sim_2df_like+nul),np.minimum(proba_w,1))*(1-nul)
	double_select = np.random.binomial(np.int_(sim_2df_like+nul),np.maximum(proba_w-1,0))*(1-nul)
	mock_w += double_select
	
	delta_r = ratio_xyz*mock_r/sel_func_ycli_xyz - 1
	delta_r[np.isnan(delta_r)]=0
	
	delta_w = mock_w/sel_func_ycli_xyz_w - 1
	delta_w[np.isnan(delta_w)]=0
	
	mask_2df = sel_func_ycli_xyz > 0
	mask = sel_func_ycli_xyz > 0.3*sel_func_ycli_xyz.max()
	
	z_phys = map_ycli_xyz.get_axis('z')
	bins_nz = (np.arange(100))*2*info['z_delta'] + z_phys[0]
	plot_bins_nz = (bins_nz[1:]+bins_nz[:-1])/2

	n_z_sel_func_r = radius_binning(sel_func_ycli_xyz/ratio_xyz*mask_2df,bins_nz,(mask_2df+0),info)
	n_z_2df_r = radius_binning(cat_ycli_xyz/ratio_xyz*mask_2df,bins_nz,(mask_2df+0),info)
	n_z_mock_r = radius_binning(mock_r*mask_2df,bins_nz,(mask_2df+0),info)
	
	n_z_sel_func_w = radius_binning(sel_func_ycli_xyz_w*mask_2df,bins_nz,(mask_2df+0),info)
	n_z_2df_w = radius_binning(cat_ycli_xyz_w*mask_2df,bins_nz,(mask_2df+0),info)
	n_z_mock_w = radius_binning(mock_w*mask_2df,bins_nz,(mask_2df+0),info)


#	al.save(path_temp+'2df_mock/mock_2df_ratio.npy',mock_r)
#	al.save(path_temp+'2df_mock/mock_2df_delta_ratio.npy',delta_r)
#	al.save(path_temp+'2df_mock/mock_2df_weighted.npy',mock_w)
#	al.save(path_temp+'2df_mock/mock_2df_delta_weighted.npy',delta_w)
	