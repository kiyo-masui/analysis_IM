from kappa.estimator import *
from kappa.prepare_map import *


path_2df = '/mnt/raid-project/gmrt/goblot/2df/'
path_cat = '/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_full_selection_10000mock/' #you have to work on prawn to access these original catalogs

params = {
	'sigma_recon': 1.25, 
	'sigma_smooth': 8., 

	'clean': True,
	'smoothing': True,
	
	'save_map_log': True,
	'save_quad_est': True,
	'save_kernel': True,
	
	'fast': False, 
	'auto_rescale': True,
	'verbose': True
	}
	
info = {
	'axes': ('z', 'x', 'y'),
	'z_centre': 494.7297656238655,
	'x_centre': 0.0, 	 
	'y_centre': 0.0, 
	'z_delta': 2.6650922751987474, 
	'x_delta': 0.7967856053857366, 
	'y_delta': 0.9240468503362674, 
	'type': 'vect'
	}

input_file = path_2df + 'delta_ycli_xyz.npy'
output_file = path_2df + 'kappa_example'


if __name__=='__main__':
	print "Loading and preprocessing maps\n"
	delta_2df = al.load(input_file)
	delta_2df = preprocess(delta_2df[3:103])
	
	sel_func = al.load(path_2df+'sel_func_ycli_xyz.npy')
	mask = sel_func > 0.3*sel_func.max()
	mask = mask[3:103]
	
	print "Running estimator :"
	kappa_est = KappaEstimator(delta_2df,params, mask=mask)
	kappa_est.run()
	
#	kappa_est.save(output_file)
#	al.save(output_file+'.npy', kappa_est.kappa)
#	al.save(output_file+'.npy', kappa_est.clean)
	