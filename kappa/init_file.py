from kappa.estimator import *
from kappa.prepare_map import *
import matplotlib.pyplot as plt


# Directory for 2dF overdensity map & selection function
path_2df = '/mnt/raid-project/gmrt/goblot/2df/'

# These maps were made based on the 2dF catalogs which can be found in :
path_cat = '/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_full_selection_10000mock/'
# (you have to work on prawn to access these original catalogs)

# Parameters of the reconstruction
params = {
    'sigma_recon': 1.25,
    'sigma_smooth': 8.,

    'clean': True,
    'smoothing': True,

    'interm_step': {'map_log': True, 'quad_est': False, 'kernel': False},
    'save': ['kappa', 'weight', 'clean', 'mask'],

    'fast': False,
    'auto_rescale': True,
    'verbose': True
}

# Example of what the info should be (this is 2df info)
#info = {
#    'axes': ('z', 'x', 'y'),
#    'z_centre': 494.7297656238655,
#    'x_centre': 0.0,
#    'y_centre': 0.0,
#    'z_delta': 2.6650922751987474,
#    'x_delta': 0.7967856053857366,
#    'y_delta': 0.9240468503362674,
#    'type': 'vect'
#    }


input_file = path_2df + 'delta_ycli_xyz.npy'

output_dir = path_2df
prefix = 'example_'
out_type = 'numpy'
#out_type = 'binary'


if __name__ == '__main__':
    print "Loading and preprocessing maps\n"
    delta_2df = al.load(input_file)
    delta_2df = preprocess(delta_2df[3:103])
    # we use delta_2df[3:103] because we don't want to use the full box (too noisy for high z)

    sel_func = al.load(path_2df+'sel_func_ycli_xyz.npy')
    mask = sel_func > 0.3*sel_func.max()
    mask = mask[3:103] # This is to have the same box as delta_2df

    print "Running estimator :"
    kappa_est = KappaEstimator(delta_2df, params, mask=mask)
    kappa_est.run()

    print ""
    kappa_est.save(output_dir, prefix, out_type)
