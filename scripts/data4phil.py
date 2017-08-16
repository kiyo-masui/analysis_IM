import numpy as np
import h5py
from core import algebra as al

out = h5py.File('forPhil.hd5', 'a')

names = ['301', '302', '303', '304', '305']

for ii in names:
    inp = h5py.File('full%s_cros_ps_10mode/'%ii + 'ps_result.hd5', 'r')
    ps_2d = al.make_vect(al.load_h5(inp, 'cros_ps_10mode_2dpow'))
    al.save_h5(out, 'cros_ps_2d_%s'%ii, ps_2d)
    weight = al.make_vect(np.loadtxt('full%s_cros_ps_10mode/'%ii + 'weights'))
    weight.info = ps_2d.info
    al.save_h5(out, 'cros_std_2d_%s'%ii, weight)

k_edges = al.make_vect(np.load('full301_cros_ps_10mode/' + 'k_vals.npy'))
k_edges.info = ps_2d.info
al.save_h5(out, 'k_edges', k_edges)
