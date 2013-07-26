import numpy as np
import core.algebra as al
import random


def shot_noise_effect(delta, n_bar, bins):
	r_collect = {}
	for n in n_bar:
		print '--- n_bar = ',n,'---'
		proba = n*delta + n
		gal = np.random.poisson(proba)
		delta_sim = al.make_vect(gal/gal.mean() - 1)
		delta_sim.info = delta.info
		dss = MapBuilder(delta_sim)
		p_collect = power_specs(dss,bins)
		r_collect[n] = p_collect['r']
	return r_collect
