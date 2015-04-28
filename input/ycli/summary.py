#! /usr/bin/env python

import os
import string
import sys

workroot = os.getenv('YICHAO_WORK_PS')

if len(sys.argv)==2:
	if sys.argv[1] == 'auto':
		map_sim = os.getenv('MAP_SIM')
		map_ssm = os.getenv('MAP_SSM')
		map_cln = os.getenv('MAP_CLN')
		mode = int(os.getenv('MODE'))
		hour = int(os.getenv('HOUR'))
		output_root = workroot + 'power_auto_%s_%d/'%(map_cln, mode)
		output_file = 'auto_%s_%d.summary'%(map_cln, mode)

		f = open( output_root + output_file, 'w')
		f.write('summary for auto power spectrum with %dmodes subtracted\n'%mode)
		f.write('for reference, using:\n')
		f.write('%s_temperature x %s_temperature\n'%(map_sim, map_sim))
		f.write('for transfre function, using:\n')
		f.write('%s_combined x %s_combined\n'%(map_ssm, map_ssm))
		f.write('for power spctrum, using:\n')
		f.write('%s x %s'%(map_cln, map_cln))

		f.close()
	elif sys.argv[1] == 'cros':
		map_sim = os.getenv('MAP_SIM')
		map_ssm = os.getenv('MAP_SSM')
		map_cln = os.getenv('MAP_CLN')
		map_wgz = os.getenv('MAP_WGZ')
		mode = int(os.getenv('MODE'))
		hour = int(os.getenv('HOUR'))
		output_root = workroot + 'power_cros_%s_%d/'%(map_cln, mode)
		output_file = 'cros_%s_%d.summary'%(map_cln, mode)

		f = open( output_root + output_file, 'w')
		f.write('summary for cross power spectrum with %dmodes subtracted\n'%mode)
		f.write('for reference, using:\n')
		f.write('%s_temperature x %s_delta\n'%(map_sim, map_sim))
		f.write('for transfre function, using:\n')
		f.write('%s_combined x %s_delta\n'%(map_ssm, map_sim))
		f.write('for power spctrum, using:\n')
		f.write('%s_combined x %s_delta_binned_data'%(map_cln, map_wgz))

		f.close()
	elif sys.argv[1] == 'wigl':
		map_wgz = os.getenv('MAP_WGZ')
		hour = int(os.getenv('HOUR'))
		output_root = workroot + 'power_auto_wiggle_z_%dhr/'%hour
		output_file = 'auto_wiggle_z_%dhr.summary'%hour

		f = open( output_root + output_file, 'w')
		f.write('summary for wigglez auto power spectrum in %dhr\n'%hour)
		f.write('for power spctrum, using:\n')
		f.write('%s x %s'%(map_wgz, map_wgz))
		f.close()
