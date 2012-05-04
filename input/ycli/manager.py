#! /usr/bin/env python
"""This module contains the pipeline manager.

It is simply a driver for scripting together many pipeline modules.
"""

import os
import sys

from kiyopy import parse_ini

params_init = {
		# A list of modules that should be executed.  All entries of the
		# list should contain an executable accepting only parameter file
		# or dict.
		'modules' : [],
		'processes' : 1
		}


def execute(pipe_file_or_dict, feedback=0) :
	"""Execute all the modules listed in the input file."""
	
	params, module_params = parse_ini.parse(pipe_file_or_dict, params_init, 
	                                   prefix='pipe_',return_undeclared=True,
	                                   feedback=feedback)
	
	#print rank
	for module in params['modules'] :
		if feedback > 1 :
			print 'Excuting analysis module: ' + str(module)
		module(module_params, feedback=feedback).execute(params['processes'])

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
	import sys
	execute(str(sys.argv[1]))

