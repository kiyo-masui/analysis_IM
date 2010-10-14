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
              }


def execute(pipe_file_or_dict) :
    """Execute all the modules listed in the input file."""

    params, module_params = parse_ini.parse(pipe_file_or_dict, params_init, 
                                       prefix='pipe_',return_undeclared=True)
    
    for module in params['modules'] :
        print 'Excuting analysis module: ' + str(module)
        module(module_params).execute()

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    execute(str(sys.argv[1]))

