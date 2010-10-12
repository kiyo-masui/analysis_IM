"""The base class for time stream steps that deal with single DataBlocks.

This class should not be used directly.  Instead any time stream processing
step where only a single DataBlock is processed at a time can be quickly
wirtten by inheriting from this class.  This class knows how to loop over files
and DataBlocks, read and write fits files, etc.

To create a valid time stream step, simply create a class that inherits from
this one and over-write the 'actions' attribute with a tuple of references to
functions.  These functions must accept DataBlock objects as arguments.  These
functions will be applied to every DataBlock processed.  Optionally, overwrite
the params_init attribute with any additional options required by your
functions.

Check out the unit tests for examples.
"""

from kiyopy import parse_ini
from kiyopy.utils import mkdir_p
#from kiyopy import custom_exceptions as ce

from core import data_block
from core import fitsGBT

# Parameters relevant for all modules that inherit from this one.
base_params = {
               # IO:
               'input_root' : './',
               'file_middles' : ["GBTdata"], # The unique part of every fname
               'input_end' : ".raw.acs.fits",
               'output_root' : "./",
               'output_end' : ".fits",
               # What data to process within each file.
               'scans' : [],
               'IFs' : []
               }

class BaseSingle(object) :
    """This is a base class from which some time stream steps should inherit.
    
    See module doc-string for detailed information.
    """
    
    actions = ()
    params_init = {}
    feedback = 2

    def __init__(self, parameter_file_or_dict) :
        
        # Add any parameters added by any class that inherited from this one.
        params_init = dict(base_params)
        params_init.update(self.params_init)
        
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init,
                                      checking=10*self.feedback + 2)
        params = self.params
        #mkdir_p(params['output_dir*'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini')

        # Loop over the files to process.
        for file_middle in params['file_middles'] :
            input_fname = (params['input_root'] + file_middle +
                           params['input_end'])
            output_fname = (params['output_root']
                            + file_middle + params['output_end'])
            Writer = fitsGBT.Writer(feedback=self.feedback)
            
            # Read in the data, and loop over data blocks.
            Reader = fitsGBT.Reader(input_fname, feedback=self.feedback)
            Blocks = Reader.read(params['scans'], params['IFs'])
            for Data in Blocks :
                # Now process the data.
                for function in self.actions :
                    function(self, Data)

                # Check that the data is valid and add it to the writer.
                Data.verify()
                Writer.add_data(Data)

            # Finally write the data back to file.
            Writer.write(output_fname)

def merge_singles(*args) :
    """Merge several time stream steps that inherit from BaseSingle.

    This takes a list of classes as an argument.  All classes must inherit from
    BaseSingle.  A new class is returned that combines all the steps in the
    passed classes.
    """

    new_actions = ()
    new_params_init = {}
    new_feedback = 0

    for SingleStep in args :
        new_actions = new_actions + Steps.actions
        new_params_ini.update(Steps.params_init)
        if Steps.feedback > new_feedback :
            new_feedback = Steps.feedback

    class NewSingle(BaseSingle) :
        pass

        
