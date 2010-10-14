"""The base class for time stream steps that deal with single DataBlocks.

This class should not be used directly.  Instead any time stream processing
step where only a single DataBlock is processed at a time can be quickly
wirtten by inheriting from this class.  This class knows how to loop over files
and DataBlocks, read and write fits files, etc.

To create a valid time stream step, simply create a class that inherits from
this one and over-write the 'actions' method a the function that acctually acts
on your data.  The function must accept and retrun a DataBlock object.  These
functions will be applied to every DataBlock processed.  Optionally, overwrite
the params_init attribute with any additional options required by your
functions.  You should also overwrite the prefix attribute to something unique
so that parameters don't get confused in parameter files.

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
               # The unique part of every fname
               'file_middles' : ("testfile_GBTfits",),
               'input_end' : ".fits",
               'output_root' : "./",
               'output_end' : ".fits",
               # What data to process within each file.
               'scans' : (),
               'IFs' : ()
               }

class BaseSingle(object) :
    """This is a base class from which some time stream steps should inherit.
    
    See module doc-string for detailed information.
    """
    
    params_init = {}
    feedback = 2
    prefix = 'bs_'
    def action(self, Data) :
        Data.add_history('Did Nothing.')
        return Data

    def __init__(self, parameter_file_or_dict) :
        
        # Add any parameters added by any class that inherits from this one.
        params_init = dict(base_params)
        params_init.update(self.params_init)
        
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init,
                                      prefix=self.prefix,
                                      checking=10*self.feedback + 2)
    
    def execute(self) :
        """Process all data."""

        params = self.params
        #mkdir_p(params['output_dir*'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
        prefix=self.prefix)

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
                NewData = self.action(Data)
                # Chances are NewData is just a reference to Data.  To avoid
                # confusion, delete the old reference.
                del Data

                # Check that the data is valid and add it to the writer.
                NewData.verify()
                Writer.add_data(NewData)

            # Finally write the data back to file.
            Writer.write(output_fname)

        
