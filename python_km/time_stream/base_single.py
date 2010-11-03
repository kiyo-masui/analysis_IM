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

from kiyopy import parse_ini, utils
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
               'output_end' : ".testout.fits",
               # What data to process within each file.
               'scans' : (),
               'IFs' : ()
               }

class BaseSingle(object) :
    """This is a base class from which some time stream steps should inherit.
    
    See module doc-string for detailed information.
    """
    
    # Here are the things that should be overwritten by a class inheriting
    # from this one.
    
    # Parameters additional to the above base_params
    params_init = {}
    # Prefixes to parameter names in input files.  Should be unique amoung
    # analysis modules.
    prefix = 'bs_'
    # Message to print after reading a file. Nothing printed if attribute is
    # missing.
    # feedback_title = 'What was done this Data Block: '
    # What will acctually be done. A function names action whose argument and
    # return value are each DataBlock objects.
    def action(self, Data) :
        Data.add_history('Did Nothing.')
        # Messages to print this iteration.
        self.block_feedback = 'Did nothing.'
        return Data
    # Don't worry too much about the feedback stuff; it's safe to ignore it.
    # If you want to see feedback in action loof at time_stream/flag_data.py.

    # Don't overwrite these in derivative classes.
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        self.feedback = feedback
        
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
        # Make parent directories if need be.
        utils.mkparents(params['output_root'])
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
            if hasattr(self, 'feedback_title') and self.feedback > 1:
                print self.feedback_title,
            # Get the number of scans if asked for all of them.
            scan_inds = params['scans']
            if len(scan_inds) == 0 or scan_inds is None :
                scan_inds = range(len(Reader.scan_set))
            # Loop over scans.
            for thisscan in scan_inds :
                Blocks = Reader.read(thisscan, params['IFs'])
                
                # Function that loops over DataBlocks within a scan.
                NewBlocks = self.scan_action(Blocks)
                del Blocks
                Writer.add_data(NewBlocks)
            
            # Go to a new line if we are printing statistics.
            if hasattr(self, 'feedback_title') and self.feedback > 1:
                    print ''
            # Finally write the data back to file.
            Writer.write(output_fname)
    
    # This part of the loop has been split off to make window stitching easier
    # since IF stitching operates on all the IFs in a scan at the same time.
    # IF stitching will overwrite this function instead of action.
    def scan_action(self, scan_blocks) :
        """scan_blocks is a tuple of DataBlocks all from the same scan."""
        
        out_data = ()
        for Data in scan_blocks :
            # Now process the data.
            NewData = self.action(Data)
            # Chances are NewData is just a reference to Data.  To avoid
            # confusion, delete the old reference.  This is just paranoia.
            del Data

            # If this module sets some feedback statistics, print the statistic
            # but don't go to a new line.
            if hasattr(self, 'block_feedback') and self.feedback > 1:
                print self.block_feedback,

            # Check that the data is valid and add it to the writer.
            NewData.verify()

            out_data = out_data + (NewData,)

        return out_data

