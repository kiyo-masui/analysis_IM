from kiyopy import parse_ini, utils
from core import data_block
from core import fitsGBT

# This is a list of input parameters that this module takes and thier default
# value.
params_init = {
               # IO:
               'input_root' : './',
               # The unique part of every fname
               'file_middles' : ("testfile_cal_name",),
               'input_end' : ".fits",
               'output_root' : "./testoutput",
               'output_end' : ".calT.fits"
               }

prefix = 'name_' # A two letter prefix identifying this module, followed by an
                 # underscore.

class Name(object) :
    """Doc-string."""
    
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        self.feedback = feedback
        
        # Read in the parameters.  Parameters can be passed as a dictionary or
        # a file name.  If the input is 'None' then all parameters revert to
        # the default.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init,
                                      prefix=prefix,
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
            
            first_scan_and_IF_DataBlock = Reader.read(scans=0,IFs=0)
            second_scan_and_first_IF_DataBlock = Reader.read(scans=1,IFs=0)
            # Kevin, after this you are on your own.




