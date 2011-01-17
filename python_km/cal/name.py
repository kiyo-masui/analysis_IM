"""This is the skeleton of a module and is here to demonstrate how a module
should be written."""

from kiyopy import parse_ini, utils
from core import data_block
from core import fitsGBT

# The first thing we do is define a list of input parameters that this module
# will read from file.
prefix = 'name_' # A two letter prefix identifying this module, followed by an
                 # underscore.  This is so several different modules can share
                 # the same input file.
# This is a list of input parameter names that this module takes and thier
# default value.
params_init = {
               # IO:
               'input_root' : './',
               # The unique part of every fname
               'file_middles' : ("testfile_cal_name",),
               'input_end' : ".fits",
               'output_root' : "./testoutput",
               'output_end' : ".calT.fits"
               }
# Now this module expects as it's only argument a file name where the above
# variables are defined (with the prefix) in python syntax.


# By convention all pipeline modules should be classes with at least 2 methods:
# __init__ and execute.  This convension makes it easier to script the module
# later.  Both these methods are called once when the module is executed.
class Name(object) :
    """Doc-string."""
    
    # The init method reads in parameters mostly, but you can put other things
    # here.
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        self.feedback = feedback
        
        # Read in the parameters.  Parameters can be passed as a dictionary or
        # a file name.  If the input is 'None' then all parameters revert to
        # the default.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init,
                                      prefix=prefix,
                                      checking=10*self.feedback + 2)
    
    # The execute method does all the work.  This is were you should write your
    # program.
    def execute(self) :
        """Process all data."""
        
        # You have access to the input parameters through the dictionary
        # self.params.
        params = self.params
        # If you have output files, make parent directories if need be.
        utils.mkparents(params['output_root'])
        # Write the input parameters to file so you can go back and look at
        # them.
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=self.prefix)

        # Loop over the files to process.
        for file_middle in params['file_middles'] :
            input_fname = (params['input_root'] + file_middle +
                           params['input_end'])
            
            # Read in the data.  The reader is an object that can read
            # DataBlock objects out of a fits file.
            Reader = fitsGBT.Reader(input_fname, feedback=self.feedback)
            
            # Some examples of how you would read DataBlock Objects:
            first_scan_and_IF_DataBlock = Reader.read(scans=0,IFs=0)
            second_scan_and_first_IF_DataBlock = Reader.read(scans=1,IFs=0)
            list_of_a_few_data_blocks = Reader.read(scans=(1,2,3),IFs=0)
            list_of_all_data_blocks = Reader.read(scans=(),IFs=())
            


# If this file is run from the command line, execute the main function.
if __name__=="__main__":
    import sys
    input_file_name = str(sys.argv[1])
    # Initialize the object.
    Object = Name(input_file_name)
    # Run the script.
    Object.execute()
