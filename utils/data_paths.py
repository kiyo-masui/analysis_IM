"""
Structure of the path specifications

path: full path to file
desc: description of file and code used to produce it
notes: notes about how the file was generated
time: auto-generated timestamp

different sorts of data:
single files
groups of files with the same desc, etc. (specified by tags in a dict)
directories for outputs

wildcards:
directories specific to the user
directories specific to the machine
"""
import os.path, time

# grab the file from the website, or optionally grab another path
# execute the file
# check sanity:
# check that all of the files exist, warning if it is not writeable
# any exceptions, environment variables not set?

# get the time stamp of the file specifying the paths, print it
# get the github sha for this version of the code
# shelve the whole object into the requested log file

# TODO: make this a class?
# TODO: recognize URLs as links
# extend the dict class, make it readonly

# functions support links through e.g. 15hrmap -> 15hrmap_run5_e1231231231_etc_Etc.

def load_pathdict():
    """Load the parameter directionry
    """

def generate_path_webpage():
    """
    """

def load_config(logfilename, filename=None):
    """Load a file path configuration
    logfilename specifies the file to save the path dictionary to;
    it is a required argument to try to force logging of the path parameters
    for each data run.
    """
    if filename:
        pathdict = 
    else:
        pathdict = 

    print "config file last modified: %s" % time.ctime(os.path.getmtime(file))

