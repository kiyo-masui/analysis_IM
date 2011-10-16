"""DataPath class specification"""
import os
import sys
import time
import datetime
import urllib2
import subprocess
import getpass


def print_path_properties(pathname):
    r"""Print whether or not a path is writable
    """
    entry = "path %s:" % pathname

    if os.access(pathname, os.W_OK):
        entry += " is writable"
    else:
        entry += " is not writable"

    print entry


def file_properties(filename):
    r"""determine if a file exists, and if so get its date and permissions,
    if not then see if the directory is writable
    """
    exists = True
    try:
        open(filename)
    except IOError:
        exists = False

    if exists:
        writable = os.access(filename, os.W_OK)
        modtime = time.ctime(os.path.getmtime(filename))
    else:
        directory = "/".join(filename.split("/")[:-1])
        writable = os.access(directory, os.W_OK)
        modtime = False

    return (modtime, writable)


def print_file_properties(filename, prepend=""):
    r"""Print the output of file_properties as an entry
    with a `prepend` string, the filename, whether it is writable, and
    what its last modification date was if the file exists.
    """
    fileprop = file_properties(filename)
    entry = "file %s %s:" % (prepend, filename)

    if fileprop[1]:
        entry += " writable"
    else:
        entry += " is not writable"

    if fileprop[0]:
        entry += ", last modified: " + fileprop[0]
    else:
        entry += ", does not exist yet"

    print entry


def print_dictionary(dict_in, key_list=None, prepend=""):
    r"""print a dictionary in a slightly nicer way
    `key_list` specifies the keys to print; None is to print all
    `prepend` is a string to prepend to each entry (useful to tag)
    """
    if key_list == None:
        key_list = dict_in.keys()

    for key in key_list:
        print "%s%s: %s" % (prepend, key, dict_in[key])


class DataPath(object):
    r"""A class to manage the data path database

    'load_pathdict' in __init__ downloads and executes the group-writable
    database specification and it executes and extracts the path database and
    file groups. __init__ also prints the user's relevant github status
    information for logging (associating a run with a code version).

    The database specification format can be found in that file's doc string:
    /cita/d/www/home/eswitzer/GBT_param/path_database.py
    or
    http://www.cita.utoronto.ca/~eswitzer/GBT_param/path_database.py

    Usage:
    # start the class, printing version information
    >>> datapath_db = DataPath()
    DataPath: opening URL
        http://www.cita.utoronto.ca/~eswitzer/GBT_param/path_database.py
    DataPath: File database date/version ...
    DataPath: Run info: ... by ...
    git_SHA: ...
    git_blame: ...
    git_date: ...

    # pick index '44' of the 15hr sims
    >>> datapath_db.fetch("sim15hr_beam", pick='44')
    file ...simulations/15hr/sim_beam_044.npy: ...
    '...simulations/15hr/sim_beam_044.npy'

    # get the 15hr sim path
    >>> datapath_db.fetch("sim15hr_path")
    path ...simulations/15hr/: ...
    '...simulations/15hr/'

    TODO: also allow dbs from local paths instead of URLs
    TODO: switch to ordered dictionaries instead of list+dictionary?
    TODO: code check that all files in the db exist, etc.

    Extensions to consider:
        -require writing to a log file; check exists; opt. overwrite
        -functions support links through e.g.
            15hrmap -> 15hrmap_run5_e1231231231_etc_Etc.
        -framework for data over http
    """

    # URL to the default path database
    _db_root = "http://www.cita.utoronto.ca/~eswitzer/GBT_param/"
    _db_file = "path_database.py"
    _db_url_default = _db_root + _db_file

    def __init__(self, db_url=_db_url_default, skip_gitlog=False):
        r"""Load a file path specification and get basic run info
        """

        self._pathdict = {}     # the main path database
        self._groups = {}       # dictionary of database key lists by group
        self._group_order = []  # order in which to print the groups
        self.version = "Empty"  # the database version
        self.db_url = db_url    # URL to load the database specification from
        self.gitlog = "Empty"   # git SHA and commit info

        # load the file database and git version info
        self.load_pathdict(db_url)
        self.check_groups()
        self.clprint("File database date/version " + self.version)

        # get user and code version information
        dateinfo = datetime.datetime.now()
        self.runinfo = (dateinfo.strftime("%Y-%m-%d %H:%M:%S"),
                        getpass.getuser())
        self.clprint("Run info: %s by %s" % self.runinfo)
        if not skip_gitlog:
            self.get_gitlog()

    def clprint(self, string_in):
        r"""print with class message; could extend to logger"""
        print "DataPath: " + string_in

    def load_pathdict(self, db_url, print_dbspec=False):
        r"""Load the parameter dictionary

        note that the class instance update()'s its dictionary, so subsequent
        calls of this function on different files will overwrite or augment
        dictionaries that have already been loaded.
        """
        self.clprint("opening URL " + self.db_url)

        resp = urllib2.urlopen(db_url)
        # urllib2 will die if this fails, but to be safe,
        if (resp.code != 200):
            print "ERROR: path database URL invalid (%s)" % resp.code

        path_db_spec = resp.read()
        if print_dbspec:
            print "-" * 80
            self.clprint("executing: ")
            print path_db_spec
            print "-" * 80

        # evaluate the database specification inside a dictionary
        bounding_box = {}
        code = compile(path_db_spec, '<string>', 'exec')
        exec code in bounding_box

        # extract only the useful database information
        self.version = bounding_box['version_tag']
        self._pathdict.update(bounding_box['pdb'])
        self._groups.update(bounding_box['groups'])
        self._group_order = bounding_box['group_order']

    def check_groups(self):
        r"""check that all the database keys are accounted for by groups
        """
        keylist_in_groups = []
        keylist = self._pathdict.keys()

        for item in self._groups:
            keylist_in_groups.extend(self._groups[item])

        for groupkey in keylist_in_groups:
            if groupkey not in keylist:
                print "ERROR: " + groupkey + " is not in the database"
                sys.exit()

        for dbkey in keylist:
            if dbkey not in keylist_in_groups:
                print "ERROR: " + dbkey + " is not in the group list"
                sys.exit()

    def print_db_item(self, dbkey, suppress_lists=6, silent=False):
        r"""print a database entry to markdown format
        suppress printing of lists of more than `suppress_lists` files
        `silent` only returns a string instead of printing
        """
        dbentry = self._pathdict[dbkey]
        retstring = "### `%s`\n" % dbkey
        #retstring += "-"*len(retstring) + "\n"
        retstring += "* __Description__: %s\n" % dbentry['desc']

        if 'notes' in dbentry:
            retstring += "* __Notes__: %s\n" % dbentry['notes']

        if ('filelist' in dbentry):
            listindex = dbentry['listindex']
            retstring += "* __List index__: " + repr(listindex) + "\n"
            if (len(listindex) <= suppress_lists):
                retstring += "* __Filelist__:\n"
                for listitem in listindex:
                    retstring += "    * " + listitem + " :`" + \
                                 dbentry['filelist'][listitem] + "`\n"

        if 'path' in dbentry:
            retstring += "* __Path__: `%s`\n" % dbentry['path']

        if 'file' in dbentry:
            retstring += "* __File__: `%s`\n" % dbentry['file']

        if not silent:
            print retstring

        return retstring

    def print_path_db(self, suppress_lists=6, fileobj=None):
        r"""print all the files in the path database; note that it is a
        dictionary and so this is un-ordered. print_path_db_by_group prints the
        database items ordered by the file group specifications.

        Suppress printing of lists of more than `suppress_lists` files.
        If given a filename, this will write to a text file in markdown format.
        """
        print "-" * 80
        for dbkey in self._pathdict:
            dbstring = self.print_db_item(dbkey, suppress_lists=suppress_lists)

            if fileobj:
                fileobj.write(dbstring + "\n")

        print "-" * 80

    def print_path_db_by_group(self, suppress_lists=6, fileobj=None):
        r"""print all the files in the path database ordered by group
        specification.

        Suppress printing of lists of more than `suppress_lists` files.
        If given a filename, this will write to a text file in markdown format.
        """
        print "-" * 80
        for groupname in self._group_order:
            print "%s\n%s\n" % (groupname, "-" * len(groupname))
            if fileobj:
                fileobj.write("****\n %s\n%s\n\n" % \
                              (groupname, "-" * len(groupname)))

            for dbkey in self._groups[groupname]:
                dbstring = self.print_db_item(dbkey,
                                              suppress_lists=suppress_lists)
                if fileobj:
                    fileobj.write(dbstring + "\n")

        print "-" * 80

    def get_gitlog(self):
        r"""parse the github status for this code version for logging
        """
        process = subprocess.Popen(["git", "log"], stdout=subprocess.PIPE)
        gitlog = process.communicate()[0]
        gitlog = gitlog.split("\n")
        gitdict = {}
        gitdict['SHA'] = " ".join(gitlog[0].split()[1:])
        gitdict['blame'] = " ".join(gitlog[1].split()[1:])
        gitdict['date'] = " ".join(gitlog[2].split()[1:])
        gitdict['note'] = " ".join(gitlog[4].split()[1:])
        self.gitlog = gitdict

        print_dictionary(self.gitlog,
                         key_list=["SHA", "blame", "date"],
                         prepend="git_")

    def fetch(self, dbkey, pick=None):
        r"""The access function for this database class:
        Fetch the data for a requested key in the db.

        `pick` takes one index from a file list
        if the database entry is a file list, return a tuple of the index
        indices and a dictionary (the list orders the dictionary)
        """
        dbentry = self._pathdict[dbkey]

        if 'file' in dbentry:
            pathout = dbentry['file']
            print_file_properties(pathout)

        if 'path' in dbentry:
            pathout = dbentry['path']
            print_path_properties(pathout)

        if 'filelist' in dbentry:
            pathout = (dbentry['listindex'], dbentry['filelist'])
            if pick:
                pathout = pathout[1][pick]
                print_file_properties(pathout)
            else:
                for item in pathout[0]:
                    filepath = pathout[1][item]
                    print_file_properties(filepath, prepend=item)

        return pathout

    def generate_path_webpage(self):
        r"""Write out a markdown file with the file database
        """
        localpath = "/cita/d/www/home/eswitzer/GBT_param/"
        localdb = localpath + "path_database.py"
        dbwebpage = localpath + "path_database.txt"

        print "writing markdown website to: " + dbwebpage
        fileobj = open(dbwebpage, "w")

        fileobj.write("Data path DB\n============\n\n")
        fileobj.write("specified by local path: `%s` and URL `%s`\n" %
                      (localdb, self.db_url))
        self.print_path_db_by_group(suppress_lists=6, fileobj=fileobj)

        fileobj.close()


if __name__ == "__main__":
    import doctest

    # generate the path db markdown website
    datapath_db = DataPath()
    datapath_db.generate_path_webpage()

    # run some tests
    OPTIONFLAGS = (doctest.ELLIPSIS |
                   doctest.NORMALIZE_WHITESPACE)
    doctest.testmod(optionflags=OPTIONFLAGS)
