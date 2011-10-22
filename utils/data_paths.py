"""DataPath class specification"""
import os
import sys
import time
import datetime
import urllib2
import subprocess
import getpass
import hashlib


def path_properties(pathname, intend_write=False, intend_read=False,
                    file_index="", prefix="", silent=False, is_file=False,
                    die_overwrite=False):
    r"""Check various properties about the path, print the result
    if `intend_write` then die if the path does not exist or is not writable
    if `intend_read` then die if the path does not exist
    `file_index` is an optional mark to put in front of the file name
    `prefix` is an optional mark to prefix the file property output
    `silent` turns off printing unless there was an intend_r/w error
    `is_file` designates a file; if it does not exist check if the path is
    writable.

    # Usage examples:
    >>> path_properties("this_does_not_exist.txt", is_file=True)
    => this_does_not_exist.txt: does not exist,
        but path: ./ exists, is writable, is not readable.
    (True, False, True, False, '...')

    >>> path_properties("/tmp/nor_does_this.txt", is_file=True)
    => /tmp/nor_does_this.txt: does not exist,
        but path: /tmp/ exists, is writable, is not readable.
    (True, False, True, False, '...')

    # here's one that you should not be able to write to
    >>> path_properties("/var/nor_does_this.txt", is_file=True)
    => /var/nor_does_this.txt: does not exist,
        but path: /var/ exists, is not writable, is not readable.
    (True, False, False, False, '...')

    # here's one that definitely exists
    >>> path_properties(__file__, is_file=True)
    => ...: exists, is writable, is readable.
    (True, True, True, False, '...')

    # test a writeable directory that exists
    >>> path_properties('/tmp/')
    => /tmp/: exists, is writable, is readable.
    (True, True, True, True, '...')
    """
    entry = "%s=> %s%s:" % (prefix, file_index, pathname)

    exists = os.access(pathname, os.F_OK)
    # save the state of the file (exist will later specify if the parent dir
    # exists, but we also care about this specific file)
    file_exists = exists
    readable = os.access(pathname, os.R_OK)
    writable = os.access(pathname, os.W_OK)
    executable = os.access(pathname, os.X_OK)
    if exists:
        modtime = time.ctime(os.path.getmtime(pathname))
    else:
        modtime = None

    # if this is a file that does not exist, check the directory
    if not exists and is_file:
        directory = "/".join(pathname.split("/")[:-1])
        if directory == "":
            directory = "."
        directory += "/"

        writable = os.access(directory, os.W_OK)
        exists = os.access(directory, os.F_OK)
        if exists:
            modtime = time.ctime(os.path.getmtime(directory))

        readable = False
        entry += " does not exist, but path: %s" % directory

    if exists:
        entry += " exists,"

        if writable:
            entry += " is writable,"
        else:
            entry += " is not writable,"

        if readable:
            entry += " is readable."
        else:
            entry += " is not readable."
    else:
        entry += " does not exist"

    if not silent:
        print entry

    if intend_read and not readable:
        print "ERROR: " + entry
        sys.exit()

    if intend_write and not writable:
        print "ERROR: " + entry
        sys.exit()

    if intend_write and file_exists:
        if die_overwrite:
            print "OVERWRITE ERROR: " + entry
            sys.exit()
        else:
            print "WARNING: you will overwrite " + entry

    return (exists, readable, writable, executable, modtime)


def hashfile(filename, hasher=hashlib.md5(), blocksize=65536, max_size=1.e8):
    r"""determine the hash of a file in blocks
    if it exceeds `max_size` just report the modification time
    if the file does not exist, report '-1'
    """
    if os.access(filename, os.F_OK):
        if (os.path.getsize(filename) < max_size):
            afile = open(filename, 'r')
            buf = afile.read(blocksize)
            while len(buf) > 0:
                hasher.update(buf)
                buf = afile.read(blocksize)

            afile.close()
            hash = hasher.hexdigest()
        else:
            hash = time.ctime(os.path.getmtime(filename))
    else:
        hash = "-1"

    return hash


def print_dictionary(dict_in, handle, key_list=None, prepend=""):
    r"""print a dictionary in a slightly nicer way
    `key_list` specifies the keys to print; None is to print all
    `prepend` is a string to prepend to each entry (useful to tag)
    """
    if key_list == None:
        key_list = dict_in.keys()

    for key in key_list:
        print >> handle, "%s%s: %s" % (prepend, key, dict_in[key])


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
    >>> datapath_db.fetch("sim_15hr_beam", pick='44')
    (sim_15hr_beam) => .../simulations/15hr/sim_beam_044.npy: ...
    '.../simulations/15hr/sim_beam_044.npy'

    # get the 15hr sim path
    >>> datapath_db.fetch("sim_15hr_path")
    (sim_15hr_path) => .../simulations/15hr/: ...
    '.../simulations/15hr/'

    TODO: also allow dbs from local paths instead of URLs
    TODO: switch to ordered dictionaries instead of list+dictionary?
    TODO: code check that all files in the db exist, etc.
    TODO: find db size in memory and total # files, print on website

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

    def print_db_item(self, dbkey, suppress_lists=6, silent=False,
                      skiphash=False):
        r"""print a database entry to markdown format
        'desc' and 'status' are requires keys in the file attributes
        suppress printing of lists of more than `suppress_lists` files
        `silent` only returns a string instead of printing
        """
        hashdict = {}  # calculate and record the hashes of the files
        dbentry = self._pathdict[dbkey]
        retstring = "### `%s`\n" % dbkey
        retstring += "* __Description__: %s\n" % dbentry['desc']

        if 'status' in dbentry:
            retstring += "* __Status__: %s\n" % dbentry['status']

        if 'notes' in dbentry:
            retstring += "* __Notes__: %s\n" % dbentry['notes']

        if ('filelist' in dbentry):
            listindex = dbentry['listindex']
            retstring += "* __List index__: `" + repr(listindex) + "`\n"

            if not skiphash:
                for listitem in listindex:
                    filename = dbentry['filelist'][listitem]
                    hashdict[filename] = hashfile(filename)

            if (len(listindex) <= suppress_lists):
                retstring += "* __File list__:\n"
                for listitem in listindex:
                    filename = dbentry['filelist'][listitem]
                    if not skiphash:
                        retstring += "    * `%s`:  (md5: %s) `%s`\n" % \
                                    (listitem, hashdict[filename], filename)
                    else:
                        retstring += "    * `%s`: `%s`\n" % \
                                    (listitem, filename)

        if 'path' in dbentry:
            retstring += "* __Path__: `%s`\n" % dbentry['path']

        if 'file' in dbentry:
            filename = dbentry['file']
            if not skiphash:
                hashdict[filename] = hashfile(filename)
                retstring += "* __File__: (md5: %s) `%s`\n" % \
                             (hashdict[filename], dbentry['file'])
            else:
                retstring += "* __File__: `%s`\n" % dbentry['file']

        if not silent:
            print retstring

        return (retstring, hashdict)

    def print_path_db(self, suppress_lists=6, fileobj=None, skiphash=False):
        r"""print all the files in the path database; note that it is a
        dictionary and so this is un-ordered. print_path_db_by_group prints the
        database items ordered by the file group specifications.

        You should only need to use this if the db groups are broken.

        Suppress printing of lists of more than `suppress_lists` files.
        If given a filename, this will write to a text file in markdown format.
        """
        print "-" * 80

        for dbkey in self._pathdict:
            (dbstring, hashdict) = self.print_db_item(dbkey,
                                          suppress_lists=suppress_lists,
                                          skiphash=skiphash)
            print dbstring

        print "-" * 80

    def print_path_db_by_group(self, suppress_lists=6, fileobj=None,
                               hashobj=None):
        r"""print all the files in the path database ordered by group
        specification.

        Suppress printing of lists of more than `suppress_lists` files.
        If given a filename, this will write to a text file in markdown format.
        """
        print "-" * 80

        if hashobj:
            skiphash = False
        else:
            skiphash = True

        for groupname in self._group_order:
            print "%s\n%s\n" % (groupname, "-" * len(groupname))
            if fileobj:
                fileobj.write("****\n %s\n%s\n\n" % \
                              (groupname, "-" * len(groupname)))

            for dbkey in self._groups[groupname]:
                (dbstring, hashdict) = self.print_db_item(dbkey,
                                              suppress_lists=suppress_lists,
                                              skiphash=skiphash)
                if fileobj:
                    fileobj.write(dbstring + "\n")

                if hashobj:
                    print_dictionary(hashdict, hashobj)

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

        print_dictionary(self.gitlog, sys.stdout,
                         key_list=["SHA", "blame", "date"],
                         prepend="git_")


    def fetch(self, dbkey, pick=None, intend_read=False, intend_write=False,
              purpose="", silent=False):
        r"""The access function for this database class:
        Fetch the data for a requested key in the db.

        `pick` takes one index from a file list
        if the database entry is a file list, return a tuple of the index
        indices and a dictionary (the list orders the dictionary)

        if `intend_write` then die if the path does not exist or not writable
        if `intend_read` then die if the path does not exist
        `purpose` inputs a purpose for this file for logging
        `silent` does not print anything upon fetch unless error
        """
        dbentry = self._pathdict[dbkey]
        prefix = "%s (%s) " % (purpose, dbkey)

        if 'file' in dbentry:
            pathout = dbentry['file']
            path_properties(pathout, intend_write=intend_write,
                            intend_read=intend_read, is_file=True,
                            prefix=prefix, silent=silent)

        if 'path' in dbentry:
            pathout = dbentry['path']
            path_properties(pathout, intend_write=intend_write,
                            intend_read=intend_read, is_file=False,
                            prefix=prefix, silent=silent)

        if 'filelist' in dbentry:
            pathout = (dbentry['listindex'], dbentry['filelist'])
            if pick:
                pathout = pathout[1][pick]
                path_properties(pathout, intend_write=intend_write,
                                intend_read=intend_read, is_file=True,
                                prefix=prefix, silent=silent)
            else:
                for item in pathout[0]:
                    filepath = pathout[1][item]
                    path_properties(filepath, intend_write=intend_write,
                                    intend_read=intend_read, is_file=True,
                                    file_index=item, prefix=prefix,
                                    silent=silent)

        return pathout

    def generate_path_webpage(self, skiphash=False):
        r"""Write out a markdown file with the file database
        """
        localpath = "/cita/d/www/home/eswitzer/GBT_param/"
        localdb = localpath + "path_database.py"
        dbwebpage = localpath + "path_database.txt"
        dbhashpage = localpath + "hashlist.txt"

        print "writing markdown website to: " + dbwebpage
        print "writing file hash list to: " + dbhashpage
        fileobj = open(dbwebpage, "w")
        if skiphash:
            hashobj = None
        else:
            hashobj = open(dbhashpage, "w")

        fileobj.write("Data path DB\n============\n\n")
        fileobj.write("specified by local path: `%s` and URL `%s`\n" %
                      (localdb, self.db_url))
        self.print_path_db_by_group(suppress_lists=30, fileobj=fileobj,
                                    hashobj=hashobj)

        fileobj.close()
        if not skiphash:
            hashobj.close()


if __name__ == "__main__":
    import doctest

    # generate the path db markdown website
    datapath_db = DataPath()
    datapath_db.generate_path_webpage()
    #datapath_db.generate_path_webpage(skiphash=True)

    # run some tests
    OPTIONFLAGS = (doctest.ELLIPSIS |
                   doctest.NORMALIZE_WHITESPACE)
    doctest.testmod(optionflags=OPTIONFLAGS)
