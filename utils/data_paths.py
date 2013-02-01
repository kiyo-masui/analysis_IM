"""DataPath class specification"""
import sys
import copy
import time
import datetime
import urllib2
import subprocess
import getpass
import ast
import re
import os
import copy
from utils import file_tools as ft
from core import algebra


def print_dictionary(dict_in, handle, key_list=None, prepend=""):
    r"""print a dictionary in a slightly nicer way
    `key_list` specifies the keys to print; None is to print all
    `prepend` is a string to prepend to each entry (useful to tag)
    """
    if key_list == None:
        key_list = dict_in.keys()

    for key in key_list:
        print >> handle, "%s%s: %s" % (prepend, key, repr(dict_in[key]))


def unique_list(listin):
    seen = set()
    seen_add = seen.add
    uniq = [ x for x in listin if x not in seen and not seen_add(x)]
    return sorted(uniq)


def tack_on_subdir(filestring, tack_on):
    r"""Given a full path to a file, tack on a new subdir:
    >>> tack_on_subdir("/path/to/file/ok.dat", "tack_on")
    '/path/to/file/tack_on/ok.dat'
    >>> tack_on_subdir("/path/to/file/", "tack_on")
    '/path/to/file/tack_on/'
    >>> tack_on_subdir("./", "tack_on")
    './tack_on/'
    >>> tack_on_subdir("ok.dat", "tack_on")
    './tack_on/ok.dat'
    """
    splitfile = list(os.path.split(filestring))
    # if the path is empty, recognize it as the local directory
    if splitfile[0] == "":
        splitfile[0] = "."

    return "%s/%s/%s" % (splitfile[0], tack_on, splitfile[1])

def unpack_cases(case_list, case_key, divider=";"):
    r"""Given a list of cases like:
    >>> case_key = "pair;map_type;treatment"
    >>> case_list = ['A;1;a', 'A;1;b', 'A;2;a', 'A;2;b']
    >>> unpack_cases(case_list, case_key)
    {'pair': ['A'], 'map_type': ['1', '2'], 'treatment': ['a', 'b']}

    """
    case_keys = case_key.split(divider)
    # python 2.7 syntax; broken on tpb, use 2.6
    #case_counter = { case_key: [] for case_key in case_keys}
    case_counter = {}
    for case_key in case_keys:
        case_counter[case_key] = []

    for entry in case_list:
        csplit = entry.split(divider)
        if len(csplit) == len(case_keys):
            for (ckey, cval) in zip(case_keys, csplit):
                case_counter[ckey].append(cval)

    for ckey in case_counter:
        case_counter[ckey] = unique_list(case_counter[ckey])

    return case_counter


def convert_dbkeydict_to_filedict(dbkeydict, datapath_db=None, silent=True):
    r"""simple caller to convert a dictionary of database keys to a dictionary
    of filenames"""
    if datapath_db is None:
        datapath_db = DataPath()

    filedict = {}
    for name in dbkeydict:
        filedict[name] = datapath_db.fetch(dbkeydict[name], silent=silent)

    return filedict


def extract_split_tag(keylist, divider=";", ignore=None):
    r"""take a list like ['A;ok', 'B;ok'], and return ['A', 'B']
    list is made unique and sorted alphabetically
    anything in `ignore` is thrown out
    TODO: use unique_list from above and possible merge with unpack_cases
    """
    taglist = []
    for key in keylist:
        taglist.append(key.split(divider)[0])

    taglist = list(set(taglist))

    if ignore is not None:
        taglist = [tag for tag in taglist if tag not in ignore]

    taglist.sort()
    return taglist


def GBTauto_cross_pairs(list1, list2, cross_sym="_with_"):
    r"""given the list of map tags, produce a list of paired tuples

    list1 = ['A_with_B', 'A_with_C', ..., 'B_with_A', ...]
    list2 = ['A_with_B', 'A_with_C', ..., 'B_with_A', ...]
    gives:
    [('A_with_B', 'B_with_A'), ('A_with_C', 'C_with_A')
     ('A_with_D', 'D_with_A'), ('B_with_C', 'C_with_B')
     ('B_with_D', 'D_with_B'), ('C_with_D', 'D_with_C')]
    """
    maplist1 = []
    for entry in list1:
        maplist1.extend(entry.split(cross_sym))
    maplist1 = list(set(maplist1))
    maplist1.sort()

    maplist2 = []
    for entry in list1:
        maplist2.extend(entry.split(cross_sym))
    maplist2 = list(set(maplist2))
    maplist2.sort()

    if (maplist1 != maplist2):
        print "ERROR: auto-power has unmatched map split"
        sys.exit()

    combolist = unique_cross_pairs(maplist1, maplist2)
    retlist = []
    for item in combolist:
        left = "%s%s%s" % (item[0], cross_sym, item[1])
        right = "%s%s%s" % (item[1], cross_sym, item[0])
        retlist.append((left, right))

    return retlist


def unique_cross_pairs(list1, list2):
    r"""given the list of map tags, produce a list of tuples of unique pairs
    of those map (recognizing AxB = BxA)

    >>> list1 = ['A', 'B', 'C', 'D']
    >>> list2 = ['A', 'B', 'C', 'D']
    >>> unique_cross_pairs(list1, list2)
    [('A', 'B'), ('A', 'C'), ('A', 'D'), ('B', 'C'), ('B', 'D'), ('C', 'D')]

    """
    retlist = []
    for ind1 in range(len(list1)):
        for ind2 in range(len(list2)):
            if ind2 > ind1:
                retlist.append((list1[ind1], list2[ind2]))

    return retlist

# TODO: write doctest for this
def cross_maps(map1key, map2key, noise_inv1key, noise_inv2key,
               map_suffix=";clean_map", noise_inv_suffix=";noise_inv",
               verbose=True, ignore=['firstpass'], cross_sym="_with_",
               pair_former="unique_cross_pairs",
               tag1prefix="", tag2prefix="", db_to_use=None):
    r"""Use the database to report all unique crossed map maps given map and
    noise_inv keys.
    """
    if db_to_use is None:
        db_to_use = DataPath()

    retpairs = {}
    retpairslist = []

    if os.path.exists(map1key):
        (map1keys, map1set) = get_mapdict(map1key)
        (noise_inv1keys, noise_inv1set) = get_mapdict(noise_inv1key)
    else:
        (map1keys, map1set) = db_to_use.fetch(map1key, intend_read=True,
                                        silent=True)
        (noise_inv1keys, noise_inv1set) = db_to_use.fetch(noise_inv1key,
                                        intend_read=True,
                                        silent=True)

    if os.path.exists(map2key):
        (map2keys, map2set) = get_mapdict(map2key)
        (noise_inv2keys, noise_inv2set) = get_mapdict(noise_inv2key)
    else:
        (map2keys, map2set) = db_to_use.fetch(map2key, intend_read=True,
                                        silent=True)
        (noise_inv2keys, noise_inv2set) = db_to_use.fetch(noise_inv2key,
                                        intend_read=True,
                                        silent=True)

    map1tags = extract_split_tag(map1keys, ignore=ignore)
    map2tags = extract_split_tag(map1keys, ignore=ignore)
    noise_inv1tags = extract_split_tag(noise_inv1keys, ignore=ignore)
    noise_inv2tags = extract_split_tag(noise_inv2keys, ignore=ignore)

    if (map1tags != noise_inv1tags) or (map2tags != noise_inv2tags):
        print "ERROR: noise_inv and maps are not matched"
        sys.exit()

    if verbose:
        print "Using map1 tags %s and map2 tags %s" % (map1tags, map2tags)

    dispatch = {'unique_cross_pairs': unique_cross_pairs,
                'GBTauto_cross_pairs': GBTauto_cross_pairs}

    if pair_former not in dispatch:
        print "ERROR: pair_former is does not match a known function"
        sys.exit()

    comblist = dispatch[pair_former](map1tags, map2tags)

    for (tag1, tag2) in comblist:
        pairname = "%s%s%s%s%s" % \
                    (tag1prefix, tag1, cross_sym, tag2prefix, tag2)
        retpairslist.append(pairname)
        pairdict = {'map1': map1set[tag1 + map_suffix],
                    'noise_inv1': noise_inv1set[tag1 + noise_inv_suffix],
                    'map2': map2set[tag2 + map_suffix],
                    'noise_inv2': noise_inv2set[tag2 + noise_inv_suffix],
                    'tag1': tag1, 'tag2': tag2}

        if verbose:
            print "-"*80
            print_dictionary(pairdict, sys.stdout,
                            key_list=['map1', 'noise_inv1',
                                      'map2', 'noise_inv2',
                                      'tag1', 'tag2'])

        retpairs[pairname] = pairdict

    return (retpairslist, retpairs)


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
        http://www.cita.utoronto.ca/~eswitzer/GBT_param/path_database.shelve
    DataPath: opening URL
        http://www.cita.utoronto.ca/~eswitzer/GBT_param/hashlist.txt
    # checksums compiled: ... by ...
    DataPath: File database date/version ...
    DataPath: Run info: ... by ...
    git_SHA: ...
    git_blame: ...
    git_date: ...
    DataPath: ... files registered; database size in memory ...

    # pick index '44' of the 15hr sims
    #>>> datapath_db.fetch("sim_15hr_oldmap_ideal", pick='44')
    #(sim_15hr_beam) => .../simulations/15hr/sim_beam_044.npy: ...
    #'.../simulations/15hr/sim_beam_044.npy'

    # get the 15hr sim path
    #>>> datapath_db.fetch("sim_15hr_oldmap_ideal_path")
    #(sim_15hr_path) => .../simulations/15hr/: ...
    #'.../simulations/15hr/'

    TODO: also allow dbs from local paths instead of URLs
    TODO: switch to ordered dictionaries instead of list+dictionary?
    TODO: code check that all files in the db exist, etc.
    TODO: make command line to regen file hash, web page
    TODO: make fetch handle file hashes and note updates

    Extensions to consider:
        -require writing to a log file; check exists; opt. overwrite
        -functions support links through e.g.
            15hrmap -> 15hrmap_run5_e1231231231_etc_Etc.
        -framework for data over http
    """

    # URL to the default path database
    _db_root = "http://www.cita.utoronto.ca/~eswitzer/GBT_param/"
    #_db_file = "path_database.py"
    _db_file = "path_database.shelve"
    _db_url_default = _db_root + _db_file
    _hash_file = "hashlist.txt"
    _hash_url_default = _db_root + _hash_file

    def __init__(self, db_url=_db_url_default, hash_url=_hash_url_default,
                 skip_gitlog=False, load_via_exec=False):
        r"""Load a file path specification and get basic run info
        """

        self._pathdict = {}       # the main path database
        self._hashdict = {}       # file hash database
        self._groups = {}         # dictionary of database key lists by group
        self._group_order = []    # order in which to print the groups
        self.version = "Empty"    # the database version
        self.db_url = db_url      # URL: database specification
        self.hash_url = hash_url  # URL: file hash specification
        self.gitlog = "Empty"     # git SHA and commit info

        # load the file database and git version info
        if load_via_exec:
            self.execload_pathdict(db_url)
        else:
            self.load_pathdict(db_url)

        self.replace_local_variables()
        self.load_hashdict(hash_url)
        self.check_groups()
        self.clprint("File database date/version " + self.version)

        # get user and code version information
        dateinfo = datetime.datetime.now()
        self.runinfo = (dateinfo.strftime("%Y-%m-%d %H:%M:%S"),
                        getpass.getuser())
        self.clprint("Run info: %s by %s" % self.runinfo)
        if not skip_gitlog:
            self.get_gitlog()
        self._db_size = (len(self._hashdict),
                        sys.getsizeof(self._pathdict) +
                        sys.getsizeof(self._hashdict))

        self.clprint("%d files registered; database size in memory = %s" %
                     self._db_size)

    def clprint(self, string_in):
        r"""print with class message; could extend to logger"""
        print "DataPath: " + string_in

    def load_pathdict(self, db_url):
        r"""Load the parameter dictionary as a shelve through http

        note that the class instance update()'s its dictionary, so subsequent
        calls of this function on different files will overwrite or augment
        dictionaries that have already been loaded.
        """
        self.clprint("opening URL " + self.db_url)

        resp = urllib2.urlopen(db_url)
        # urllib2 will die if this fails, but to be safe,
        if (resp.code != 200):
            print "ERROR: path database URL invalid (%s)" % resp.code

        bounding_box = ft.load_shelve_over_http(db_url)

        # extract only the useful database information
        self.version = bounding_box['version_tag']
        self._pathdict.update(bounding_box['pdb'])
        self._groups.update(bounding_box['groups'])
        self._group_order = bounding_box['group_order']

    def execload_pathdict(self, db_url, print_dbspec=False):
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

    def replace_local_variables(self):
        r"""go through the entries and replace the enviornment variables with
        their actual values
        """
        for groupname in self._group_order:
            for db_key in self._groups[groupname]:
                dbentry = self._pathdict[db_key]

                if 'filelist' in dbentry:
                    listindex = dbentry['listindex']
                    for listitem in listindex:
                        filename = dbentry['filelist'][listitem]
                        filename = ft.repl_bracketed_env(filename)
                        dbentry['filelist'][listitem] = \
                                    re.sub('/+', '/', filename)

                if 'file' in dbentry:
                    filename = ft.repl_bracketed_env(dbentry['file'])
                    dbentry['file'] = re.sub('/+', '/', filename)

                if 'path' in dbentry:
                    pathname = ft.repl_bracketed_env(dbentry['path'])
                    dbentry['path'] = re.sub('/+', '/', pathname)

                self._pathdict[db_key] = dbentry

    def load_hashdict(self, hash_url):
        r"""Load the file hash library
        """
        self.clprint("opening URL " + self.hash_url)

        resp = urllib2.urlopen(hash_url)
        if (resp.code != 200):
            print "ERROR: path database URL invalid (%s)" % resp.code
        hash_spec = resp.read()
        hash_spec = hash_spec.split("\n")
        for entry in hash_spec:
            if len(entry) > 0:
                if (entry[0] == "#"):
                    print entry
                else:
                    parser = entry.split(": ")
                    filename = parser[0]
                    filename = filename.strip()
                    fileinfo = ast.literal_eval(parser[1])
                    self._hashdict[filename] = fileinfo

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

        for db_key in keylist:
            if db_key not in keylist_in_groups:
                print "ERROR: " + db_key + " is not in the group list"
                sys.exit()

    def print_db_item_desc(self, db_key):
        r"""print just the key: desc as a short summary"""
        dbentry = self._pathdict[db_key]
        return "* `%s`: %s\n" % (db_key, dbentry['desc'])

    def print_db_item(self, db_key, suppress_lists=90, silent=False):
        r"""print a database entry to markdown format
        'desc' and 'status' are requires keys in the file attributes
        suppress printing of lists of more than `suppress_lists` files
        `silent` only returns a string instead of printing
        """
        dbentry = self._pathdict[db_key]
        retstring = "### `%s`\n" % db_key
        retstring += "* __Description__: %s\n" % dbentry['desc']

        if 'parent' in dbentry:
            retstring += "* __Parent path key__: `%s`\n" % dbentry['parent']
            retstring += "* __Parent path__: `%s`\n" % \
                        self.fetch(dbentry['parent'], silent=True)

        if 'group_key' in dbentry:
            retstring += "* __Group key__: `%s`\n" % dbentry['group_key']

        if 'status' in dbentry:
            retstring += "* __Status__: %s\n" % dbentry['status']

        if 'notes' in dbentry:
            retstring += "* __Notes__: %s\n" % dbentry['notes']

        if 'filelist' in dbentry:
            listindex = dbentry['listindex']
            retstring += "* __List index__: `" + repr(listindex) + "`\n"

            if (len(listindex) <= suppress_lists):
                retstring += "* __File list__:\n"
                for listitem in listindex:
                    filename = dbentry['filelist'][listitem]
                    if filename in self._hashdict:
                        retstring += "    * `%s`: `%s` `%s`\n" % \
                                    (listitem, filename, self._hashdict[filename])
                    else:
                        retstring += "    * `%s`: `%s`\n" % \
                                    (listitem, filename)

        if 'path' in dbentry:
            retstring += "* __Path__: `%s`\n" % dbentry['path']

        if 'file' in dbentry:
            filename = dbentry['file']
            if filename in self._hashdict:
                retstring += "* __File__: `%s` `%s`\n" % \
                             (filename, self._hashdict[filename])
            else:
                retstring += "* __File__: `%s`\n" % dbentry['file']

        if not silent:
            print retstring

        return retstring

    def print_path_db(self, suppress_lists=90, fileobj=None):
        r"""print all the files in the path database; note that it is a
        dictionary and so this is un-ordered. print_path_db_by_group prints the
        database items ordered by the file group specifications.

        You should only need to use this if the db groups are broken.

        Suppress printing of lists of more than `suppress_lists` files.
        If given a filename, this will write to a text file in markdown format.
        """
        print "-" * 80

        for db_key in self._pathdict:
            dbstring = self.print_db_item(db_key,
                                          suppress_lists=suppress_lists)
            print dbstring

        print "-" * 80

    def print_path_db_by_group(self, suppress_lists=90, fileobj=None):
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

            for db_key in self._groups[groupname]:
                dbstring = self.print_db_item(db_key,
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

        print_dictionary(self.gitlog, sys.stdout,
                         key_list=["SHA", "blame", "date"],
                         prepend="git_")


    def fetch_parent(self, db_key, return_path=False):
        r"""Fetch the parent key associated with a file database entry
        if `return_path` return the directory instead of the key
        """
        dbentry = self._pathdict[db_key]

        if 'parent' in dbentry:
            parent_key = dbentry['parent']
            if return_path:
                return self.fetch(parent_key, silent=True)
            else:
                return parent_key
        else:
            print "%s has no parent directory field" % db_key

    def fetch_multi(self, data_obj, db_token="db:", silent=False,
                    intend_read=True):
        r"""Handle various sorts of file pointers/data
        if `data_obj`
            is an array, return a deep copy of it
            is a string:
                if it begins with "db:" -- string after db is a db key
                otherwise assume it is a file and try to open it
        """
        if isinstance(data_obj, str):
            if data_obj[0:len(db_token)] == db_token:
                db_key = data_obj[len(db_token):]
                filename = self.fetch(db_key, intend_read=intend_read,
                                      silent=silent)
            else:
                filename = data_obj
                prefix = "non-db filename "
                ft.path_properties(filename, intend_read=intend_read,
                                   is_file=True,
                                   prefix=prefix, silent=silent)

            ret_data = algebra.make_vect(algebra.load(filename))
        else:
            ret_data = copy.deepcopy(data_obj)

        return ret_data

    def fileset_cases(self, db_key, case_key, divider=";"):
        r"""Parse the list of files in a db entry into unique identifiers
        see behavior of unpack_cases(case_list, case_key, divider=";")
        """
        fdb = self.fetch(db_key, silent=True)
        return unpack_cases(fdb[0], case_key, divider=divider)

    def fetch(self, db_key, pick=None, intend_read=False, intend_write=False,
              purpose="", silent=False, tack_on=None):
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
        # if not given pick explicitly, try to extract it from the given key
        # "key_to_db:pick_entry"
        if not pick:
            db_key = db_key.split(":")
            if len(db_key) == 2:
                main_db = db_key[0]
                pick = db_key[1]
            else:
                main_db = db_key[0]
                pick = None
            db_key = main_db

        dbentry = self._pathdict[db_key]
        prefix = "%s (%s) " % (purpose, db_key)

        if 'file' in dbentry:
            pathout = copy.deepcopy(dbentry['file'])
            if tack_on:
                pathout = tack_on_subdir(pathout, tack_on)

            ft.path_properties(pathout, intend_write=intend_write,
                               intend_read=intend_read, is_file=True,
                               prefix=prefix, silent=silent)

        if 'path' in dbentry:
            pathout = copy.deepcopy(dbentry['path'])
            if tack_on:
                pathout = tack_on_subdir(pathout, tack_on)

            ft.path_properties(pathout, intend_write=intend_write,
                               intend_read=intend_read, is_file=False,
                               prefix=prefix, silent=silent)

        if 'filelist' in dbentry:
            pathout = (copy.deepcopy(dbentry['listindex']), \
                       copy.deepcopy(dbentry['filelist']))

            if tack_on:
                for item in pathout[0]:
                    pathout[1][item] = tack_on_subdir(pathout[1][item],
                                                      tack_on)

            if pick:
                pathout = pathout[1][pick]
                ft.path_properties(pathout, intend_write=intend_write,
                                   intend_read=intend_read, is_file=True,
                                   prefix=prefix, silent=silent)
            else:
                for item in pathout[0]:
                    filepath = pathout[1][item]
                    ft.path_properties(filepath, intend_write=intend_write,
                                       intend_read=intend_read, is_file=True,
                                       file_index=item, prefix=prefix,
                                       silent=silent)

        return pathout

    def generate_path_webpage(self, suppress_lists=90):
        r"""Write out a markdown file with the file database
        """
        localpath = "/cita/d/www/home/eswitzer/GBT_param/"
        localdb = localpath + "path_database.py"
        dbwebpage = localpath + "path_database.txt"

        print "writing markdown website to: " + dbwebpage
        fileobj = open(dbwebpage, "w")

        fileobj.write("Data path DB\n============\n\n")
        fileobj.write("* specified by local path: `%s` and URL `%s`\n" %
                      (localdb, self.db_url))
        fileobj.write("* file hash table specified by `%s`\n" % self.hash_url)
        fileobj.write("* website, checksums compiled: %s by %s\n" % \
                       self.runinfo)
        fileobj.write("* %d files registered; database size in memory = %s\n" %
                     self._db_size)
        self.print_path_db_by_group(suppress_lists=suppress_lists,
                                    fileobj=fileobj)

        fileobj.close()

    def generate_path_webpage_by_group(self, suppress_lists=45):
        r"""Write out a markdown file with the file database
        """
        localpath = "/cita/d/www/home/eswitzer/GBT_param/"
        localdb = localpath + "path_database.py"
        dbwebpage = localpath + "path_database.txt"

        print "writing markdown website to: " + dbwebpage
        fileobj = open(dbwebpage, "w")

        fileobj.write("Data path DB\n============\n\n")
        fileobj.write("* specified by local path: `%s` and URL `%s`\n" %
                      (localdb, self.db_url))
        fileobj.write("* file hash table specified by `%s`\n" % self.hash_url)
        fileobj.write("* website, checksums compiled: %s by %s\n" % \
                       self.runinfo)
        fileobj.write("* %d files registered; database size in memory = %s\n" %
                     self._db_size)

        fileobj.write("****\n")
        for groupname in self._group_order:
            groupfile = localpath + "path_database_group_" + groupname + ".txt"
            grouplink = "[" + groupname + "](path_database_group_" + \
                        groupname + ".html)"
            fileobj.write("* %s\n" % grouplink)

            print groupfile, grouplink

            groupobj = open(groupfile, "w")

            print "%s\n%s\n" % (groupname, "-" * len(groupname))
            groupobj.write("## %s ##\n\n" % groupname)

            for db_key in self._groups[groupname]:
                dbdesc = self.print_db_item_desc(db_key)
                groupobj.write(dbdesc)

            groupobj.write("\n---------\n\n")

            for db_key in self._groups[groupname]:
                dbstring = self.print_db_item(db_key,
                                              suppress_lists=suppress_lists)
                groupobj.write(dbstring + "\n")

            groupobj.close()

        print "-" * 80
        fileobj.close()

    def generate_hashtable(self):
        r"""Write out the file hash table
        """
        localpath = "/cita/d/www/home/eswitzer/GBT_param/"
        dbhashpage = localpath + "hashlist.txt"

        print "writing file hash list to: " + dbhashpage
        hashobj = open(dbhashpage, "w")
        hashdict = {}
        hashlist = []

        hashobj.write("# checksums compiled: %s by %s\n" % \
                       self.runinfo)

        for groupname in self._group_order:
            for db_key in self._groups[groupname]:
                dbentry = self._pathdict[db_key]

                if 'filelist' in dbentry:
                    listindex = dbentry['listindex']
                    for listitem in listindex:
                        filename = dbentry['filelist'][listitem]
                        hashdict[filename] = ft.hashfile(filename)
                        hashlist.append(filename)
                        print "%s: %s" % (filename, hashdict[filename])

                if 'file' in dbentry:
                    filename = dbentry['file']
                    hashdict[filename] = ft.hashfile(filename)
                    hashlist.append(filename)
                    print "%s: %s" % (filename, hashdict[filename])

        print_dictionary(hashdict, hashobj, key_list=hashlist)
        hashobj.close()

def get_mapdict(dir, selection=None):
    r"""
    Generate a map dict according to the map file in a dir
    """
    maplist = os.listdir(dir)
    mapdict = {}
    for map in maplist:
        if os.path.isfile(dir+map) and map.split('.')[-1]=='npy':
            mapsplit = map.split('.')[0].split('_')
            if mapsplit[0] == 'sec':
                #print map
                key1 = mapsplit[1] + '_with_' + mapsplit[7]
                if mapsplit[2] == 'modes':
                    key2 = mapsplit[2]
                else:
                    key2 = mapsplit[4]
                if key2 == 'inv':
                    key2 = mapsplit[3] + '_' + key2
                key3 = mapsplit[-1]

                mapdict['%s;%s;%s'%(key1, key2, key3)] = dir + map
            if mapsplit[0] == 'combined':
                key1 = mapsplit[2]
                key2 = mapsplit[3]

                mapdict['%s;%s'%(key1, key2)] = dir + map

            if mapsplit[0] == 'secA' or\
               mapsplit[0] == 'secB' or\
               mapsplit[0] == 'secC' or\
               mapsplit[0] == 'secD':
                key1 = mapsplit[0][-1]
                key2 = mapsplit[3] + '_' + mapsplit[4]
                #if key2=='noise_inv' and mapsplit[5] == 'diag':
                #    key2='noise_weight'

                mapdict['%s;%s'%(key1, key2)] = dir + map

            if  mapsplit[0] == 'sim' and mapsplit[1] == selection:
                key1 = int(mapsplit[2])

                mapdict['%d'%key1] = dir + map

        if os.path.isfile(dir+map) and map.split('.')[-1]=='pkl':
            mapsplit = map.split('.')[0].split('_')
            if mapsplit[0] == 'SVD':
                #print map
                key1 = mapsplit[2] + '_with_' + mapsplit[4]
                key2 = mapsplit[0]
                #print key1, key2

                mapdict['%s;%s'%(key1, key2)] = dir + map
    
    maps = [mapdict.keys(), mapdict]
    return maps



if __name__ == "__main__":
    import doctest

    # generate the path db markdown website
    datapath_db = DataPath()
    #datapath_db.generate_hashtable()
    #datapath_db.generate_path_webpage()
    datapath_db.generate_path_webpage_by_group()

    # run some tests
    OPTIONFLAGS = (doctest.ELLIPSIS |
                   doctest.NORMALIZE_WHITESPACE)
    doctest.testmod(optionflags=OPTIONFLAGS)
