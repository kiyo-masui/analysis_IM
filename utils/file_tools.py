import os
import cPickle
import hashlib
import sys
import time
import shelve
# TODO: use path_properties to check on files other functions here try to write


def save_pickle(pickle_data, filename):
    r"""wrap cPickle; useful for general class load/save
    note that if you load a pickle outside of the class that saved itself, you
    must fully specify the class space as if you were the class, e.g.:
    from correlate.freq_slices import * to import the freq_slices object
    This clobbers anything in the requested file.
    """
    pickle_out_root = os.path.dirname(filename)
    if not os.path.isdir(pickle_out_root):
        os.mkdir(pickle_out_root)

    pickle_handle = open(filename, 'wb')
    cPickle.dump(pickle_data, pickle_handle)
    pickle_handle.close()


def load_pickle(filename):
    r"""Return `pickle_data` saved in file with name `filename`.
    """
    pickle_handle = open(filename, 'r')
    pickle_data = cPickle.load(pickle_handle)
    pickle_handle.close()
    return pickle_data


def extract_rootdir(pathname):
    r"""superceded by os.path.dirname(filename), but we still use it in
    path_properties because "./" is nice to indicate local paths.
    """
    directory = "/".join(pathname.split("/")[:-1])
    if directory == "":
        directory = "."
    directory += "/"
    return directory


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
        directory = extract_rootdir(pathname)

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
        print "ERROR: no file to read"
        sys.exit()

    if intend_write and not writable:
        print "ERROR: can not write this file"
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
            hash_digest = hasher.hexdigest()
        else:
            hash_digest = "too_big"

        modtime = time.ctime(os.path.getmtime(filename))
    else:
        hash_digest = "not_exist"
        modtime = "not_exist"

    return (hash_digest, modtime)


class ClassPersistence(object):
    r"""note that pickle files are convenient but inflexible in that they
    depend on the context of the class instance, so are a pain to open later
    in different contexts -- they are fine to purely save intermediate data.
    """
    def __init__(self):
        print "This is your ClassPersistence base class talking"

    def save_pickle(self, filename):
        """note that to withold certain attributes (like open file handles)
        from the pickle file to save, use __getstate__ and __setstate__
        as in http://docs.python.org/library/pickle.html#example to delete
        certain items from the class __dict__."""
        print "save_pickle: to file " + filename
        save_pickle(self, filename)

    @classmethod
    def load_pickle(cls, filename):
        r"""reinvigorate a class from a pickle which has saved everything out:
        ok = ClassPersistence.load_pickle(filename)
        """
        print "load_pickle: from file " + filename
        return load_pickle(filename)

    def shelve_variables(self, filename, varlist=None):
        r"""save the variables in a class, optionally just those named in
        `varlist`. Note that __dict__ of an instance is just user-provided
        variables, but Class.__dict__ here is everything in the class.
        This clobbers anything in the requested file open as 'n'
        """
        shelve_out_root = os.path.dirname(filename)
        if not os.path.isdir(shelve_out_root):
            os.mkdir(shelve_out_root)

        shelveobj = shelve.open(filename, 'n')

        if varlist is None:
            print "shelve_variables: shelving all var. %s to file %s" % \
                  (repr(self.__dict__.keys()), filename)
            shelveobj.update(self.__dict__)
        else:
            print "shelve_variables: shelving %s to file %s" % \
                  (varlist, filename)
            for key in varlist:
                shelveobj[key] = self.__dict__[key]

        shelveobj.close()

    def load_variables(self, filename, varlist=None):
        r"""load variables from a shelve directly into class attributes"""
        shelveobj = shelve.open(filename, 'r')
        if varlist is None:
            varlist = shelveobj.keys()

        for key in varlist:
            setattr(self, key, shelveobj[key])
        shelveobj.close()


class TestClassPersistence(ClassPersistence):
    r"""Example class for ClassPersistence object
    One option the load_variables provides is to split the __init__ into to two
    cases; 1) usual initialization with all variables handed in, 2) a shelve
    file initialization where some key variables are refreshed from the shelve
    file.

    # test the full recovery via pickle files
    >>> test = TestClassPersistence('test', [[0,0],[1,1]])
    >>> print "original class __dict__:" + repr(test.__dict__)
    original class __dict__:{'var1': 'test', 'var2': [[0, 0], [1, 1]]}
    >>> pklfile = "/tmp/testClassPersistence.pkl"
    >>> shelvefile = "/tmp/testClassPersistence.shelve"

    # test the pickle
    >>> test.save_pickle(pklfile)
    save_pickle: to file /tmp/testClassPersistence.pkl
    >>> test2 = TestClassPersistence.load_pickle(pklfile)
    load_pickle: from file /tmp/testClassPersistence.pkl
    >>> test2.print_var()
    print_var(): 'test' [[0, 0], [1, 1]]

    # test the shelve
    >>> test.shelve_variables(shelvefile)
    shelve_variables: shelving all var. ['var1', 'var2'] to file
        /tmp/testClassPersistence.shelve
    >>> testr = shelve.open(shelvefile)
    >>> print "recovered shelve :" + repr(testr)
    recovered shelve :{'var1': 'test', 'var2': [[0, 0], [1, 1]]}
    >>> testr.close()

    >>> test2 = TestClassPersistence(shelve_filename=shelvefile)
    >>> print "shelve-loaded class __dict__:" + repr(test2.__dict__)
    shelve-loaded class __dict__:{'var1': 'test', 'var2': [[0, 0], [1, 1]]}

    # with only one variable
    >>> test.shelve_variables(shelvefile, varlist=['var1'])
    shelve_variables: shelving ['var1'] to file
        /tmp/testClassPersistence.shelve
    >>> testr = shelve.open(shelvefile)
    >>> print "recovered shelve: " + repr(testr)
    recovered shelve: {'var1': 'test'}
    >>> testr.close()
    >>> test3 = TestClassPersistence(shelve_filename=shelvefile)
    >>> print "reduced shelve-loaded class __dict__:" + repr(test3.__dict__)
    reduced shelve-loaded class __dict__:{'var1': 'test'}
    >>> os.remove(pklfile)
    >>> os.remove(shelvefile)
    """
    def __init__(self, *args, **kwargs):
        if ((len(args) == 0) and ("shelve_filename" in kwargs)):
            self.shelve_init(*args, **kwargs)
        else:
            self.standard_init(*args, **kwargs)

    def shelve_init(self, *args, **kwargs):
        shelve_filename = kwargs['shelve_filename']
        self.load_variables(shelve_filename)

    def standard_init(self, *args, **kwargs):
        self.var1 = args[0]
        self.var2 = args[1]

    def print_var(self):
        print "print_var(): " + repr(self.var1) + " " + repr(self.var2)


if __name__ == "__main__":
    import doctest

    # run some tests
    OPTIONFLAGS = (doctest.ELLIPSIS |
                   doctest.NORMALIZE_WHITESPACE)
    doctest.testmod(optionflags=OPTIONFLAGS)
