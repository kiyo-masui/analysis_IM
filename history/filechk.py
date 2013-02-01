import os
import stat
import sys
import time

def path_properties(pathname, intend_write=False, intend_read=False,
                    prepend="", silent=False, file=True):
    r"""Check various properties about the path, print the result
    if `intend_write` then die if the path does not exist or is not writable
    if `intend_read` then die if the path does not exist
    with a `prepend` string the filename, whether it is writable, and
    what its last modification date was if the file exists.

    # Usage examples:
    >>> path_properties("this_does_not_exist.txt")
    file this_does_not_exist.txt: does not exist,
        but path: ./ exists, is writable and not readable.
    (True, False, True, False, '...')

    >>> path_properties("/tmp/nor_does_this.txt")
    file /tmp/nor_does_this.txt: does not exist,
    but path: /tmp/ exists, is writable and not readable.
    (True, False, True, False, '...')

    # here's one that you should not be able to write to
    >>> path_properties("/var/nor_does_this.txt")
    file /var/nor_does_this.txt: does not exist,
        but path: /var/ exists, is not writable and not readable.
    (True, False, False, False, '...')

    # here's one that definitely exists
    >>> path_properties(__file__)
    file ...: exists, is writable and readable.
    (True, True, True, False, '...')

    # test a writeable directory that exists
    >>> path_properties('/tmp/')
    file /tmp/: exists, is writable and readable.
    (True, True, True, True, '...')
    """
    entry = "file %s%s:" % (prepend, pathname)

    exists = os.access(pathname, os.F_OK)
    readable = os.access(pathname, os.R_OK)
    writable = os.access(pathname, os.W_OK)
    executable = os.access(pathname, os.X_OK)
    if exists:
        modtime = time.ctime(os.path.getmtime(pathname))
    else:
        modtime = None

    # if this is a file that does not exist, check the directory
    pathexists = False
    if not exists and file:
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
            entry += " is writable"
        else:
            entry += " is not writable"

        if readable:
            entry += " and readable."
        else:
            entry += " and not readable."
    else:
        entry += " does not exist"

    if not silent:
        print entry

    if intend_read and not readable:
        sys.exit()

    if intend_write and not writable:
        sys.exit()

    return (exists, readable, writable, executable, modtime)


if __name__ == "__main__":
    import doctest
    OPTIONFLAGS = (doctest.ELLIPSIS |
                   doctest.NORMALIZE_WHITESPACE)
    doctest.testmod(optionflags=OPTIONFLAGS)

