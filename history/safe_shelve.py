import shelve
import os
import fcntl
import new
import __builtin__
import anydbm
import time


def _close(self):
    shelve.Shelf.close(self)
    fcntl.flock(self.lockfile.fileno(), fcntl.LOCK_UN)
    self.lockfile.close()


# TODO: consider new class like
# http://stackoverflow.com/questions/486490/python-shelve-module-question
# rather than overriding shelve
def open_dumb(filename, flag='c', protocol=-1, writeback=False, block=True,
                lockfilename=None):
    """Open the sheve file, creating a lockfile at filename.lock.  If
    block is False then a IOError will be raised if the lock cannot
    be acquired
    based on activestate: 576591-simple-shelve-with-linux-file-locking
    """
    if lockfilename == None:
        lockfilename = filename + ".lock"

    lockfile = __builtin__.open(lockfilename, 'w')

    # accquire the lock
    if flag == 'r':
        lockflags = fcntl.LOCK_SH
    else:
        lockflags = fcntl.LOCK_EX
    if not block:
        lockflags = fcntl.LOCK_NB
    fcntl.flock(lockfile.fileno(), lockflags)
    print "ACQUIRED LOCK"

    # open the shelf and override its standard methods
    shelf = shelve.open(filename, flag=flag, protocol=protocol, writeback=writeback)
    shelf.close = new.instancemethod(_close, shelf, shelve.Shelf)
    shelf.lockfile = lockfile

    print "RETURNING SHELVE"
    return shelf


