import shelve
import fcntl
import anydbm
import time
import os

# this doesn't work so well, but doing this sort of locking is fraught with
# problems
class SafeShelve(shelve.Shelf):
    def __init__(self, filename, flag='c', protocol=-1, writeback=False,
                 block=True):

        self.lockfilename = filename + ".lock"

        # accquire the lock
        if flag == 'r':
            lockflags = fcntl.LOCK_SH
            lockfilemode = "r"
        else:
            lockflags = fcntl.LOCK_EX
            lockfilemode = "a"
        if not block:
            lockflags = fcntl.LOCK_NB

        self.lockfile = open(self.lockfilename, lockfilemode)
        fcntl.flock(self.lockfile.fileno(), lockflags)
        print "ACQUIRED LOCK"

        #for retry in (2,1,0):
        #    try:
        #        db = anydbm.open(filename, flag=flag)
        #    except anydbm.error, val:
        #        raise ValueError(val)
        #    except:
        #        if not retry:
        #            raise
        #        time.sleep(2 - (0.9*retry))

        time.sleep(0.1)
        db = anydbm.open(filename, flag=flag)

        # this does not work with super() ?
        shelve.Shelf.__init__(self, db, protocol=protocol, writeback=writeback)

    def close(self):
        print "closing this"

        shelve.Shelf.close(self)
        fcntl.flock(self.lockfile.fileno(), fcntl.LOCK_UN)
        self.lockfile.close()
        os.remove(self.lockfilename)
