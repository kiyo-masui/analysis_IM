"""manage batch runs and memoization
The MemoizeBatch class is best for trivial parallelization
The memoize_persistent decorator is best for memoized serial calls
"""
import multiprocessing
import cPickle as pickle
import hashlib
import re
import shelve
import os
import fcntl
import new
import functools
import __builtin__


# TODO: do not print long arguments
def print_call(args_package):
    r"""Print the function and its arguments and the file in which the outputs
    are saved in the cache directory."""
    (signature, directory, funcname, args, kwargs) = args_package

    kwlist = ["%s=%s" % (item, repr(kwargs[item])) for item in kwargs]
    kwstring = ", ".join(kwlist)
    argstring = ", ".join([repr(item) for item in args])

    filename = "%s/%s.shelve" % (directory, signature)
    filename = re.sub('/+', '/', filename)

    # TODO: sometimes there are spurious commas
    print "%s(%s, %s) -> %s" % (funcname, argstring, kwstring, filename)

    return filename


def _safe_close(self):
    shelve.Shelf.close(self)
    fcntl.flock(self.lockfile.fileno(), fcntl.LOCK_UN)
    self.lockfile.close()


# TODO: consider new class like
# http://stackoverflow.com/questions/486490/python-shelve-module-question
# rather than overriding shelve
def safe_shelve(filename, flag='c', protocol=-1, writeback=False, block=True,
                lockfilename=None):
    """Open the sheve file, createing a lockfile at filename.lock.  If
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

    # open the shelf and override its standard methods
    shelf = shelve.open(filename, flag=flag, protocol=protocol, writeback=writeback)
    shelf.close = new.instancemethod(_safe_close, shelf, shelve.Shelf)
    shelf.lockfile = lockfile

    return shelf


memoize_directory = "/mnt/raid-project/gmrt/eswitzer/persistent_memoize/"
def memoize_persistent(func):
    def memoize(*args, **kwargs):
        funcname = func.__name__
        argpkl = pickle.dumps((funcname, args,
                               tuple(sorted(kwargs.items()))), -1)
        arghash = hashlib.md5(argpkl).hexdigest()

        filename = print_call((arghash, memoize_directory,
                               funcname, args, kwargs))

        try:
            #input_shelve = shelve.open(filename, "r")
            input_shelve = safe_shelve(filename, "r")
            retval = input_shelve['result']
            input_shelve.close()
            print "used cached value %s" % filename
        except:
            print "no cache, recalculating %s" % filename
            retval = func(*args, **kwargs)
            #outfile = shelve.open(filename, "n")
            outfile = safe_shelve(filename, "n")
            outfile["signature"] = arghash
            outfile["filename"] = filename
            outfile["funcname"] = funcname
            outfile["args"] = args
            outfile["kwargs"] = kwargs
            outfile["result"] = retval
            outfile.close()

        return retval

    return functools.update_wrapper(memoize, func)


def function_wrapper(args_package):
    r"""A free-standing function wrapper that supports MemoizeBatch
    (Why? Multiprocessing pool's map can not handle class functions or generic
    function arguments to calls in that pool.)
    Data are saved here rather than handed back to avoid the scenario where
    all of the output from a batch run is held in memory.
    """
    (signature, directory, funcname, args, kwargs) = args_package

    filename = print_call(args_package)

    funcattr = None
    funcsplit = funcname.split(".")

    # consider http://docs.python.org/dev/library/importlib.html
    if len(funcsplit) > 1:
        mod = __import__(".".join(funcsplit[0:-1]))
        for comp in funcsplit[1:-1]:
            mod = getattr(mod, comp)

        funcattr = getattr(mod, funcsplit[-1])
    else:
        funcattr = globals()[funcsplit[0]]

    result = funcattr(*args, **kwargs)

    #outfile = shelve.open(filename, 'n')
    outfile = safe_shelve(filename, 'n')
    outfile["signature"] = signature
    outfile["filename"] = filename
    outfile["funcname"] = funcname
    outfile["args"] = args
    outfile["kwargs"] = kwargs
    outfile["result"] = result
    outfile.close()

    return signature


class MemoizeBatch(object):
    r"""A class to manage a large batch of function calls and outputs.

    When to use this: you have a function that runs slowly and would like to
    run it many times, trivially parallelized over many parameters. A common
    problem is that one code typically does all the loops for the batch
    processing -> dump to files, and the a second code reads the files to make
    a final product. There is overhead in naming the files and also in the fact
    that the batch caller loops are often structurally different from the loops
    over parameters in the final consumer.

    This class circumvents these issues by assigning an md5 hash to the
    arguments and writing a file by that name. If those arguments have been
    called, it can retrieve the output. Note that these have different md5
    signatures for func(arg=True): call func(arg=True) and call func()
    (even though the result will be the same.) Also, dictionaries should not
    be given as arguments because they are not sorted. Simplejson has
    sort_keys=True but can not serialize numpy. Consider a more general
    argument serialization that uses numpy .tolist() and json, or does a
    nested sort and uses pickle.

    To specify the batch of parameters to run over, initialize the class with
    `generate=True`. Then when execute() is called, it stacks up the list of
    parameters in the call in `self.call_stack`. To perform the computation,
    call one of the call_stack handlers (either multiprocess_stack or
    pbsprocess_stack) to parallelize the function call.

    To use the cached values, initialize the class with `generate=False`
    (default) and then call execute in the same way as above. The cached
    results will be returned by looking up the hash for the arguments.

    Example to generate the cache table:
    >>> caller1 = MemoizeBatch("trial_function", "./testdata/", generate=True)
    >>> caller1.execute("ok",3, arg1=True, arg2=3)
    '1d6e885be47e6f350befb10fd4b28318'
    >>> import numpy as np
    >>> bins = np.arange(0,1, 0.2)
    >>> caller1.execute(5, "ok2", arg1="stringy", arg2=bins)
    'fc3eb86900a1f325fb3369d7f686d475'
    >>> caller1.execute("ok2",4)
    '5b81c329717a428261f3bebe2b4b5bba'
    >>> caller1.multiprocess_stack()
    ['1d6e885be47e6f350befb10fd4b28318',
     'fc3eb86900a1f325fb3369d7f686d475',
     '5b81c329717a428261f3bebe2b4b5bba']

    Example to use the cache in later computation:
    >>> caller2 = MemoizeBatch("trial_function", "./testdata/")
    >>> caller2.execute("ok",3, arg1=True, arg2=3)
    ('ok', 3, True, 3)
    >>> caller2.execute(5, "ok2", arg1="stringy", arg2=bins)
    (5, 'ok2', 'stringy', array([ 0. ,  0.2,  0.4,  0.6,  0.8]))
    >>> caller2.execute("ok2",4)
    ('ok2', 4, 'e', 'r')

    """
    def __init__(self, funcname, directory, generate=False, regenerate=False,
                 verbose=False):
        r"""
        Parameters:
        -----------
        funcname: string
            the function name pointer in one of the forms
            [etc].[package].[module].[function_name]
        directory: string
            directory in which to write the output database
        generate: boolean
            choose this to generate the result cache
        """
        self.funcname = funcname
        self.directory = directory
        self.call_stack = []
        self.generate = generate
        self.regenerate = regenerate
        self.verbose = verbose

    def execute(self, *args, **kwargs):
        # also see stackoverflow's: computing-an-md5-hash-of-a-data-structure
        # based on: how-to-memoize-kwargs
        argpkl = pickle.dumps((self.funcname, args,
                               tuple(sorted(kwargs.items()))), -1)
        arghash = hashlib.md5(argpkl).hexdigest()

        if self.generate or self.regenerate:
        # TODO: if not regenerate, check to see if the result exists
            args_package = (arghash, self.directory,
                            self.funcname, args, kwargs)
            self.call_stack.append(args_package)
            if self.verbose:
                print_call(args_package)

            retval = arghash
        else:
            # TODO: raise NotCalculated?
            # TODO: check args, kwargs with requested?
            filename = "%s/%s.shelve" % (self.directory, arghash)
            #input_shelve = shelve.open(filename, "r")
            input_shelve = safe_shelve(filename, "r")
            retval = input_shelve['result']
            input_shelve.close()

        return retval

    def multiprocess_stack(self, save_cpu=4, debug=False):
        r"""process the call stack built up by 'execute' calls using
        multiprocessing.
        `save_cpu` is the number of CPUs to leave free
        `debug` runs one process at a time because of funny logging/exception
        handling in multiprocessing
        """
        if debug:
            results = []
            for item in self.call_stack:
                results.append(function_wrapper(item))
        else:
            num_cpus = multiprocessing.cpu_count() - save_cpu
            pool = multiprocessing.Pool(processes=num_cpus)
            result = pool.map(function_wrapper, self.call_stack)
            pool.close()

        print result

    def pbsprocess_stack(self):
        r"""process the call stack built up by 'execute' call using the PBS
        batch processing protocol
        """
        print "PBS batch not implemented yet"

    def report_cache_table(self):
        r"""read through the cache directory and make a report of the arguments
        of function for which its results are cached.
        """
        print "report cache is not implemented yet"


def trial_function(thing1, thing2, arg1="e", arg2="r"):
    r"""a test function (won't go in docstring"""
    return thing1, thing2, arg1, arg2


@memoize_persistent
def useless_loop(x, arg1='a', arg2=2):
    return (x, arg1, arg2)


if __name__ == "__main__":
    import doctest

    OPTIONFLAGS = (doctest.ELLIPSIS |
                   doctest.NORMALIZE_WHITESPACE)
    doctest.testmod(optionflags=OPTIONFLAGS)

    print useless_loop(10, arg1=3, arg2='b')
    print useless_loop(10, arg2='b', arg1=4)
    print useless_loop("w", arg2="ok")
    print useless_loop(10)
    print useless_loop("w")

