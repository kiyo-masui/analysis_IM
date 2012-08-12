"""manage batch runs and memoization
The MemoizeBatch class is best for trivial parallelization
The memoize_persistent decorator is best for memoized serial calls
"""
from core import algebra
import multiprocessing
import cPickle as pickle
import hashlib
import re
import random
import time
import shelve
import functools
import os
import copy
import logging as log


def log_timing_func():
    '''Decorator generator that logs the time it takes a function to execute'''
    def decorator(func_to_decorate):
        def wrapper(*args, **kwargs):
            start = time.time()
            result = func_to_decorate(*args, **kwargs)
            elapsed = (time.time() - start)

            log.debug("[TIMING] %s: %s" % (func_to_decorate.__name__, elapsed))

            return result
        wrapper.__doc__ = func_to_decorate.__doc__
        wrapper.__name__ = func_to_decorate.__name__
        return wrapper
    return decorator


# for classes
def log_timing(func_to_decorate):
    '''Decorator generator that logs the time it takes a function to execute'''
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func_to_decorate(*args, **kwargs)
        elapsed = (time.time() - start)

        #log.debug("[TIMING] %s: %s" % (func_to_decorate.__name__, elapsed))
        print "[TIMING] %s: %s" % (func_to_decorate.__name__, elapsed)

        return result

    wrapper.__doc__ = func_to_decorate.__doc__
    wrapper.__name__ = func_to_decorate.__name__
    return wrapper


def short_repr(input_var, maxlen=None):
    r"""return a repr() for something if it is short enough
    could also use len(pickle.dumps()) to see if the object and not just its
    repr is too big to bother printing.
    """
    repr_out = repr(input_var)
    if maxlen:
        if len(rep_rout) > maxlen:
            repr_out = "BIG_ARG"

    return repr_out


def rehasher(item):
    r"""classes derived on numpy have a non-uniform way of being hashed, so
    just call their tolist member"""
    try:
        infofield = tuple(sorted(item.info.items()))
        fullarray = (item.tolist(), infofield)
        rehashed = hashlib.sha224(pickle.dumps(fullarray)).hexdigest()
    except:
        rehashed = item

    return rehashed


# TODO: do not print long arguments
def print_call(args_package):
    r"""Print the function and its arguments and the file in which the outputs
    are saved in the cache directory."""
    (signature, directory, funcname, args, kwargs) = args_package

    kwlist = ["%s=%s" % (item, short_repr(kwargs[item])) for item in kwargs]
    kwstring = ", ".join(kwlist)
    argstring = ", ".join([short_repr(item) for item in args])

    filename = "%s/%s.shelve" % (directory, signature)
    filename = re.sub('/+', '/', filename)

    # TODO: sometimes there are spurious commas
    print "%s(%s, %s) -> %s" % (funcname, argstring, kwstring, filename)

    return filename


def repackage_kiyo(kiyo_map_split):
    combined_kiyo = algebra.make_vect(kiyo_map_split[0], axis_names=('freq', 'ra', 'dec'))
    combined_kiyo.info = kiyo_map_split[1]

    return combined_kiyo


memoize_directory = "/mnt/raid-project/gmrt/eswitzer/persistent_memoize/"
def memoize_persistent(func):
    r"""threadsafe shelve using flock was too annoying; here, just wait
    note that shelve protocol=-1 does not handle Kiyo-style array metadata.
    The repackage_kiyo function above handles this.

    This should really be transferred to a database than can handle concurrency
    """
    def memoize(*args, **kwargs):
        funcname = func.__name__
        rehashed = [rehasher(item) for item in args]
        argpkl = pickle.dumps((funcname, rehashed,
                               tuple(sorted(kwargs.items()))), -1)
        arghash = hashlib.sha224(argpkl).hexdigest()

        filename = print_call((arghash, memoize_directory,
                               funcname, rehashed, kwargs))
        done_filename = filename + ".done"
        busy_filename = filename + ".busy"

        # prevent the race condition where many threads want to
        # start writing a new cache file at the same time
        time.sleep(random.uniform(0, 0.5))
        if (os.access(done_filename, os.F_OK) or \
            os.access(busy_filename, os.F_OK)):
            printed = False
            while(not os.access(done_filename, os.F_OK)):
                time.sleep(1.)
                if not printed:
                    print "waiting for %s" % filename
                    printed = True

            # if many threads were waiting to read, space their reads
            # so that they don't attack at once.
            time.sleep(random.uniform(0, 0.5))
            print "ready to read %s" % filename
            try:
                input_shelve = shelve.open(filename, "r", protocol=-1)
                retval = input_shelve['result']
                input_shelve.close()
                print "used cached value %s" % filename
            except:
                raise ValueError
        else:
            # first flag the cachefile as busy so other threads wait
            busyfile = open(busy_filename, "w")
            busyfile.write("working")
            busyfile.close()

            # recalculate the function
            print "no cache, recalculating %s" % filename
            start = time.time()
            retval = func(*args, **kwargs)

            outfile = shelve.open(filename, "n", protocol=-1)
            outfile["signature"] = arghash
            outfile["filename"] = filename
            outfile["funcname"] = funcname
            outfile["args"] = rehashed
            outfile["kwargs"] = kwargs
            outfile["result"] = retval
            outfile.close()
            time.sleep(0.2)     # TODO: this is a trial; remove?

            # indicate that the function is done being recalculated
            donefile = open(done_filename, "w")
            donefile.write("%10.15f\n" % (time.time() - start))
            donefile.close()
            # and remove the busy flag
            os.remove(busy_filename)

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

    # some voodoo to prevent lockfile collisions
    time.sleep(random.uniform(0, 2.))
    result = funcattr(*args, **kwargs)

    outfile = shelve.open(filename, 'n', protocol=-1)
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

    This class circumvents these issues by assigning an SHA hash to the
    arguments and writing a file by that name. If those arguments have been
    called, it can retrieve the output. Note that these have different SHA
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
    '4a1e33dbab4b660fcddbe98d604d87ed28a3dce9c4844190df3c05c5'
    >>> import numpy as np
    >>> bins = np.arange(0,1, 0.2)
    >>> caller1.execute(5, "ok2", arg1="stringy", arg2=bins)
    'cc9dd30e2b656e96bfdec98a219c030d5dc4792391bdfa54ad209d3c'
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
    # TODO: add regenerate option for generate mode:
    # (check to see if the result exists)
    def __init__(self, funcname, directory, generate=False, verbose=False):
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
        self.verbose = verbose

    def execute(self, *args, **kwargs):
        r"""if one of the arguments is "inifile" expand that file reference
        into a string which is added to the argument checksum
        also see stackoverflow's: computing-an-md5-hash-of-a-data-structure
        based on: how-to-memoize-kwargs
        """
        kwini = copy.deepcopy(kwargs)
        if "inifile" in kwargs:
            #if self.verbose:
            #    print "MemoizeBatch: expanding inifile into argument checksum"

            open_inifile = open(kwargs["inifile"], "r")
            kwini["inifile"] = "".join(open_inifile.readlines())
            open_inifile.close()

        # TODO: add support for more generic arguments (numpy derivatives)
        argpkl = pickle.dumps((self.funcname, args,
                               tuple(sorted(kwini.items()))), -1)

        arghash = hashlib.sha224(argpkl).hexdigest()

        args_package = (arghash, self.directory,
                        self.funcname, args, kwargs)

        if self.verbose:
            print_call(args_package)

        if self.generate:
            self.call_stack.append(args_package)
            retval = arghash
        else:
            # TODO: raise NotCalculated?
            # TODO: check args, kwargs with requested?
            filename = "%s/%s.shelve" % (self.directory, arghash)
            input_shelve = shelve.open(filename, "r", protocol=-1)
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
            result = []
            for item in self.call_stack:
                print_call(item)
                result.append(function_wrapper(item))
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
def useless_loop(input_var, arg1='a', arg2=2):
    r"""Useless function to test the memoize decorator"""
    return (input_var, arg1, arg2)


class memoized(object):
    """Decorator that caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned, and
    not re-evaluated.

    e.g.:
    @memoized
    def fibonacci(n):
        "Return the nth fibonacci number."
        if n in (0, 1):
            return n
        return fibonacci(n-1) + fibonacci(n-2)

     print fibonacci(12)
    """
    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        try:
            return self.cache[args]
        except KeyError:
            value = self.func(*args)
            self.cache[args] = value
            return value
        except TypeError:
            # uncachable -- for instance, passing a list as an argument.
            # Better to not cache than to blow up entirely.
            return self.func(*args)

    def __repr__(self):
        """Return the function's docstring."""
        return self.func.__doc__

    def __get__(self, obj, objtype):
        """Support instance methods."""
        return functools.partial(self.__call__, obj)


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
