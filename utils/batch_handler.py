"""manage batch runs and memoization
The MemoizeBatch class is best for trivial parallelization
The memoize_persistent decorator is best for memoized serial calls
"""
import multiprocessing
import cPickle as pickle
import hashlib
import re
import random
import time
import shelve
import functools
import anydbm
import os
import copy

def short_repr(input, maxlen=256):
    r"""return a repr() for something if it is short enough
    could also use len(pickle.dumps()) to see if the object and not just its
    repr is too big to bother printing.
    """
    reprout = repr(input)
    #return "BIG_ARG" if len(reprout) > maxlen else reprout
    return repr(input)

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


memoize_directory = "/mnt/raid-project/gmrt/eswitzer/persistent_memoize/"
def memoize_persistent(func):
    r"""threadsafe shelve using flock was too annoying; here, just wait
    note that shelve protocol=-1 does not handle Kiyo-style array metadata
    There seems to be a race condition where two threads can open the same
    shelve because of the finite time it takes to write the shelve. Currently
    trying a simpler lock.
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

        #if os.access(filename, os.F_OK):
        if (os.access(done_filename, os.F_OK) or \
            os.access(busy_filename, os.F_OK)):
            printed = False
            while(not os.access(done_filename, os.F_OK)):
                time.sleep(1.)
                if not printed:
                    print "waiting for %s" % filename
                    printed = True

            time.sleep(random.uniform(0, 2.))
            print "ready to read %s" % filename
            try:
                input_shelve = shelve.open(filename, "r", protocol=-1)
                retval = input_shelve['result']
                input_shelve.close()
                print "used cached value %s" % filename
            except:
                raise ValueError
        else:
            busyfile = open(busy_filename, "w")
            busyfile.write("working")
            busyfile.close()
            print "no cache, recalculating %s" % filename
            outfile = shelve.open(filename, "n", protocol=-1)
            start = time.time()
            retval = func(*args, **kwargs)
            outfile["signature"] = arghash
            outfile["filename"] = filename
            outfile["funcname"] = funcname
            outfile["args"] = rehashed
            outfile["kwargs"] = kwargs
            outfile["result"] = retval
            outfile.close()
            donefile = open(done_filename, "w")
            donefile.write("%10.15f\n" % (time.time() - start))
            donefile.close()
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
        r"""if one of the arguments is "inifile" expand that file reference
        into a string which is added to the argument checksum
        also see stackoverflow's: computing-an-md5-hash-of-a-data-structure
        based on: how-to-memoize-kwargs
        """
        kwini = copy.deepcopy(kwargs)
        if "inifile" in kwargs:
            print "MemoizeBatch: expanding inifile into argument checksum"
            f = open(kwargs["inifile"], "r")
            kwini["inifile"] = "".join(f.readlines())

        # TODO: add support for more generic arguments (numpy derivatives)
        argpkl = pickle.dumps((self.funcname, args,
                               tuple(sorted(kwini.items()))), -1)

        arghash = hashlib.sha224(argpkl).hexdigest()

        args_package = (arghash, self.directory,
                        self.funcname, args, kwargs)

        if self.verbose:
            print_call(args_package)

        if self.generate or self.regenerate:
        # TODO: if not regenerate, check to see if the result exists
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

