import multiprocessing
import cPickle as pickle
import hashlib
import numpy as np
import re
import shelve
import sys

# TODO: incorporate class persistence to have this save itself in the directory
# save itself using presistent class (selected variables)
# list directory -> list of hashes known?
# using persist .load to reconstruct class and make hash -> variables list
# function: regenerate regardless of whether hash exists
# initialization: generation mode or consumer mode
# generation mode: always return md5; only add to stack if not known
# consumer mode: raise exception if not known, return result if known
# TODO: write PBS caller


def function_wrapper(args):
    r"""this wrapper supports MemoizeBatch and circumvents the fact that
    multiprocessing pool's map can not handle class functions or generic
    function arguments to calls in that pool. Note that we also save the data
    here as opposed to handing it back to multiprocessing. This is to prevent
    the outputs from a potentially large number of function calls from being
    saved in memory.
    """
    (hash, directory, funcname, args, kwargs) = args

    funcattr = None
    funcsplit = funcname.split(".")

    kwlist = []
    for item in kwargs:
        kwlist.append("%s=%s" % (item, repr(kwargs[item])))

    kwstring = ", ".join(kwlist)

    arglist = []
    for item in args:
        arglist.append(repr(item))

    argstring = ", ".join(arglist)

    filename = "%s/%s.shelve" % (directory, hash)
    filename = re.sub('/+', '/', filename)

    # TODO: sometimes there are spurious commas
    print "%s(%s, %s) -> %s" % (funcname, argstring, kwstring, filename)

    if len(funcsplit) > 1:
        mod = __import__(".".join(funcsplit[0:-1]))
        for comp in funcsplit[1:-1]:
            mod = getattr(mod, comp)

        funcattr = getattr(mod, funcsplit[-1])
    else:
        funcattr = globals()[funcsplit[0]]

    result = funcattr(*args, **kwargs)
    print result
    # TODO now save the outputs

    outfile = shelve.open(filename, 'n')
    outfile["hash"] = hash
    outfile["filename"] = filename
    outfile["funcname"] = funcname
    outfile["args"] = args
    outfile["kwargs"] = kwargs
    outfile["result"] = result
    outfile.close()

    return hash


class MemoizeBatch(object):
    r"""
    initialize a class that points to the function you want to call
    call execute of that class
       -if calculated, fetch value based on hash
       -if not: put the parameters on a heap and assign a hash

    call the parallel exec: either multiproc or PBS

    use the same expressions to recall the data
    """
    def __init__(self, funcname, directory, generate=False, regenerate=False):
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

    def execute(self, *args, **kwargs):
        r"""
        """
        argpkl = pickle.dumps((self.funcname, args,
                               tuple(sorted(kwargs.items()))), -1)
        arghash = hashlib.md5(argpkl).hexdigest()

        if self.generate or self.regenerate:
        # TODO: if not regenerate, check to see if the result exists

            self.call_stack.append((arghash, self.directory,
                                    self.funcname, args, kwargs))

            retval = arghash
        else:
            #TODO: raise NotCalculated
            filename = "%s/%s.shelve" % (self.directory, arghash)
            input_shelve = shelve.open(filename, "r")
            # check the values for sanity's sake (can prob. remove)
            if ((input_shelve['hash'] != arghash) or
               (input_shelve['funcname'] != self.funcname) or
               (input_shelve['args'] != args) or
               (input_shelve['kwargs'] != kwargs)):
                print "ERROR: cache file is corrupted"
                sys.exit()

            retval = input_shelve['result']
            input_shelve.close()

        return retval


    def multiprocess_stack(self, save_cpu=4):
        r"""save_cpu is the number of CPUs to leave free
        """
        num_cpus = multiprocessing.cpu_count() - save_cpu
        pool = multiprocessing.Pool(processes=num_cpus)
        result = pool.map(function_wrapper, self.call_stack)
        pool.close()
        print result

    def pbsprocess_stack(self):
        print "PBS batch not implemented yet"


def basic_function(thing1, thing2, arg1="e", arg2="r"):
    return thing1, thing2, arg1, arg2


if __name__ == "__main__":
    import doctest

    caller1 = MemoizeBatch("basic_function", "./testdata/", generate=True)
    caller1.execute("ok",3, arg1=True, arg2=3)
    bins = np.arange(0,1, 0.2)
    caller1.execute("ok2",4, arg1="weee", arg2=bins)
    caller1.execute("ok2",4)
    caller1.multiprocess_stack()

    print "Now testing a function from another class -----------"
    caller2 = MemoizeBatch("utils.data_paths.extract_split_tag", "./testdata/",
                           generate=True)
    caller2.execute(['A;ok', 'B;ok'], divider=';')
    caller2.execute(['B;oke', 'C;okd'], divider=';')
    caller2.multiprocess_stack()

    print "Now calling a new run -----------"
    caller3 = MemoizeBatch("utils.data_paths.extract_split_tag", "./testdata/")
    ret1 = caller3.execute(['A;ok', 'B;ok'], divider=';')
    ret2 = caller3.execute(['B;oke', 'C;okd'], divider=';')
    print "wow", ret1
    print "wow", ret2

    #print bins
    #wrapper_function("ok", "ok2", arg1="that", arg2=bins)
