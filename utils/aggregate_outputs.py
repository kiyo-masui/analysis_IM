import multiprocessing
import shelve

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


def print_call(args_package):
    (execute_key, funcname, args, kwargs) = args_package

    kwlist = ["%s=%s" % (item, short_repr(kwargs[item])) for item in kwargs]
    kwstring = ", ".join(kwlist)
    argstring = ", ".join([short_repr(item) for item in args])

    print "%s: %s(%s, %s)" % (execute_key, funcname, argstring, kwstring)


def function_wrapper(args_package):
    (execute_key, funcname, args, kwargs) = args_package
    print_call(args_package)

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

    return (args_package, funcattr(*args, **kwargs))


class AggregateOutputs(object):
    r"""
    test_agg = AggregateOutputs("test_function")
    test_agg.execute("ok", "ok2", kwarg="that", execute_key="one")
    test_agg.execute("ok3", "ok4", kwarg="that", execute_key="two")
    test_agg.execute("ok5", "ok6", kwarg="that", execute_key="three")
    test_agg.multiprocess_stack("test.shelve", debug=True)

    map returns a list of results and we would like to associate a tag to each
    result and build up a dictionary of results. Because we're handing args and
    kwargs to the function wrapper, there is not a convenient
    """
    def __init__(self, funcname, verbose=False):
        self.call_stack = []
        self.funcname = funcname
        self.verbose = verbose

    def execute(self, *args, **kwargs):
        execute_key = kwargs["execute_key"]
        # delete the kwarg for this function to associate an ID to the output
        # so it does not interfere with the function call.
        del kwargs["execute_key"]

        args_package = (execute_key, self.funcname, args, kwargs)
        if self.verbose:
            print_call(args_package)

        self.call_stack.append(args_package)

    def multiprocess_stack(self, filename, save_cpu=4, debug=False, ncpu=8):
        r"""process the call stack built up by 'execute' calls using
        multiprocessing.
        filename is the output to write the shelve to
        `save_cpu` is the number of CPUs to leave free
        `debug` runs one process at a time because of funny logging/exception
        handling in multiprocessing
        """
        if debug:
            results = []
            for item in self.call_stack:
                print item
                results.append(function_wrapper(item))
        else:
            num_cpus = multiprocessing.cpu_count() - save_cpu
            if ncpu:
                num_cpus = ncpu
            pool = multiprocessing.Pool(processes=num_cpus)
            results = pool.map(function_wrapper, self.call_stack)
            pool.close()

        # after running the jobs reset the batch
        self.call_stack = []

        print "multiprocessing_stack: jobs finished"
        outshelve = shelve.open(filename, "n", protocol=-1)
        for result_item in results:
            args_package = result_item[0]
            execute_key = args_package[0]
            outshelve[execute_key] = result_item

        outshelve.close()

def test_function(arg1, arg2, kwarg=None):
    return arg1 + arg2 + repr(kwarg)
