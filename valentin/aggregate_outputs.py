import multiprocessing
import shelve
import valentin.agg_utils as utils
import time
# TODO make this work with hdf5 files across multiple nodes
# each function call required to write one hd5 "outfile"


def _function_wrapper(args_package):
    """helper to wrap function evaluation
    TODO: use multiprocessing logger in case of error?
    multiprocessing.get_logger().error("f%r failed" % (arg_kwargs,))
    """
    (execute_key, funcname, args, kwargs) = args_package
    return (args_package, utils.func_exec(funcname, args, kwargs))


def example_function(arg1, arg2, kwarg=None):
    time.sleep(1)
    return arg1 + arg2, kwarg


class AggregateOutputs(object):
    r"""Spin off a stack of function calls in parallel

    >>> test_agg = AggregateOutputs("aggregate_outputs.example_function")
    >>> test_agg.execute("a1", "a2", execute_key="one")
    >>> test_agg.execute("a3", "a4", kwarg="ok", execute_key="two")
    >>> test_agg.execute("a5", "a6", kwarg="no", execute_key="three")
    >>> print test_agg.multiprocess_stack()
    multiprocessing_stack: jobs finished
    {'three': ('a5a6', 'no'), 'two': ('a3a4', 'ok'), 'one': ('a1a2', None)}
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
            print (args_package)

        self.call_stack.append(args_package)

    def multiprocess_stack(self, filename=None, save_cpu=None,
                           debug=False, ncpu=8):
        r"""process the call stack built up by 'execute' calls using
        multiprocessing.

        `filename` is the output to write the shelve to; in this case, save all
        the function arguments for recordkeeping.

        `debug` runs one process at a time because of funny logging/exception
        handling in multiprocessing

        Can do either:
        `save_cpu` is the number of CPUs to leave free
        `ncpu` is the number of CPUs to leave free

        good for many CPU-heavy processes on one node
        all results stay in memory
        """
        if debug:
            results = []
            for item in self.call_stack:
                print item
                results.append(_function_wrapper(item))
        else:
            if ncpu:
                num_cpus = ncpu
            else:
                num_cpus = multiprocessing.cpu_count() - save_cpu

            pool = multiprocessing.Pool(processes=num_cpus)
            results = pool.map(_function_wrapper, self.call_stack)
            pool.close()
            pool.join()

        # after running the jobs reset the batch
        self.call_stack = []

        print "multiprocessing_stack: jobs finished"

        if filename:
            outshelve = shelve.open(filename, "n", protocol=-1)
            for result_item in results:
                args_package = result_item[0]
                execute_key = args_package[0]
                outshelve[execute_key] = result_item[1]

            outshelve.close()
        else:
            outdict = {}
            for result_item in results:
                args_package = result_item[0]
                execute_key = args_package[0]
                outdict[execute_key] = result_item[1]

            return outdict


if __name__ == "__main__":
    import doctest

    OPTIONFLAGS = (doctest.ELLIPSIS |
                   doctest.NORMALIZE_WHITESPACE)
    doctest.testmod(optionflags=OPTIONFLAGS)
