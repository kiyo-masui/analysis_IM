import time
import string
import os
import logging as log


def get_env(name):
    r"""retrieve a selected environment variable"""
    try:
        env_var = os.environ[name]
    except KeyError:
        print "your environment variable %s is not set" % name
        sys.exit()

    return env_var


def repl_bracketed_env(string_in):
    r"""replace any []-enclosed all caps subtring [ENV_VAR] with its
    environment vairable

    >>> repl_bracketed_env("this is your path: [PATH]")
    'this is your path: ...'
    """
    ttable = string.maketrans("[]", "[[")
    breakup = string_in.translate(ttable).split("[")
    retval = []
    for item in breakup:
        if item.isupper():
            item = get_env(item)

        retval.append(item)

    retval = ''.join(retval)

    return retval


def log_timing_func():
    '''Decorator that logs the time it takes a function to execute
    example for timing your_function():

    from process_tools import utils as proc_util

    @proc_util.log_timing_func
    def your_function():
    '''
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
    '''Decorator that logs the time it takes a function to execute
    example for timing __init__:

    from process_tools import utils as proc_util

    @proc_util.log_timing
    def __init__(self, ...):
    '''
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


def short_repr(input_var, maxlen=None, repl="BIG_ARG"):
    r"""return a repr() for something if it is short enough

    >>> short_repr(3, maxlen=30)
    '3'
    >>> short_repr("This string is too long", maxlen=5, repl="tl;dr")
    'tl;dr'

    could also use len(pickle.dumps()) to see if the object and not just its
    repr is too big to bother printing.
    """
    repr_out = repr(input_var)
    if maxlen:
        if len(repr_out) > maxlen:
            repr_out = repl

    return repr_out


def readable_call(funcname, args, kwargs, maxlen=None, repl="BIG_ARG"):
    """Print a human-readable expression for a function call: arguments and
    keyword arguments.

    >>> readable_call("test_func", [3,4,5], {'kw1': "this", 'kw2': "that"})
    "test_func(3, 4, 5, kw1='this', kw2='that')"

    >>> readable_call("test_func", [3,4,5], None)
    'test_func(3, 4, 5)'

    >>> readable_call("test_func", None, {'kw1': "this", 'kw2': "that"})
    "test_func(kw1='this', kw2='that')"

    >>> readable_call("test_func", None, None)
    'test_func()'

    Previous ways this was used:
    (execute_key, funcname, args, kwargs) = args_package
    print "%s: %s" % (execute_key, readable)
    TODO: deal with spurious commas
    """
    kwlist = []
    try:
        for item in kwargs:
            short_form = short_repr(kwargs[item], maxlen=maxlen, repl=repl)
            kwlist.append("%s=%s" % (item, short_form))
    except TypeError:
        kwlist = []

    arglist = []
    try:
        for item in args:
            short_form = short_repr(item, maxlen=maxlen, repl=repl)
            arglist.append(short_form)
    except TypeError:
        arglist = []

    outstring = "%s(" % funcname
    if arglist:
        outstring += ", ".join(arglist)

    if arglist and kwlist:
        outstring += ", "

    if kwlist:
        outstring += ", ".join(kwlist)

    return outstring + ")"


def _test_func(arg1, kwarg1="test"):
    print arg1, kwarg1


def func_exec(funcname, args, kwargs, printcall=True):
    """Execute a function based on its module.name with a given set of
    arguments and kwargs.

    >>> func_exec("_test_func", ["this works"], {"kwarg1": "ok"})
    _test_func('this works', kwarg1='ok')
    this works ok
    """
    if printcall:
        print readable_call(funcname, args, kwargs)

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

    return funcattr(*args, **kwargs)


if __name__ == "__main__":
    import doctest

    OPTIONFLAGS = (doctest.ELLIPSIS |
                   doctest.NORMALIZE_WHITESPACE)
    doctest.testmod(optionflags=OPTIONFLAGS)
