"""
Function pickling

REFERENCE: The python cookbook.
"""

import new, types, copy_reg, cPickle

def code_ctor(*args):
    """
    EXAMPLES:
    This indirectly tests this function.
        sage: def foo(a,b,c=10): return a+b+c
        sage: sage.misc.fpickle.reduce_code(foo.func_code)
        (<built-in function code_ctor>, ...)
        sage: unpickle_function(pickle_function(foo))
        <function foo at ...>
    """
    return new.code(*args)

def reduce_code(co):
    """
    EXAMPLES:
        sage: def foo(N): return N+1
        sage: sage.misc.fpickle.reduce_code(foo.func_code)
        (<built-in function code_ctor>, ...)
    """
    if co.co_freevars or co.co_cellvars:
        raise ValueError, "Cannot pickle code objects from closures"
    return code_ctor, (co.co_argcount, co.co_nlocals, co.co_stacksize,
                       co.co_flags, co.co_code, co.co_consts, co.co_names,
                       co.co_varnames, co.co_filename, co.co_name,
                       co.co_firstlineno, co.co_lnotab)

copy_reg.pickle(types.CodeType, reduce_code)

def pickle_function(func):
    """
    Pickle the Python function func.  This is not a normal pickle; you
    must use the unpickle_function method to unpickle the pickled
    function.

    NOTE: This does not work on all functions, but does work on
    'surprisingly' many functions.  In particular, it does not
    work on functions that includes nested functions.

    INPUT:
        func -- a Python function
    OUTPUT:
        a string

    EXAMPLES:
        sage: def f(N): return N+1
        ...
        sage: g = pickle_function(f)
        sage: h = unpickle_function(g)
        sage: h(10)
        11
    """
    return cPickle.dumps(func.func_code)

def unpickle_function(pickled):
    """
    Unpickle a pickled function.
    EXAMPLES:
        sage: def f(N,M): return N*M
        ...
        sage: unpickle_function(pickle_function(f))(3,5)
        15
    """
    recovered = cPickle.loads(pickled)
    return new.function(recovered, globals())



def call_pickled_function(fpargs):
    import sage.all
    from sage.misc.fpickle import unpickle_function
    (fp, (args, kwds)) = fpargs
    f = eval("unpickle_function(fp)", sage.all.__dict__, {'fp':fp})
    res = eval("f(*args, **kwds)",sage.all.__dict__, {'args':args, 'kwds':kwds, 'f':f})
    return ((args, kwds), res)
