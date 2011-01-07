r"""
Decorate interface for parallel computation
"""

import os
import types

from sage.rings.all import Integer

from reference import parallel_iter as p_iter_reference
from use_fork import p_iter_fork
import multiprocessing_sage

def normalize_input(a):
    r"""
    Convert a to a pair (args, kwds) using some rules:

        * if already of that form, leave that way.
        * if a is a tuple make (a,{})
        * if a is a dict make (tuple([]),a)
        * otherwise make ((a,),{})

    INPUT:

     - ``a`` -- object

    OUTPUT:

    args -- tuple
    kwds -- dictionary

    EXAMPLES:

        sage: sage.parallel.decorate.normalize_input( (2, {3:4}) )
        ((2, {3: 4}), {})
        sage: sage.parallel.decorate.normalize_input( (2,3) )
        ((2, 3), {})
        sage: sage.parallel.decorate.normalize_input( {3:4} )
        ((), {3: 4})
        sage: sage.parallel.decorate.normalize_input( 5 )
        ((5,), {})
    """
    if isinstance(a, tuple) and len(a) == 2 and isinstance(a[0],tuple) and isinstance(a[1],dict):
        return a
    elif isinstance(a, tuple):
        return (a, {})
    elif isinstance(a, dict):
        return (tuple([]), a)
    else:
        return ((a,), {})


class Parallel:
    r"""
    Create parallel decorated function.

    """
    def __init__(self, p_iter = 'fork', ncpus=None, **kwds):
        # The default p_iter is currently the reference implementation.
        # This may change.
        self.p_iter = None

        if isinstance(p_iter, (int, long, Integer)):
            p_iter, ncpus = 'fork', p_iter

        if ncpus is None:
            from ncpus import ncpus as compute_ncpus
            ncpus = compute_ncpus()

        if p_iter == 'fork':
            self.p_iter = p_iter_fork(ncpus, **kwds)
        elif p_iter == 'multiprocessing':
            self.p_iter = multiprocessing_sage.pyprocessing(ncpus)
        elif p_iter == 'reference':
            self.p_iter = p_iter_reference
        elif isinstance(p_iter, str):
            raise ValueError, "unknown iterator '%s'" % p_iter
        else:
            if self.p_iter is None:
                self.p_iter = p_iter

    def __call__(self, f):
        r"""
        Create a function that wraps f and that when called with a
        list of inputs returns an iterator over pairs
             (x, f(x))
        in possibly random order. Here x is replaced by its
        normalized form (args, kwds) using normalize_inputs.

        INPUT:

         - ``f`` -- Python callable object or function

        OUTPUT:

        decorated version of f

        EXAMPLES:

            sage: from sage.parallel.decorate import Parallel
            sage: p = Parallel()
            sage: f = x^2-1
            sage: p(f)
            <function g at ...
            sage: p(f)(x=5)
            24

        """
        # Construct the wrapper parallel version of the function we're wrapping.
        # We may rework this so g is a class instance, which has the plus that
        # we can query g for how it works, etc.
        def g(*args, **kwds):
            if len(args) > 0 and isinstance(args[0], (list, types.GeneratorType)):
                return self.p_iter(f, (normalize_input(a) for a in args[0]))
            else:
                return f(*args, **kwds)
        return g

def parallel(p_iter = 'fork', ncpus=None, **kwds):
    r"""
    This is a decorator that gives a function a parallel interface,
    allowing it to be called with a list of inputs, whose values will
    be computed in parallel.

    INPUT:

     - ``p_iter`` -- parallel iterator function or string:
            - 'fork'            -- (default) use a new fork for each input
            - 'multiprocessing' -- use multiprocessing library
            - 'reference'       -- use a fake serial reference implementation
     - ``ncpus`` -- integer, number of cpus
     - ``timeout`` -- number of seconds until task is killed (only supported by 'fork')

    EXAMPLES:

    We create a simple decoration for a simple function. The number
    of cpus is automatically detected::

        sage: @parallel
        ... def f(n): return n*n
        sage: f(10)
        100
        sage: sorted(list(f([1,2,3])))
        [(((1,), {}), 1), (((2,), {}), 4), (((3,), {}), 9)]

    We use exactly 2 cpus::

        sage: @parallel(2)
        ... def f(n): return n*n


    We create a decorator that uses 3 processes, and times out
    individual processes after 10 seconds::

        sage: @parallel(ncpus=3, timeout=10)
        ... def fac(n): return factor(2^n-1)
        sage: for X, Y in sorted(list(fac([101,119,151,197,209]))): print X,Y
        ((101,), {}) 7432339208719 * 341117531003194129
        ((119,), {}) 127 * 239 * 20231 * 131071 * 62983048367 * 131105292137
        ((151,), {}) 18121 * 55871 * 165799 * 2332951 * 7289088383388253664437433
        ((197,), {}) 7487 * 26828803997912886929710867041891989490486893845712448833
        ((209,), {}) 23 * 89 * 524287 * 94803416684681 * 1512348937147247 * 5346950541323960232319657

        sage: @parallel('multiprocessing')
        ... def f(N): return N^2
        sage: v = list(f([1,2,4])); v.sort(); v
        [(((1,), {}), 1), (((2,), {}), 4), (((4,), {}), 16)]
        sage: @parallel('reference')
        ... def f(N): return N^2
        sage: v = list(f([1,2,4])); v.sort(); v
        [(((1,), {}), 1), (((2,), {}), 4), (((4,), {}), 16)]
    """
    import types
    if isinstance(p_iter, types.FunctionType):
        return Parallel()(p_iter)
    return Parallel(p_iter, ncpus, **kwds)

