"""
Decorate interface for parallel computation
"""

import os
import types

from sage.rings.all import Integer

from reference import parallel_iter as p_iter_reference
from sage.dsage.interface import dsage_interface
import multiprocessing

def normalize_input(a):
    """
    Convert a to a pair (args, kwds) using some rules:

        * if already of that form, leave that way.
        * if a is a tuple make (a,{})
        * if a is a dict make (tuple([]),a)
        * otherwise make ((a,),{})

    INPUT:
        a -- object
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


class parallel:
    """
    Create paralleled functions.

    INPUT:
        p_iter -- parallel iterator function or string:
                  'multiprocessing' -- use multiprocessing (aka pyprocessing)
                  'reference'       -- use a fake serial reference
                                       implementation
                  DSage instance    -- use dsage
    """
    def __init__(self, p_iter = 'multiprocessing'):
        """
        Create a parallel iterator decorator object.

        EXAMPLES:
            sage: @parallel()
            ... def f(N): return N^2
            sage: v = list(f([1,2,4])); v.sort(); v
            [(((1,), {}), 1), (((2,), {}), 4), (((4,), {}), 16)]
            sage: @parallel('reference')
            ... def f(N): return N^2
            sage: v = list(f([1,2,4])); v.sort(); v
            [(((1,), {}), 1), (((2,), {}), 4), (((4,), {}), 16)]
        """
        # The default p_iter is currently the reference implementation.
        # This may change.
        if isinstance(p_iter, (int, long, Integer)):
            self.p_iter = multiprocessing.pyprocessing(p_iter)
        elif p_iter == 'multiprocessing':
            self.p_iter = multiprocessing.pyprocessing()
        elif isinstance(p_iter, dsage_interface.BlockingDSage):
            self.p_iter = p_iter.parallel_iter
        elif p_iter == 'reference':
            self.p_iter = p_iter_reference
        else:
            if isinstance(p_iter, str):
                raise ValueError, "unknown iterator '%s'" % p_iter
            self.p_iter = p_iter

    def __call__(self, f):
        """
        Create a function that wraps f and that when called with a
        list of inputs returns an iterator over pairs
             (x, f(x))
        in possibly random order.   Here x is replaced by its
        normalized form (args, kwds) using normalize_inputs.

        INPUT:
            f -- Python callable object or function
        OUTPUT:
            decorated version of f

        EXAMPLES:

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

