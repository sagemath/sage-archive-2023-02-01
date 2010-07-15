"""
Decorate interface for parallel computation
"""

import os
import types

from sage.rings.all import Integer

from reference import parallel_iter as p_iter_reference
from use_fork import p_iter_fork
import multiprocessing_sage

def normalize_input(a):
    """
    Convert a to a pair (args, kwds) using some rules:

        - if already of that form, leave that way.
        - if a is a tuple make (a,{})
        - if a is a dict make (tuple([]),a)
        - otherwise make ((a,),{})

    INPUT:
        - ``a`` -- object
    OUTPUT:
        - ``args`` -- tuple
        - ``kwds`` -- dictionary

    EXAMPLES::

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
    """
    Create parallel decorated function.
    """
    def __init__(self, p_iter = 'fork', ncpus=None, **kwds):
        """
        EXAMPLES::

            sage: P = sage.parallel.decorate.Parallel(); P
            <sage.parallel.decorate.Parallel instance at 0x...>
        """
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
        """
        Create a function that wraps f and that when called with a
        list of inputs returns an iterator over pairs
             (x, f(x))
        in possibly random order.   Here x is replaced by its
        normalized form (args, kwds) using normalize_inputs.

        INPUT:
            - `f` -- Python callable object or function

        OUTPUT:
            decorated version of `f`

        EXAMPLES::

            sage: P = sage.parallel.decorate.Parallel()
            sage: def g(n,m): return n+m
            sage: h = P(g)          # indirect doctest
            sage: list(h([(2,3)]))
            [(((2, 3), {}), 5)]
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
    """
    This is a decorator that gives a function a parallel interface,
    allowing it to be called with a list of inputs, whose valuaes will
    be computed in parallel.

    .. warning::

         The parallel subprocesses will not have access to data
         created in pexpect interfaces.  This behavior with respect to
         pexpect interfaces is very important to keep in mind when
         setting up certain computations. It's the one big limitation
         of this decorator.

    INPUT:

        - ``p_iter`` -- parallel iterator function or string:
            - ``fork``            -- (default) use a new fork for each input
            - ``multiprocessing`` -- use multiprocessing library
            - ``reference``       -- use a fake serial reference implementation
        - ``ncpus`` -- integer, number of cpus
        - ``timeout`` -- number of seconds until task is killed (only supported by 'fork')

    .. warning::

         If you use anything but 'fork' above, then a whole new
         subprocess is spawned, so all your local state (variables,
         new functions, etc.,) is not available.

    EXAMPLES::

    We create a simple decoration for a simple function.  The nummber
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



###################################################################
# The @fork decorator -- run a function with no sideeffects in memory
# (so only side effects are on disk).
# We have both a function and a class below, so that the decorator
# can be used with or without options:
#   @fork
#   def f(...): ...
# and
#   @fork(...options...):
#   def f(...): ...
###################################################################

class Fork:
    """
    A fork decorator class.
    """
    def __init__(self, timeout=0):
        """
        INPUT:
            - ``timeout`` -- (default: 0) kills subrocess after this many
              seconds, or if timeout=0, do not kill the subprocess.

        EXAMPLES::

            sage: sage.parallel.decorate.Fork()
            <sage.parallel.decorate.Fork instance at 0x...>
            sage: sage.parallel.decorate.Fork(timeout=3)
            <sage.parallel.decorate.Fork instance at 0x...>

        """
        self.timeout = timeout

    def __call__(self, f):
        """
        INPUT:
            - f -- a function
        OUTPUT:
            - a decorated function

        EXAMPLES::

            sage: F = sage.parallel.decorate.Fork(timeout=3)
            sage: def g(n,m): return n+m
            sage: h = F(g)     # indirect doctest
            sage: h(2,3)
            5
        """
        P = Parallel(p_iter='fork', ncpus=1, timeout=self.timeout)
        g = P(f)
        def h(*args, **kwds):
            return list(g([(args, kwds)]))[0][1]
        return h

def fork(f=None, timeout=0):
    """
    Decorate a function so that when called it runs in a forked
    subprocess.  This means that it won't have any in-memory
    side-effects on the master Sage process.  The pexpect interfaces
    are all reset.

    INPUT:
        - ``f`` -- a function
        - ``timeout`` -- (default: 0) kills subrocess after this many
          seconds, or if timeout=0, do not kill the subprocess.

    .. warning::

        The forked subprocesses will not have access to data created
        in pexpect interfaces.  This behavior with respect to pexpect
        interfaces is very important to keep in mind when setting up
        certain computations. It's the one big limitation of this
        decorator.

    EXAMPLES::

    We create a function and run it with the fork decorator.  Note that
    it does not have a side effect.  Despite trying to change that global
    variable "a" below, it does not get changed.::

        sage: a = 5
        sage: @fork
        ... def g(n, m):
        ...       global a
        ...       a = 10
        ...       return factorial(n).ndigits() + m
        sage: g(5, m=5)
        8
        sage: a
        5

    We use fork to make sure that the function dies after 1 second no matter what::

        sage: @fork(timeout=1)
        ... def g(n, m): return factorial(n).ndigits() + m
        sage: g(5, m=5)
        8
        sage: g(10^7, m=5)
        Killing subprocess ... with input ((10000000,), {'m': 5}) which took too long
        'NO DATA (timed out)'

    We illustrate that pexpect interface state is not effected by
    forked functions (they get their own new pexpect interfaces!)::

        sage: gp.eval('a = 5')
        '5'
        sage: @fork()
        ... def g():
        ...       gp.eval('a = 10')
        ...       return gp.eval('a')
        sage: g()
        '10'
        sage: gp.eval('a')
        '5'

    We illustrate that the forked function has its own pexpect interface.::

        sage: gp.eval('a = 15')
        '15'
        sage: @fork()
        ... def g(): return gp.eval('a')
        sage: g()
        'a'

    We illustrate that segfaulting subprocesses are no trouble at all::

        sage: cython('def f(): print <char*>0')
        sage: @fork
        ... def g(): f()
        sage: g()
        'NO DATA'
    """
    F = Fork(timeout=timeout)
    return F(f) if f else F
