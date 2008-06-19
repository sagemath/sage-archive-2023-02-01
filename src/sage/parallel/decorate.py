"""
Decorate interface for parallel computation
"""

import cPickle
import os

from sage.structure.sage_object import save, load
from sage.misc.fpickle import pickle_function

from reference import parallel_iter as p_iter_reference
from dsage_p_iter import parallel_iter as p_iter_dsage

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

def hashargs(a):
    """
    INPUT:
        a -- object; could e list, tuple, etc.

    OUTPUT:
        nonnegative integer

    EXAMPLES:
        sage: sage.parallel.decorate.hashargs( ((2,), {'a':4}) )
        450573205
    """
    b = normalize_input(a)
    return abs(hash((b[0],tuple(b[1]))))

def easy_parallel_iter(f, inputs, p_iter=p_iter_reference):
    """
    Iterate over the inputs and pass them to f using p_iter.
    Returns an iterator over (input[i], f(...)).

    INPUT:
        f -- function
        inputs -- list
        p_iter -- function

    OUTPUT:
        generator

    EXAMPLES:
        sage: def f(N, m=2): return N*m
        sage: p_iter = sage.parallel.decorate.easy_parallel_iter
        sage: set_random_seed(0)
        sage: list(p_iter(f, [(2,3), ((5,), {'m':10}), 8], sage.parallel.reference.parallel_iter))
        [(8, 16), (((5,), {'m': 10}), 50), ((2, 3), 6)]
    """
    v = []
    argmap = {}

    for a in inputs:
        # we store the (args,kwds) version of each input argument,
        # along with a mapping back, which we pickle to make hashable.
        b = normalize_input(a)
        v.append(b)
        argmap[(b[0],tuple(b[1]))] = a

    # Now we iterate over the answers yielding the original inputs and values.
    for b, z in p_iter(f, v):
        yield argmap[(b[0],tuple(b[1]))], z

def parallel_eval(f, inputs, p_iter, dir=None, compress=True, threads=2):
    """
    INPUT:
        f      -- a function
        inputs -- a list of tuples, dicts, or objects
        p_iter -- parallel iterator function
        dir    -- string (optional); if given, make
                  a directory with a name made by pickling the
                  function f and inputs[i].  The pickled file
                  contains
                        (pf, inputs[i], f(inputs[i]))
                  where pf is the pickle of f.
        compress -- bool (default: True) if True and
                  dir is not None, then cached objects are
                  saved compressed

    OUTPUT:
        list of values of f on the tuples (all positional arguments),
        dicts (all keyword arguments), or on objects.  Also, files are
        created in dir if specified.

    EXAMPLES:
        sage: def f(N, m=2): return N*m
        sage: p_iter = sage.parallel.decorate.easy_parallel_iter
        sage: sage.parallel.decorate.parallel_eval(f, [(2,3), 5, (8,18), {'N':5,'m':3}], p_iter)
        [((2, 3), 6), (5, 10), ((8, 18), 144), ({'m': 3, 'N': 5}, 15)]
        sage: tmpdir = tmp_dir()
        sage: sage.parallel.decorate.parallel_eval(f, [{'N':5,'m':7}, (2,10)], p_iter, dir=tmpdir)
        [({'m': 7, 'N': 5}, 35), ((2, 10), 20)]
        sage: len(os.listdir(tmpdir))
        4
        sage: sage.parallel.decorate.parallel_eval(f, [(-1,3), {'N':5,'m':7}], p_iter)
        [((-1, 3), -3), ({'m': 7, 'N': 5}, 35)]
        sage: sage.parallel.decorate.parallel_eval(f, [1..4], p_iter)
        [(1, 2), (2, 4), (3, 6), (4, 8)]
        sage: sage.parallel.decorate.parallel_eval(f, (1..4), p_iter)
        [(1, 2), (2, 4), (3, 6), (4, 8)]
    """
    # Create the directory dir if it does not exist.
    # Note that we do *not* delete existing files in that directory.
    if dir is not None:
        dir = str(dir)
        if not os.path.exists(dir):
           os.makedirs(dir)

    # make inputs into a list.
    inputs = list(inputs)

    # v will be the list of (input,output) values that this function will return
    v = [None]*len(inputs)

    # We pickle the input function f
    pf = pickle_function(f)

    # Iterator over the parallel iterator, which can return values in
    # *any* order.
    to_compute = []  # non-cached values that we will need to compute
    for i, a in enumerate(inputs):
        # filename where we will write answer
        file = '%s/%s'%(dir, hashargs(a))

        # if it exists, just use cached value
        if dir is not None and os.path.exists(file+'.sobj'):
            pg, b, z = load(file + '.sobj')
            if b == a and pf == pg:
                v[i] = (a,z)
                continue
        # no cached value (or hash collision); we put
        # in the queue for running through p_iter
        to_compute.append(a)

    # Now for each value that wasn't already cached, we use p_iter to
    # find its value.  Note that the order of output from p_iter is random.
    # As we get results we optionally write them to a file.
    if len(to_compute) > 0:
        for a, z in easy_parallel_iter(f, to_compute, p_iter):
            v[inputs.index(a)] = (a,z)
            if dir is not None:
                file = '%s/%s'%(dir, hashargs(a))
                save((pf,a,z), file + '.sobj', compress=compress)
                open(file+'.txt','w').write(str(a))

    # Finally we return the answer list.
    return v

class parallel:
    """
    Create parallelizable functions.

    INPUT:
        p_iter -- parallel iterator function or string:
                   'reference'
                   'dsage'
        dir -- string (default: None)
        threads -- integer (default: 2)
        compress -- bool (default: True)

    EXAMPLES:
        sage: @parallel()
        ... def f(N): return N^2
        sage: f(10)
        100
        sage: P = parallel()
        sage: def g(N,M): return N*M
        sage: P(g)([(1,2), (4,5), (8,3)])
        [((1, 2), 2), ((4, 5), 20), ((8, 3), 24)]

    TESTS:
        sage: a = parallel(threads=10)
        sage: loads(dumps(a)) == a
        True
    """
    def __init__(self,
                 p_iter   = None,
                 dir       = None,
                 threads  = 2,
                 compress = True):
        """
        Create a parallel iterator decorator object.

        EXAMPLES:
            sage: @parallel(p_iter=None, dir=tmp_dir(), threads=4, compress=False)
            ... def f(N): return N^2
            sage: f([1,2,4])
            [(1, 1), (2, 4), (4, 16)]
        """
        if dir is not None:
            if not isinstance(dir, str):
                raise TypeError, "dir must be a string or None"
        # The default p_iter is currently the reference implementation.
        # This may change.
        if p_iter is None:
            self.p_iter = p_iter_reference
        elif p_iter == 'dsage':
            self.p_iter = p_iter_dsage
        elif p_iter == 'reference':
            self.p_iter = p_iter_reference
        else:
            if isinstance(p_iter, str):
                raise ValueError, "unknown iterator '%s'"%p_iter
            self.p_iter = p_iter
        try:
            t = int(threads)
        except TypeError:
            raise TypeError, "use @parallel(opts)\ndef foo(...):"
        if t != threads: # care that threads wasn't says 3.2 with int(3.2) --> 3.
            raise ValueError, "number of threads must be an integer"
        self.threads  = t
        self.dir      = dir
        self.compress = compress

    def __cmp__(self, right):
        """
        Compare self and right.

        EXAMPLES:
            sage: cmp(parallel(threads=2), parallel(threads=3))
            -1
            sage: cmp(parallel(threads=3), parallel(threads=2))
            1
        """
        if not isinstance(right, parallel):
            return cmp(type(self), type(right))
        return cmp((self.dir, self.threads, self.p_iter, self.compress),
                   (right.dir, right.threads, right.p_iter, right.compress))

    def __call__(self, f):
        """
        Create a function that wraps f and that when called with a
        list of inputs computes all values in parallel.

        INPUT:
            f -- Python callable object or function
        OUTPUT:
            decorated version of f

        EXAMPLES:
            sage: P = parallel(dir=tmp_dir())
            sage: def g(N,M): return N*M
            sage: P(g)([(1,2),(17,4)])
            [((1, 2), 2), ((17, 4), 68)]
            sage: len(os.listdir(P.dir))
            4
        """
        # Construct the wrapper parallel version of the function we're wrapping.
        # We may rework this so g is a class instance, which has the plus that
        # we can query g for how it works, etc.
        def g(*args, **kwds):
            if len(args) > 0 and isinstance(args[0], list):
                return parallel_eval(f, args[0], dir=self.dir, compress=self.compress,
                                     threads=self.threads, p_iter=self.p_iter)
            else:
                return f(*args, **kwds)
        return g


