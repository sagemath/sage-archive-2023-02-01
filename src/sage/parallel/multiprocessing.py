from processing import Pool
from functools import partial
from sage.misc.fpickle import pickle_function

def pyprocessing(processes):
    """
    EXAMPLES:
        sage: from sage.parallel.multiprocessing import pyprocessing
        sage: p_iter = pyprocessing(4)
        sage: P = parallel(p_iter=p_iter)
        sage:
        sage: def f(x): return x+x
        ...
        sage: P(f)(range(10))

    """
    return partial(parallel_iter, processes)

def parallel_iter(processes, f, inputs):
    p = Pool(processes)
    fp = pickle_function(f)

    result = p.imapUnordered(call_pickled_function, [ (fp, t) for t in inputs ])
    for res in result:
        yield res

def call_pickled_function(fpargs):
    import sage.all
    from sage.misc.fpickle import unpickle_function
    (fp, (args, kwds)) = fpargs
    f = eval("unpickle_function(fp)", sage.all.__dict__, locals())
    res = eval("f(*args, **kwds)",sage.all.__dict__, locals())
    return ((args, kwds), res)
