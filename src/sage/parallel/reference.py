"""
Reference Parallel Primitives.

These are reference implementations of basic parallel
primitives. These are not actually parallel, but work the same way.
They are good for testing.
"""

def parallel_iter(f, inputs):
    """
    Reference parallel iterator implementation.

    INPUT:
        f -- a Python function that can be pickled using
             the pickle_function command.
        inputs -- a list of pickleable pairs (args, kwds), where args
             is a tuple and kwds is a dictionary.

    OUTPUT:
        iterator over 2-tuples (inputs[i], f(inputs[i])),
        where the order may be completely random

    EXAMPLES:
        sage: def f(N,M=10): return N*M
        sage: for a, val in sage.parallel.reference.parallel_iter(f, [((2,3),{}),  (tuple([]), {'N':3,'M':5}), ((2,),{})]):
        ...       print a, val
        ((2, 3), {}) 6
        ((), {'M': 5, 'N': 3}) 15
        ((2,), {}) 20
    """
    for args, kwds in inputs:
        yield ((args, kwds), f(*args, **kwds))
