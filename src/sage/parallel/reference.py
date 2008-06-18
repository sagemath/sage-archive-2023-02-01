"""
Reference Parallel Primitives.

These are reference implementations of basic parallel
primitives. These are not actually parallel, but work the same way.
They are good for testing.
"""

def parallel_iter(f, inputs, threads=2, blocking=True):
    """
    Reference parallel iterator.

    INPUT:
        f -- a function
        inputs -- a list of tuples, dicts, or objects
        threads -- integer (default: 2)
        blocking -- bool (default: True)

    OUTPUT:
        iterator over 2-tuples (inputs[i], f(inputs[i])),
        where the order may be completely random

    EXAMPLES:
        sage: def f(N,M=10): return N*M
        sage: for a, val in sage.parallel.reference.parallel_iter(f, [(2,3), {'N':3,'M':5}, 2], threads=3):
        ...       print a, val
        (2, 3) 6
        {'M': 5, 'N': 3} 15
        2 20
    """
    for a in inputs:
        if isinstance(a, tuple):
            z = f(*a)
        elif isinstance(a, dict):
            z = f(**a)
        else:
            z = f(a)
        yield a, z

