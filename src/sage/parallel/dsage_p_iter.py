"""
Implementation of parallel iterator using Dsage.
"""

import base64

from sage.misc.fpickle import pickle_function
from sage.structure.sage_object import loads

dsage_client = None
def parallel_iter(f, inputs, threads=2, blocking=True):
    """
    Parallel iterator implemented using DSage.

    INPUT:
        f -- a function
        inputs -- a list of tuples, dicts, or objects
        threads -- integer (default: 2)
        blocking -- bool (default: True)

    OUTPUT:
        iterator over 2-tuples (inputs[i], f(inputs[i])),
        where the order may be completely random

    EXAMPLES:
        sage: @parallel(p_iter='dsage')
        ... def foo(n,m):
        ...     return n+m
        sage: foo([(1,2), (5,10/3)])
        Going into testing mode...
        [((1, 2), 3), ((5, 10/3), 25/3)]
    """
    from sage.dsage.all import dsage
    import time

    fp = pickle_function(f)

    global dsage_client
    if dsage_client is None:
        dsage_client = dsage.start_all(workers=threads)
        time.sleep(1)  # HACK

    v = []
    for a in inputs:
        if isinstance(a, tuple):
            z = 'f(*a)'
        elif isinstance(a, dict):
            z = 'f(**a)'
        else:
            z = 'f(a)'
        # We do this base64 encoding, since dsage's DSAGE_RESULT just sucks.
        cmd = 'import base64;f=unpickle_function(fp);print "|"+base64.standard_b64encode(dumps(%s))+"|"'%z
        v.append((a, dsage_client.eval(cmd, user_vars={'fp':fp, 'a':a})))

    while True:
        if len(v) == 0:
            return
        for i in range(len(v)):
            a, j = v[i]
            #print j
            j.get_job()
            if j.status == 'completed':
                del v[i]
                out = j.output
                n = out.find('|')
                out = out[n+1:]
                n = out.find('|')
                out = out[:n]
                yield (a, loads(base64.standard_b64decode(out)))
                break
        time.sleep(0.5)   # HACK?


