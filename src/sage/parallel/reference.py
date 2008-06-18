import cPickle
import os
from sage.structure.sage_object import save, load
from sage.misc.fpickle import pickle_function

def parallel_eval(f, inputs, db=None, compress=True, threads=2):
    """
    INPUT:
        f -- a function
        inputs -- a list of tuples, dicts, or objects
        db -- string (optional); if given, make
              a directory with a name made by pickling the
              function f and inputs[i].  The pickled file
              contains
                    (pf, inputs[i], f(inputs[i]))
              where pf is the pickle of f.
        compress -- bool (default: True) if True and
              db is not None, then cached objects are
              saved compressed

    OUTPUT:
        list of values of f on the tuples (all positional
        arguments), dicts (all keyword arguments), or
        on objects.

    EXAMPLES:
        sage: def f(N, m=2): return N*m
        sage: parallel_eval(f, [(2,3), 5, (8,18), {'N':5,'m':3}])
        [((2, 3), 6), (5, 10), ((8, 18), 144), ({'m': 3, 'N': 5}, 15)]
        sage: tmpdir = tmp_dir()
        sage: print parallel_eval(f, [{'N':5,'m':7}, (2,10)], tmpdir)
        sage: os.listdir(tmpdir)
        [({'m': 7, 'N': 5}, 35), ((2, 10), 20)]
        ['0.sobj', '1.sobj']
        sage: parallel_eval(f, [(-1,3), {'N':5,'m':7}])
        [((-1, 3), -3), ({'m': 7, 'N': 5}, 35)]
        sage: parallel_eval(f, [1,2,3,4])
        [(1, 2), (2, 4), (3, 6), (4, 8)]
    """
    if db is not None:
        db = str(db)
        if not os.path.exists(db):
           os.makedirs(db)
    v = []
    pf = pickle_function(f)
    for i, a in enumerate(inputs):
        n = abs(hash((pf, cPickle.dumps(a, 2))))
        file = '%s/%s'%(db, n)
        file_sobj = file + '.sobj'
        if os.path.exists(file_sobj):
            pg, b, z = load(file_sobj)
            if b == a and pf == pg:
                v.append((a,z))
                continue
        if isinstance(a, tuple):
            z = f(*a)
        elif isinstance(a, dict):
            z = f(**a)
        else:
            z = f(a)
        if db is not None:
            save((pf,a,z), file_sobj, compress=compress)
            open(file + '.txt','w').write(str(a))
        v.append((a,z))
    return v

