"""
Cremona's mwrank C++ library
"""

__doc_exclude = []  # to include everything

from interface import (mwrank_EllipticCurve, mwrank_MordellWeil,
                       set_precision)


def mwrank_initprimes(filename, verb=False):
    """
    mwrank_initprimes(filename, verb=False):

    INPUT:


    -  ``filename`` - (string) the name of a file of
       primes

    -  ``verb`` - (bool: default False) verbose or not?


    EXAMPLES::

        sage: file= Sage_TMP + '/PRIMES'
        sage: open(file,'w').write(' '.join([str(p) for p in prime_range(10^6)]))
        sage: mwrank_initprimes(file, verb=False)
    """
    from mwrank import initprimes as mwrank_initprimes
    return mwrank_initprimes(filename, verb)



