import os
from sage.misc.misc import SAGE_ROOT
from sage.misc.cachefunc import cached_function
from sage.rings.integer import Integer
from sage.structure.sage_object import load

PATH = os.path.join(SAGE_ROOT,"data","cunningham_tables")

@cached_function
def cunningham_prime_factors():
    """
    List of all the prime numbers occuring in the so called Cunningham table.
    They occur in the factorization of numbers of type $b^n+1$ or $b^n-1$ with $b \in \{2,3,5,6,7,10,11,12\}$.
    Data from http://cage.ugent.be/~jdemeyer/cunningham/
    """
    file = os.path.join(PATH,"cunningham_prime_factors.sobj")
    if os.path.exists(file):
        return map(Integer,load(file))
    else:
        from warnings import warn
        warn("You might consider installing the optional package for factoring Cunningham numbers with the following command: ``sage -i cunningham_tables-1.0``")
        return []


