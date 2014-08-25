r"""
Difference Matrices

This module gathers code related to difference matrices. One can build those
objects (or know if they can be built) with :func:`difference_matrix`::

    sage: G,DM = designs.difference_matrix(9,5,1)

Functions
---------
"""
from sage.misc.unknown import Unknown
from sage.categories.sets_cat import EmptySetError
from sage.rings.finite_rings.constructor import FiniteField
from sage.rings.arith import is_prime_power
from designs_pyx import is_difference_matrix
from database import DM as DM_constructions

def difference_matrix(g,k,lmbda=1,existence=False,check=True):
    r"""
    Return a `(g,k,\lambda)`-difference matrix

    A matrix `M` is a `(g,k,\lambda)`-difference matrix if its entries are
    element of a group `G` of cardinality `g`, and if for any two rows `R,R'` of
    `M` and `x\in G` there are exactly `\lambda` values `i` such that
    `R_i-R'_i=x`.

    INPUT:

    - ``k`` -- (integer) number of columns. If ``k=None`` it is set to the
      largest value available.

    - ``g`` -- (integer) cardinality of the group `G`

    - ``lmbda`` -- (integer; default: 1) -- number of times each element of `G`
      appears as a difference.

    - ``check`` -- (boolean) Whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to
      ``True`` by default.

    - ``existence`` (boolean) -- instead of building the design, return:

        - ``True`` -- meaning that Sage knows how to build the design

        - ``Unknown`` -- meaning that Sage does not know how to build the
          design, but that the design may exist (see :mod:`sage.misc.unknown`).

        - ``False`` -- meaning that the design does not exist.

      .. NOTE::

          When ``k=None`` and ``existence=True`` the function returns an
          integer, i.e. the largest `k` such that we can build a
          `(g,k,\lambda)`-DM.

    EXAMPLES::

        sage: G,M = designs.difference_matrix(25,10); G
        Finite Field in x of size 5^2
        sage: designs.difference_matrix(993,None,existence=1)
        32

    TESTS::

        sage: designs.difference_matrix(10,12,1,existence=True)
        False
        sage: designs.difference_matrix(10,12,1)
        Traceback (most recent call last):
        ...
        EmptySetError: No (10,12,1)-Difference Matrix exists as k(=12)>g(=10)
        sage: designs.difference_matrix(10,9,1,existence=True)
        Unknown
        sage: designs.difference_matrix(10,9,1)
        Traceback (most recent call last):
        ...
        NotImplementedError: I don't know how to build a (10,9,1)-Difference Matrix!
    """

    if lmbda == 1 and k is not None and k>g:
        if existence:
            return False
        raise EmptySetError("No ({},{},{})-Difference Matrix exists as k(={})>g(={})".format(g,k,lmbda,k,g))

    # Prime powers
    elif lmbda == 1 and is_prime_power(g):
        if k is None:
            if existence:
                return g
            else:
                k = g
        elif existence:
            return True
        F       = FiniteField(g,'x')
        F_set   = list(F)
        F_k_set = F_set[:k]

        G = F
        M = [[x*y for y in F_k_set] for x in F_set]

    # From the database
    elif (g,lmbda) in DM_constructions and (k is None or DM_constructions[g,lmbda][0]>=k):
        if k is None:
            k = DM_constructions[g,lmbda][0]
            if existence:
                return k
        elif existence:
            return True
        _,f = DM_constructions[g,lmbda]
        G, M = f()
        M = [R[:k] for R in M]

    else:
        if existence:
            return Unknown
        raise NotImplementedError("I don't know how to build a ({},{},{})-Difference Matrix!".format(g,k,lmbda))

    if check:
        assert is_difference_matrix(M,G,k,lmbda,1), "Sage built something which is not a ({},{},{})-DM!".format(g,k,lmbda)

    return G,M
