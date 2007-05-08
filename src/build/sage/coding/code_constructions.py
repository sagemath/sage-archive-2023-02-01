r"""
AUTHOR:
    -- David Joyner (2007-05): initial version

This file contains contructions of error-correcting codes which are not
obtained from wrapping GUAVA functions. The GUAVA wrappers are in guava.py.
"""

import copy
import sage.modules.free_module as fm
import sage.modules.module as module
import sage.modules.free_module_element as fme
from sage.databases.lincodes import linear_code_bound
from sage.interfaces.all import gap
from sage.misc.preparser import *
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.finite_field import *
from sage.groups.perm_gps.permgroup import *
from sage.misc.sage_eval import sage_eval
from sage.misc.misc import prod, add
from sage.misc.functional import log
from sage.rings.rational_field import QQ
from sage.structure.parent_gens import ParentWithGens
from linear_code import *
from sage.modules.free_module import span
from sage.misc.functional import rank

def ToricCode(P,F):
    """
    Let $P$ denote a list of lattice points in $\Z^d$ and let $T$ denote the
    set of all points in $(F^x )^d$ (ordered in some fixed way). Put $n=|T|$
    and let $k$ denote the dimension of the vector space of functions
    $V = Span \{x^e \ |\ e \\in P\}$. The associated {\it toric code} $C$ is the
    evaluation code which is the image of the evaluation map
    $$
                        eval_T : V \\rightarrow F^n,
    $$
    where $x^e$ is the multi-index notation ($x=(x_1,...,x_d)$, $e=(e_1,...,e_d)$, and
    $x^e = x_1^{e_1}...x_d^{e_d}$), where $eval_T (f(x)) = (f(t_1),...,f(t_n))$, and
    where $T=\{t_1,...,t_n\}$. This function returns the toric codes discussed in [J].

    INPUT:
        P -- all the integer lattice points in a polytope defining the toric variety.
        F -- a finite field.

    OUTPUT:
        Returns toric code with length n = , dimension k over field F.

    EXAMPLES:
        sage: C = ToricCode([[0,0],[1,0],[2,0],[0,1],[1,1]],GF(7))
        sage: C
        Linear code of length 36, dimension 5 over Finite Field of size 7
        sage: C.minimum_distance()
        24
        sage: C.minimum_distance_upper_bound()  # optional (requires internet)
        28
        sage: C = ToricCode([[-2,-2],[-1,-2],[-1,-1],[-1,0],[0,-1],[0,0],[0,1],[1,-1],[1,0]],GF(5))
        sage: C
        Linear code of length 16, dimension 9 over Finite Field of size 5
        sage: C.minimum_distance()
        6
        sage: C.minimum_distance_upper_bound()   # optional -- uses internet
        6
        sage: C = ToricCode([ [0,0],[1,1],[1,2],[1,3],[1,4],[2,1],[2,2],[2,3],[3,1],[3,2],[4,1]],GF(8,"a"))
        sage: C
        Linear code of length 49, dimension 11 over Finite Field in a of size 2^3
        sage: C.minimum_distance()  ## long time -- very time consuming
        28
        sage: print linear_code_bound(8,49,11)[0]    # optional -- uses internet
        28


    AUTHOR: David Joyner (07-2006)

    REFERENCES:
        [J] D. Joyner, {\it Toric codes over finite fields}, Applicable Algebra in Engineering,
            Communication and Computing, 15, (2004), p. 63--79
    """
    from sage.combinat.combinat import tuples
    mset = [x for x in F if x!=0]
    d = len(P[0])
    pts = tuples(mset,d)
    n = len(pts) ## (q-1)^d
    k = len(P)
    e = P[0]
    B = []
    for e in P:
       tmpvar = [prod([t[i]**e[i] for i in range(d)]) for t in pts]
       B.append(tmpvar)
    ## now B0 *should* be a full rank matrix
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(B))

def TrivialCode(F,n):
    MS = MatrixSpace(F,1,n)
    return LinearCode(MS(0))

def BinaryGolayCode():
    """
    BinaryGolayCode() returns a binary Golay code. This is a perfect [23,12,7] code.
    It is also (equivalent to) a cyclic code, with generator polynomial
    $g(x)=1+x^2+x^4+x^5+x^6+x^{10}+x^{11}$.
    Extending it yields the extended Golay code (see ExtendedBinaryGolayCode).

    EXAMPLE:
        sage: C = BinaryGolayCode()
        sage: C
        Linear code of length 23, dimension 12 over Finite Field of size 2
        sage: C.minimum_distance()               # long time
        7

    AUTHOR: David Joyner (2007-05)
    """
    F = GF(2)
    M = [[1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
          [0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
          [0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
          [0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],\
          [0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0],\
          [0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],\
          [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],\
          [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0],\
          [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0],\
          [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0],\
          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0],\
          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1]]
    #MS = MatrixSpace(F,12,23)
    #V = VectorSpace(F,23)
    V = span(F, M)
    return LinearCodeFromVectorSpace(V)

def ExtendedBinaryGolayCode():
    """
    ExtendedBinaryGolayCode() returns the extended binary Golay code. This is a
    perfect [24,12,8] code. This code is self-dual.

    EXAMPLES:
        sage: C = ExtendedBinaryGolayCode()
        sage: C
        Linear code of length 24, dimension 12 over Finite Field of size 2
        sage: C.minimum_distance()               # long time
        8

    AUTHOR: David Joyner (2007-05)
    """
    C = BinaryGolayCode()
    return C.extended_code()

def TernaryGolayCode():
    """
    TernaryGolayCode returns a ternary Golay code. This is a perfect [11,6,5] code.
    It is also equivalenet to a cyclic code, with generator polynomial $g(x)=2+x^2+2x^3+x^4+x^5$.

    EXAMPLES:
        sage: C = TernaryGolayCode()
        sage: C
        Linear code of length 11, dimension 6 over Finite Field of size 3
        sage: C.minimum_distance()
        5

    AUTHOR: David Joyner (2007--5)
    """
    F = GF(3)
    M = [[2, 0, 1, 2, 1, 1, 0, 0, 0, 0, 0],\
         [0, 2, 0, 1, 2, 1, 1, 0, 0, 0, 0],\
         [0, 0, 2, 0, 1, 2, 1, 1, 0, 0, 0],\
         [0, 0, 0, 2, 0, 1, 2, 1, 1, 0, 0],\
         [0, 0, 0, 0, 2, 0, 1, 2, 1, 1, 0],\
         [0, 0, 0, 0, 0, 2, 0, 1, 2, 1, 1]]
    V = span(F, M)
    return LinearCodeFromVectorSpace(V)

def ExtendedTernaryGolayCode():
    """
    ExtendedTernaryGolayCode returns a ternary Golay code. This is a self-dual perfect [12,6,6] code.

    EXAMPLES:
        sage: C = ExtendedTernaryGolayCode()
        sage: C
        Linear code of length 12, dimension 6 over Finite Field of size 3
        sage: C.minimum_distance()
        6

    AUTHOR: David Joyner (11-2005)
    """
    C = TernaryGolayCode()
    return C.extended_code()

def RandomLinearCode(n,k,F):
    """
    The method used is to first construct a $k \\times n$ matrix using SAGE's
    random_element method for the MatrixSpace class. The construction is
    probabilistic but should only fail extremely rarely.

    INPUT:
        Integers n,k, with n>k>1, and a finite field F

    OUTPUT:
        Returns a "random" linear code with length n, dimension k over field F.

    EXAMPLES:
        sage: C = RandomLinearCode(30,15,GF(2))
        sage: C                                        # random output
        Linear code of length 30, dimension 15 over Finite Field of size 2
        sage: C = RandomLinearCode(10,5,GF(4,'a'))
        sage: C                                       # random output
        Linear code of length 10, dimension 5 over Finite Field in a of size 2^2

    AUTHOR: David Joyner (2007-05)
    """
    MS = MatrixSpace(F,k,n)
    for i in range(10):
        G = MS.random_element()
        if G.rank() == k:
            V = span(F, G.rows())
            return LinearCodeFromVectorSpace(V)
    return "fail"
