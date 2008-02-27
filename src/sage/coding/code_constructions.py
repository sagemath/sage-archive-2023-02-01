r"""
AUTHOR:
    -- David Joyner (2007-05): initial version
    --    "         (2008-02): added cyclic codes, Hamming codes

This file contains contructions of error-correcting codes which are not
obtained from wrapping GUAVA functions. The GUAVA wrappers are in guava.py.
"""

import copy
import sage.modules.free_module as fm
import sage.modules.module as module
import sage.modules.free_module_element as fme
from sage.interfaces.all import gap
#from sage.misc.preparser import *
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.finite_field import FiniteField as GF
#from sage.groups.perm_gps.permgroup import *
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.misc.sage_eval import sage_eval
from sage.misc.misc import prod, add
from sage.misc.functional import log
from sage.rings.rational_field import QQ
from sage.structure.parent_gens import ParentWithGens
from linear_code import LinearCodeFromVectorSpace, LinearCode
from sage.modules.free_module import span
from sage.misc.functional import rank
from sage.schemes.generic.projective_space import ProjectiveSpace
from sage.structure.sequence import Sequence
from sage.rings.arith import gcd

############### utility functions ################

def permutation_action(g,v):
    """
    Returns permutation of rows g*v. Works on lists, matrices,
    sequences and vectors (by permuting coordinates). The code requires
    switching from i to i+1 (and back again) since the SymmetricGroup is,
    by convention, the symmetric group on the "letters"
    {1, 2, ..., n} (not {0, 1, ..., n-1}).

    EXAMPLES:
        sage: V = VectorSpace(GF(3),5)
        sage: v = V([0,1,2,0,1])
        sage: G = SymmetricGroup(5)
        sage: g = G([(1,2,3)])
        sage: permutation_action(g,v)
        (1, 2, 0, 0, 1)
        sage: g = G([()])
        sage: permutation_action(g,v)
        (0, 1, 2, 0, 1)
        sage: g = G([(1,2,3,4,5)])
        sage: permutation_action(g,v)
        (1, 2, 0, 1, 0)
        sage: L = Sequence([1,2,3,4,5])
        sage: permutation_action(g,L)
        [2, 3, 4, 5, 1]
        sage: MS = MatrixSpace(GF(3),3,7)
        sage: A = MS([[1,0,0,0,1,1,0],[0,1,0,1,0,1,0],[0,0,0,0,0,0,1]])
        sage: S5 = SymmetricGroup(5)
        sage: g = S5([(1,2,3)])
        sage: A; permutation_action(g,A)
        <BLANKLINE>
        [1 0 0 0 1 1 0]
        [0 1 0 1 0 1 0]
        [0 0 0 0 0 0 1]
        <BLANKLINE>
        [0 1 0 1 0 1 0]
        [0 0 0 0 0 0 1]
        [1 0 0 0 1 1 0]

    It also works on lists and is a "left action":

        sage: v = [0,1,2,0,1]
        sage: G = SymmetricGroup(5)
        sage: g = G([(1,2,3)])
        sage: gv = permutation_action(g,v); gv
        [1, 2, 0, 0, 1]
        sage: permutation_action(g,v) == g(v)
        True
        sage: h = G([(3,4)])
        sage: gv = permutation_action(g,v)
        sage: hgv = permutation_action(h,gv)
        sage: hgv == permutation_action(h*g,v)
        True

    AUTHOR: David Joyner, licensed under the GPL v2 or greater.
    """
    v_type_list = False
    if type(v) == list:
        v_type_list = True
        v = Sequence(v)
    V = v.parent()
    n = len(list(v))
    gv = []
    for i in range(n):
        gv.append(v[g(i+1)-1])
    if v_type_list:
        return gv
    return V(gv)

##################### main constructions #####################

def HammingCode(r,F):
    """
    Implements the Hamming codes.

    The $r^{th}$ Hamming code over $F=GF(q)$ is an $[n,k,d]$ code
    with length $n=(q^r-1)/(q-1)$, dimension $k=(q^r-1)/(q-1) - r$ and
    minimum distance $d=3$.
    The parity check matrix of a Hamming code has rows consisting of
    all nonzero vectors of length r in its columns, modulo a scalar factor
    so no parallel columns arise. A Hamming code is a single error-correcting
    code.

    INPUT:
        r -- an integer > 2
        F -- a finite field.

    OUTPUT:
        Returns the r-th q-ary Hamming code.

    EXAMPLES:
        sage: HammingCode(3,GF(2))
        Linear code of length 7, dimension 4 over Finite Field of size 2
        sage: C = HammingCode(3,GF(3)); C
        Linear code of length 13, dimension 10 over Finite Field of size 3
        sage: C.minimum_distance()
        3
        sage: C = HammingCode(3,GF(4,'a')); C
        Linear code of length 21, dimension 18 over Finite Field in a of size 2^2

    """
    q = F.order()
    n =  (q**r-1)/(q-1)
    k = n-r
    MS = MatrixSpace(F,n,r)
    X = ProjectiveSpace(r-1,F)
    PFn = [list(p) for p in X.point_set(F).points(F)]
    H = MS(PFn).transpose()
    Cd = LinearCode(H)
    return Cd.dual_code()

def ToricCode(P,F):
    r"""
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
        sage: C = ToricCode([[-2,-2],[-1,-2],[-1,-1],[-1,0],[0,-1],[0,0],[0,1],[1,-1],[1,0]],GF(5))
        sage: C
        Linear code of length 16, dimension 9 over Finite Field of size 5
        sage: C.minimum_distance()
        6
       sage: C = ToricCode([ [0,0],[1,1],[1,2],[1,3],[1,4],[2,1],[2,2],[2,3],[3,1],[3,2],[4,1]],GF(8,"a"))
        sage: C
        Linear code of length 49, dimension 11 over Finite Field in a of size 2^3

    This is in fact a [49,11,28] code over GF(8). If you type next
    \code{C.minimum_distance()} and wait overnight (!), you should get 28.


    AUTHOR: David Joyner (07-2006)

    REFERENCES:
        [J] D. Joyner, {\it Toric codes over finite fields}, Applicable Algebra in Engineering,
            Communication and Computing, 15, (2004), p. 63--79
    """
    from sage.combinat.all import Tuples
    mset = [x for x in F if x!=0]
    d = len(P[0])
    pts = Tuples(mset,d).list()
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
        sage: C.minimum_distance()
        8

    AUTHOR: David Joyner (2007-05)
    """
    v = [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]]
    V = span(GF(2), v)
    return LinearCodeFromVectorSpace(V)
    #C = BinaryGolayCode()
    #return C.extended_code()

def TernaryGolayCode():
    """
    TernaryGolayCode returns a ternary Golay code. This is a perfect
    [11,6,5] code.  It is also equivalenet to a cyclic code, with
    generator polynomial $g(x)=2+x^2+2x^3+x^4+x^5$.

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
    v = [[1, 0, 0, 0, 0, 0, 2, 0, 1, 2, 1, 2], [0, 1, 0, 0, 0, 0, 1, 2, 2, 2, 1, 0], [0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1], [0, 0, 0, 1, 0, 0, 1, 1, 0, 2, 2, 2], [0, 0, 0, 0, 1, 0, 2, 1, 2, 2, 0, 1], [0, 0, 0, 0, 0, 1, 0, 2, 1, 2, 2, 1]]
    V = span(GF(3), v)
    return LinearCodeFromVectorSpace(V)
    #C = TernaryGolayCode()
    #return C.extended_code()

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

def LinearCodeFromCheckMatrix(H):
    """
    A linear [n,k]-code C is uniquely determined by its generator
    matrix G and check matrix H. We have the following short
    exact sequence

    \[
    0 \rightarrow
    {\mathbf{F}}^k \stackrel{G}{\rightarrow}
    {\mathbf{F}}^n \stackrel{H}{\rightarrow}
    {\mathbf{F}}^{n-k} \rightarrow
    0.
    \end{equation}
    (``Short exact'' means (a) the arrow $G$ is injective,
    i.e., $G$ is a full-rank $k\times n$ matrix, (b) the arrow $H$ is
    surjective, and (c) ${\rm image}(G)={\rm kernel}(H)$.)

    EXAMPLES:
        sage: C = HammingCode(3,GF(2))
        sage: H = C.check_mat(); H
        <BLANKLINE>
        [1 0 0 1 1 0 1]
        [0 1 0 1 0 1 1]
        [0 0 1 1 1 1 0]
        sage: LinearCodeFromCheckMatrix(H) == C
        True
        sage: C = HammingCode(2,GF(3))
        sage: H = C.check_mat(); H
        <BLANKLINE>
        [1 0 2 2]
        [0 1 2 1]
        sage: LinearCodeFromCheckMatrix(H) == C
        True
        sage: C = RandomLinearCode(10,5,GF(4,"a"))
        sage: H = C.check_mat()
        sage: LinearCodeFromCheckMatrix(H) == C
        True

    """
    Cd = LinearCode(H)
    return Cd.dual_code()

def CyclicCodeFromGeneratingPolynomial(n,g,ignore=True):
    """
    If g is a polynomial over GF(q) which divides x^n-1 then this
    constructs the code "generated by g" (ie, the code associated with
    the principle ideal (g) in the ring GF(q)[x]/(x^n-1) in the usual way).

    The option "ignore" says to ignore the condition that
    (a) the characteristic of the base field does not divide the length
    (the usual assumtion in the theory of cyclic codes), and
    (b) $g$ must divide $x^n-1$. If ignore=True, instead of returning an error, a
    code generated by $gcd(x^n-1,g)$ is created.

    EXAMPLES:
        sage: P.<x> = PolynomialRing(GF(3),"x")
        sage: g = x-1
        sage: C = CyclicCodeFromGeneratingPolynomial(4,g); C
        Linear code of length 4, dimension 3 over Finite Field of size 3
        sage: P.<x> = PolynomialRing(GF(4,"a"),"x")
        sage: g = x^3+1
        sage: C = CyclicCodeFromGeneratingPolynomial(9,g); C
        Linear code of length 9, dimension 6 over Finite Field in a of size 2^2
        sage: P.<x> = PolynomialRing(GF(2),"x")
        sage: g = x^3+x+1
        sage: C = CyclicCodeFromGeneratingPolynomial(7,g); C
        Linear code of length 7, dimension 4 over Finite Field of size 2
        sage: C.gen_mat()
        <BLANKLINE>
        [1 1 0 1 0 0 0]
        [0 1 1 0 1 0 0]
        [0 0 1 1 0 1 0]
        [0 0 0 1 1 0 1]
        sage: g = x+1
        sage: C = CyclicCodeFromGeneratingPolynomial(4,g); C
        Linear code of length 4, dimension 3 over Finite Field of size 2
        sage: C.gen_mat()
        <BLANKLINE>
        [1 1 0 0]
        [0 1 1 0]
        [0 0 1 1]

    On the other hand, CyclicCodeFromPolynomial(4,x) will produce
    a ValueError including a traceback error message: "x must divide x^4 - 1".
    You will also get a ValueError if you type
        sage: P.<x> = PolynomialRing(GF(4,"a"),"x")
        sage: g = x^2+1

    followed by CyclicCodeFromGeneratingPolynomial(6,g). You will also
    get a ValueError if you type
        sage: P.<x> = PolynomialRing(GF(3),"x")
        sage: g = x^2-1
        sage: C = CyclicCodeFromGeneratingPolynomial(5,g); C
        Linear code of length 5, dimension 4 over Finite Field of size 3

    followed by C = CyclicCodeFromGeneratingPolynomial(5,g,False), with
    a traceback message including "x^2 + 2 must divide x^5 - 1".

    """
    P = g.parent()
    x = P.gen()
    F = g.base_ring()
    p = F.characteristic()
    if not(ignore) and p.divides(n):
        raise ValueError, 'The characteristic %s must not divide %s'%(p,n)
    if not(ignore) and not(g.divides(x**n-1)):
        raise ValueError, '%s must divide x^%s - 1'%(g,n)
    gn = gcd(x**n-1,g)
    d = gn.degree()
    coeffs = Sequence(gn.list())
    r1 = Sequence(coeffs+[0]*(n - d - 1))
    Sn = SymmetricGroup(n)
    s = Sn.gens()[0] ## assumes 1st gen of S_n is (1,2,...,n)
    rows = [permutation_action(s**(-i),r1) for i in range(n-d)]
    MS = MatrixSpace(F,n-d,n)
    return LinearCode(MS(rows))

def CyclicCodeFromCheckPolynomial(n,h,ignore=True):
    """
    If h is a polynomial over GF(q) which divides x^n-1 then this
    constructs the code "generated by g=(x^n-1)/h" (ie, the code associated with
    the principle ideal (g) in the ring GF(q)[x]/(x^n-1) in the usual way).
    The option "ignore" says to ignore the condition that the
    characteristic of the base field does nto divide the length
    (the usual assumtion in the theory of cyclic codes).

    EXAMPLES:
        sage: P.<x> = PolynomialRing(GF(3),"x")
        sage: C = CyclicCodeFromCheckPolynomial(4,x + 1); C
        Linear code of length 4, dimension 1 over Finite Field of size 3
        sage: C = CyclicCodeFromCheckPolynomial(4,x^3 + x^2 + x + 1); C
        Linear code of length 4, dimension 3 over Finite Field of size 3
        sage: C.gen_mat()
        [2 1 0 0]
        [0 2 1 0]
        [0 0 2 1]

    """
    P = h.parent()
    x = P.gen()
    d = h.degree()
    F = h.base_ring()
    p = F.characteristic()
    if not(ignore) and p.divides(n):
        raise ValueError, 'The characteristic %s must not divide %s'%(p,n)
    if not(h.divides(x**n-1)):
        raise ValueError, '%s must divide x^%s - 1'%(h,n)
    g = P((x**n-1)/h)
    return CyclicCodeFromGeneratingPolynomial(n,g)
