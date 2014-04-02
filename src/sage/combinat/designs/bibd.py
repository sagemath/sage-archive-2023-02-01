r"""
Balanced Incomplete Block Designs (BIBD)

This module implements two constructions of Balanced Incomplete Block Designs:

* Steiner Triple Systems, i.e. `(v,3,1)`-BIBD.
* `K_4`-decompositions of `K_v`, i.e. `(v,4,1)`-BIBD.

These BIBD can be obtained through the :meth:`BalancedIncompleteBlockDesign`
method, available in Sage as ``designs.BalancedIncompleteBlockDesign``.

EXAMPLES::

    sage: designs.BalancedIncompleteBlockDesign(7,3)
    Incidence structure with 7 points and 7 blocks
    sage: designs.BalancedIncompleteBlockDesign(7,3).blocks()
    [[0, 1, 3], [0, 2, 4], [0, 5, 6], [1, 2, 6], [1, 4, 5], [2, 3, 5], [3, 4, 6]]
    sage: designs.BalancedIncompleteBlockDesign(13,4).blocks()
    [[0, 1, 2, 12], [0, 3, 6, 9], [0, 4, 8, 11], [0, 5, 7, 10], [1, 3, 8, 10],
     [1, 4, 7, 9], [1, 5, 6, 11], [2, 3, 7, 11], [2, 4, 6, 10], [2, 5, 8, 9],
     [3, 4, 5, 12], [6, 7, 8, 12], [9, 10, 11, 12]]

`K_4`-decompositions of `K_v`
-----------------------------

Decompositions of `K_v` into `K_4` (i.e. `(v,4,1)`-BIBD) are built following
Douglas Stinson's construction as presented in [Stinson2004]_ page 167. It is
based upon the construction of `(v\{4,5,8,9,12\})`-PBD (see the doc of
:meth:`PBD_4_5_8_9_12`), knowing that a `(v\{4,5,8,9,12\})`-PBD on `v` points
can always be transformed into a `((k-1)v+1,4,1)`-BIBD, which covers all
possible cases of `(v,4,1)`-BIBD.

Functions
---------
"""
from design_catalog import transversal_design
from block_design import BlockDesign
from sage.rings.arith import binomial
from sage.rings.arith import is_prime_power

def BalancedIncompleteBlockDesign(v,k,use_LJCR=False):
    r"""
    Returns a BIBD of parameters `v,k`.

    A Balanced Incomplete Block Design of parameters `v,k` is a collection
    `\mathcal C` of `k`-subsets of `V=\{0,\dots,v-1\}` such that for any two
    distinct elements `x,y\in V` there is a unique element `S\in \mathcal C`
    such that `x,y\in S`.

    For more information on BIBD, see the
    :wikipedia:`corresponding Wikipedia entry <Block_design#Definition_of_a_BIBD_.28or_2-design.29>`.

    INPUT:

    - ``v,k`` (integers)

    - ``use_LJCR`` (boolean) -- whether to query the La Jolla Covering
      Repository for the design when Sage does not know how to build it (see
      :meth:`~sage.combinat.designs.covering_design.best_known_covering_design_www`). This
      requires internet.

    .. SEEALSO::

        * :meth:`steiner_triple_system`
        * :meth:`v_4_1_bibd`

    TODO:

        * Implement `(v,5,1)`-BIBD using `this text <http://www.argilo.net/files/bibd.pdf>`_.
        * Implement other constructions from the Handbook of Combinatorial
          Designs.

    EXAMPLES::

        sage: designs.BalancedIncompleteBlockDesign(7,3).blocks()
        [[0, 1, 3], [0, 2, 4], [0, 5, 6], [1, 2, 6], [1, 4, 5], [2, 3, 5], [3, 4, 6]]
        sage: B = designs.BalancedIncompleteBlockDesign(21,5, use_LJCR=True) # optional - internet
        sage: B                                                              # optional - internet
        Incidence structure with 21 points and 21 blocks
        sage: B.blocks()                                                     # optional - internet
        [[0, 1, 2, 3, 20], [0, 4, 8, 12, 16], [0, 5, 10, 15, 19],
         [0, 6, 11, 13, 17], [0, 7, 9, 14, 18], [1, 4, 11, 14, 19],
         [1, 5, 9, 13, 16], [1, 6, 8, 15, 18], [1, 7, 10, 12, 17],
         [2, 4, 9, 15, 17], [2, 5, 11, 12, 18], [2, 6, 10, 14, 16],
         [2, 7, 8, 13, 19], [3, 4, 10, 13, 18], [3, 5, 8, 14, 17],
         [3, 6, 9, 12, 19], [3, 7, 11, 15, 16], [4, 5, 6, 7, 20],
         [8, 9, 10, 11, 20], [12, 13, 14, 15, 20], [16, 17, 18, 19, 20]]
        sage: designs.BalancedIncompleteBlockDesign(20,5, use_LJCR=True) # optional - internet
        Traceback (most recent call last):
        ...
        ValueError: No such design exists !

    TESTS:

    A BIBD from a Finite Projective Plane::

        sage: _ = designs.BalancedIncompleteBlockDesign(21,5)
    """
    if ((binomial(v,2)%binomial(k,2) != 0) or
        (v-1)%(k-1) != 0):
        raise ValueError("No such design exists !")

    if k == 2:
        from itertools import combinations
        return BlockDesign(v, combinations(range(v),2), test = False)
    if k == 3:
        return steiner_triple_system(v)
    if k == 4:
        return BlockDesign(v, v_4_1_BIBD(v), test = False)
    if v == (k-1)**2+k and is_prime_power(k-1):
        from block_design import ProjectivePlaneDesign
        return ProjectivePlaneDesign(k-1)
    if use_LJCR:
        from covering_design import best_known_covering_design_www
        B = best_known_covering_design_www(v,k,2)

        # Is it a BIBD or just a good covering ?
        expected_n_of_blocks = binomial(v,2)/binomial(k,2)
        if B.low_bd() > expected_n_of_blocks:
            raise ValueError("No such design exists !")
        B = B.incidence_structure()
        if len(B.blcks) == expected_n_of_blocks:
            return B

    raise ValueError("I don't know how to build this design.")

def steiner_triple_system(n):
    r"""
    Returns a Steiner Triple System

    A Steiner Triple System (STS) of a set `\{0,...,n-1\}`
    is a family `S` of 3-sets such that for any `i \not = j`
    there exists exactly one set of `S` in which they are
    both contained.

    It can alternatively be thought of as a factorization of
    the complete graph `K_n` with triangles.

    A Steiner Triple System of a `n`-set exists if and only if
    `n \equiv 1 \pmod 6` or `n \equiv 3 \pmod 6`, in which case
    one can be found through Bose's and Skolem's constructions,
    respectively [AndHonk97]_.

    INPUT:

    - ``n`` returns a Steiner Triple System of `\{0,...,n-1\}`

    EXAMPLE:

    A Steiner Triple System on `9` elements ::

        sage: sts = designs.steiner_triple_system(9)
        sage: sts
        Incidence structure with 9 points and 12 blocks
        sage: list(sts)
        [[0, 1, 5], [0, 2, 4], [0, 3, 6], [0, 7, 8], [1, 2, 3], [1, 4, 7], [1, 6, 8], [2, 5, 8], [2, 6, 7], [3, 4, 8], [3, 5, 7], [4, 5, 6]]

    As any pair of vertices is covered once, its parameters are ::

        sage: sts.parameters(t=2)
        (2, 9, 3, 1)

    An exception is raised for invalid values of ``n`` ::

        sage: designs.steiner_triple_system(10)
        Traceback (most recent call last):
        ...
        ValueError: Steiner triple systems only exist for n = 1 mod 6 or n = 3 mod 6

    REFERENCE:

    .. [AndHonk97] A short course in Combinatorial Designs,
      Ian Anderson, Iiro Honkala,
      Internet Editions, Spring 1997,
      http://www.utu.fi/~honkala/designs.ps
    """

    name = "Steiner Triple System on "+str(n)+" elements"

    if n%6 == 3:
        t = (n-3)/6
        Z = range(2*t+1)

        T = lambda (x,y) : x + (2*t+1)*y

        sts = [[(i,0),(i,1),(i,2)] for i in Z] + \
            [[(i,k),(j,k),(((t+1)*(i+j)) % (2*t+1),(k+1)%3)] for k in range(3) for i in Z for j in Z if i != j]

    elif n%6 == 1:

        t = (n-1)/6
        N = range(2*t)
        T = lambda (x,y) : x+y*t*2 if (x,y) != (-1,-1) else n-1

        L1 = lambda i,j : (i+j) % (int((n-1)/3))
        L = lambda i,j : L1(i,j)/2 if L1(i,j)%2 == 0 else t+(L1(i,j)-1)/2

        sts = [[(i,0),(i,1),(i,2)] for i in range(t)] + \
            [[(-1,-1),(i,k),(i-t,(k+1) % 3)] for i in range(t,2*t) for k in [0,1,2]] + \
            [[(i,k),(j,k),(L(i,j),(k+1) % 3)] for k in [0,1,2] for i in N for j in N if i < j]

    else:
        raise ValueError("Steiner triple systems only exist for n = 1 mod 6 or n = 3 mod 6")

    from sage.sets.set import Set
    sts = Set(map(lambda x: Set(map(T,x)),sts))

    return BlockDesign(n, sts, name=name)

def v_4_1_BIBD(v, check=True):
    r"""
    Returns a `(v,4,1)`-BIBD.

    A `(v,4,1)`-BIBD is an edge-decomposition of the complete graph `K_v` into
    copies of `K_4`. For more information, see
    :meth:`BalancedIncompleteBlockDesign`. It exists if and only if `v\equiv 1,4
    \pmod {12}`.

    See page 167 of [Stinson2004]_ for the construction details.

    .. SEEALSO::

        * :meth:`BalancedIncompleteBlockDesign`

    INPUT:

    - ``v`` (integer) -- number of points.

    - ``check`` (boolean) -- whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to ``True``
      by default.

    EXAMPLES::

        sage: from sage.combinat.designs.bibd import v_4_1_BIBD  # long time
        sage: for n in range(13,100):                            # long time
        ....:    if n%12 in [1,4]:                               # long time
        ....:       _ = v_4_1_BIBD(n, check = True)              # long time
    """
    from sage.rings.finite_rings.constructor import FiniteField
    k = 4
    if v == 0:
        return []
    if v <= 12 or v%12 not in [1,4]:
        raise ValueError("A K_4-decomposition of K_v exists iif v=2,4 mod 12, v>12 or v==0")

    # Step 1. Base cases.
    if v == 13:
        from block_design import ProjectivePlaneDesign
        return ProjectivePlaneDesign(3).blocks()
    if v == 16:
        from block_design import AffineGeometryDesign
        return AffineGeometryDesign(2,1,FiniteField(4,'x')).blocks()
    if v == 25:
        return [[0, 1, 17, 22], [0, 2, 11, 21], [0, 3, 15, 18], [0, 4, 7, 13],
                [0, 5, 12, 14], [0, 6, 19, 23], [0, 8, 16, 24], [0, 9, 10, 20],
                [1, 2, 3, 4], [1, 5, 6, 7], [1, 8, 12, 15], [1, 9, 13, 16],
                [1, 10, 11, 14], [1, 18, 20, 23], [1, 19, 21, 24], [2, 5, 15, 24],
                [2, 6, 9, 17], [2, 7, 14, 18], [2, 8, 22, 23], [2, 10, 12, 13],
                [2, 16, 19, 20], [3, 5, 16, 22], [3, 6, 11, 20], [3, 7, 12, 19],
                [3, 8, 9, 14], [3, 10, 17, 24], [3, 13, 21, 23], [4, 5, 10, 23],
                [4, 6, 8, 21], [4, 9, 18, 24], [4, 11, 15, 16], [4, 12, 17, 20],
                [4, 14, 19, 22], [5, 8, 13, 20], [5, 9, 11, 19], [5, 17, 18, 21],
                [6, 10, 15, 22], [6, 12, 16, 18], [6, 13, 14, 24], [7, 8, 11, 17],
                [7, 9, 15, 23], [7, 10, 16, 21], [7, 20, 22, 24], [8, 10, 18, 19],
                [9, 12, 21, 22], [11, 12, 23, 24], [11, 13, 18, 22], [13, 15, 17, 19],
                [14, 15, 20, 21], [14, 16, 17, 23]]
    if v == 28:
        return [[0, 1, 23, 26], [0, 2, 10, 11], [0, 3, 16, 18], [0, 4, 15, 20],
                [0, 5, 8, 9], [0, 6, 22, 25], [0, 7, 14, 21], [0, 12, 17, 27],
                [0, 13, 19, 24], [1, 2, 24, 27], [1, 3, 11, 12], [1, 4, 17, 19],
                [1, 5, 14, 16], [1, 6, 9, 10], [1, 7, 20, 25], [1, 8, 15, 22],
                [1, 13, 18, 21], [2, 3, 21, 25], [2, 4, 12, 13], [2, 5, 18, 20],
                [2, 6, 15, 17], [2, 7, 19, 22], [2, 8, 14, 26], [2, 9, 16, 23],
                [3, 4, 22, 26], [3, 5, 7, 13], [3, 6, 14, 19], [3, 8, 20, 23],
                [3, 9, 15, 27], [3, 10, 17, 24], [4, 5, 23, 27], [4, 6, 7, 8],
                [4, 9, 14, 24], [4, 10, 16, 21], [4, 11, 18, 25], [5, 6, 21, 24],
                [5, 10, 15, 25], [5, 11, 17, 22], [5, 12, 19, 26], [6, 11, 16, 26],
                [6, 12, 18, 23], [6, 13, 20, 27], [7, 9, 17, 18], [7, 10, 26, 27],
                [7, 11, 23, 24], [7, 12, 15, 16], [8, 10, 18, 19], [8, 11, 21, 27],
                [8, 12, 24, 25], [8, 13, 16, 17], [9, 11, 19, 20], [9, 12, 21, 22],
                [9, 13, 25, 26], [10, 12, 14, 20], [10, 13, 22, 23], [11, 13, 14, 15],
                [14, 17, 23, 25], [14, 18, 22, 27], [15, 18, 24, 26], [15, 19, 21, 23],
                [16, 19, 25, 27], [16, 20, 22, 24], [17, 20, 21, 26]]
    if v == 37:
        return [[0, 1, 3, 24], [0, 2, 23, 36], [0, 4, 26, 32], [0, 5, 9, 31],
                [0, 6, 11, 15], [0, 7, 17, 25], [0, 8, 20, 27], [0, 10, 18, 30],
                [0, 12, 19, 29], [0, 13, 14, 16], [0, 21, 34, 35], [0, 22, 28, 33],
                [1, 2, 4, 25], [1, 5, 27, 33], [1, 6, 10, 32], [1, 7, 12, 16],
                [1, 8, 18, 26], [1, 9, 21, 28], [1, 11, 19, 31], [1, 13, 20, 30],
                [1, 14, 15, 17], [1, 22, 35, 36], [1, 23, 29, 34], [2, 3, 5, 26],
                [2, 6, 28, 34], [2, 7, 11, 33], [2, 8, 13, 17], [2, 9, 19, 27],
                [2, 10, 22, 29], [2, 12, 20, 32], [2, 14, 21, 31], [2, 15, 16, 18],
                [2, 24, 30, 35], [3, 4, 6, 27], [3, 7, 29, 35], [3, 8, 12, 34],
                [3, 9, 14, 18], [3, 10, 20, 28], [3, 11, 23, 30], [3, 13, 21, 33],
                [3, 15, 22, 32], [3, 16, 17, 19], [3, 25, 31, 36], [4, 5, 7, 28],
                [4, 8, 30, 36], [4, 9, 13, 35], [4, 10, 15, 19], [4, 11, 21, 29],
                [4, 12, 24, 31], [4, 14, 22, 34], [4, 16, 23, 33], [4, 17, 18, 20],
                [5, 6, 8, 29], [5, 10, 14, 36], [5, 11, 16, 20], [5, 12, 22, 30],
                [5, 13, 25, 32], [5, 15, 23, 35], [5, 17, 24, 34], [5, 18, 19, 21],
                [6, 7, 9, 30], [6, 12, 17, 21], [6, 13, 23, 31], [6, 14, 26, 33],
                [6, 16, 24, 36], [6, 18, 25, 35], [6, 19, 20, 22], [7, 8, 10, 31],
                [7, 13, 18, 22], [7, 14, 24, 32], [7, 15, 27, 34], [7, 19, 26, 36],
                [7, 20, 21, 23], [8, 9, 11, 32], [8, 14, 19, 23], [8, 15, 25, 33],
                [8, 16, 28, 35], [8, 21, 22, 24], [9, 10, 12, 33], [9, 15, 20, 24],
                [9, 16, 26, 34], [9, 17, 29, 36], [9, 22, 23, 25], [10, 11, 13, 34],
                [10, 16, 21, 25], [10, 17, 27, 35], [10, 23, 24, 26], [11, 12, 14, 35],
                [11, 17, 22, 26], [11, 18, 28, 36], [11, 24, 25, 27], [12, 13, 15, 36],
                [12, 18, 23, 27], [12, 25, 26, 28], [13, 19, 24, 28], [13, 26, 27, 29],
                [14, 20, 25, 29], [14, 27, 28, 30], [15, 21, 26, 30], [15, 28, 29, 31],
                [16, 22, 27, 31], [16, 29, 30, 32], [17, 23, 28, 32], [17, 30, 31, 33],
                [18, 24, 29, 33], [18, 31, 32, 34], [19, 25, 30, 34], [19, 32, 33, 35],
                [20, 26, 31, 35], [20, 33, 34, 36], [21, 27, 32, 36]]

    # Step 2 : this is function PBD_4_5_8_9_12
    PBD = PBD_4_5_8_9_12((v-1)/(k-1),check=False)

    # Step 3 : Theorem 7.20
    bibd = BIBD_from_PBD(PBD,v,k,check=False)

    if check:
        _check_pbd(bibd,v,[k])

    return bibd

def BIBD_from_PBD(PBD,v,k,check=True,base_cases={}):
    r"""
    Returns a `(v,k,1)`-BIBD from a `(r,K)`-PBD where `r=(v-1)/(k-1)`.

    This is Theorem 7.20 from [Stinson2004]_.

    INPUT:

    - ``v,k`` -- integers.

    - ``PBD`` -- A PBD on `r=(v-1)/(k-1)` points, such that for any block of
      ``PBD`` of size `s` there must exist a `((k-1)s+1,k,1)`-BIBD.

    - ``check`` (boolean) -- whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to ``True``
      by default.

    - ``base_cases`` -- caching system, for internal use.

    EXAMPLES::

        sage: from sage.combinat.designs.bibd import PBD_4_5_8_9_12
        sage: from sage.combinat.designs.bibd import BIBD_from_PBD
        sage: from sage.combinat.designs.bibd import _check_pbd
        sage: PBD = PBD_4_5_8_9_12(17)
        sage: bibd = _check_pbd(BIBD_from_PBD(PBD,52,4),52,[4])
    """
    r = (v-1)/(k-1)
    bibd = []
    for X in PBD:
        n = len(X)
        N = (k-1)*n+1
        if not (n,k) in base_cases:
            base_cases[n,k] = _relabel_bibd(BalancedIncompleteBlockDesign(N,k).blcks,N)

        for XX in base_cases[n,k]:
            if N-1 in XX:
                continue
            bibd.append([X[x//(k-1)] + (x%(k-1))*r for x in XX])

    for x in range(r):
        bibd.append([x+i*r for i in range(k-1)]+[v-1])

    if check:
        _check_pbd(bibd,v,[k])

    return bibd

def _check_pbd(B,v,S):
    r"""
    Checks that ``B`` is a PBD on `v` points with given block sizes.

    INPUT:

    - ``bibd`` -- a list of blocks

    - ``v`` (integer) -- number of points

    - ``S`` -- list of integers

    EXAMPLE::

        sage: designs.BalancedIncompleteBlockDesign(40,4).blocks() # indirect doctest
        [[0, 1, 2, 12], [0, 3, 6, 9], [0, 4, 8, 11], [0, 5, 7, 10],
         [0, 13, 26, 39], [0, 14, 28, 38], [0, 15, 25, 27],
         [0, 16, 32, 35], [0, 17, 34, 37], [0, 18, 33, 36],
        ...
    """
    from itertools import combinations
    from sage.graphs.graph import Graph

    if not all(len(X) in S for X in B):
        raise RuntimeError("This is not a nice honest PBD from the good old days !")

    g = Graph()
    m = 0
    for X in B:
        g.add_edges(list(combinations(X,2)))
        if g.size() != m+binomial(len(X),2):
            raise RuntimeError("This is not a nice honest PBD from the good old days !")
        m = g.size()

    if not (g.is_clique() and g.vertices() == range(v)):
        raise RuntimeError("This is not a nice honest PBD from the good old days !")

    return B

def _relabel_bibd(B,n):
    r"""
    Relabels the BIBD on `n` points and blocks of size k such that
    `\{0,...,k-2,n-1\},\{k-1,...,2k-3,n-1\},...,\{n-k,...,n-2,n-1\}` are blocks
    of the BIBD.

    INPUT:

    - ``B`` -- a list of blocks.

    - ``n`` (integer) -- number of points.

    EXAMPLE::

        sage: designs.BalancedIncompleteBlockDesign(40,4).blocks() # indirect doctest
        [[0, 1, 2, 12], [0, 3, 6, 9], [0, 4, 8, 11], [0, 5, 7, 10],
         [0, 13, 26, 39], [0, 14, 28, 38], [0, 15, 25, 27],
         [0, 16, 32, 35], [0, 17, 34, 37], [0, 18, 33, 36],
        ...
    """
    found = 0
    last = n-1
    d = {}
    for X in B:
        if last in X:
            for x in X:
                if x == last:
                    continue
                d[x] = found
                found += 1
            if found == n-1:
                break
    d[n-1] = n-1
    return [[d[x] for x in X] for X in B]

def PBD_4_5_8_9_12(v, check=True):
    """
    Returns a `(v,\{4,5,8,9,12\})-PBD` on `v` elements.

    A `(v,\{4,5,8,9,12\})`-PBD exists if and only if `v\equiv 0,1 \pmod 4`. The
    construction implemented here appears page 168 in [Stinson2004]_.

    INPUT:

    - ``v`` (integer)

    - ``check`` (boolean) -- whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to ``True``
      by default.

    EXAMPLES::

        sage: designs.BalancedIncompleteBlockDesign(40,4).blocks() # indirect doctest
        [[0, 1, 2, 12], [0, 3, 6, 9], [0, 4, 8, 11], [0, 5, 7, 10],
         [0, 13, 26, 39], [0, 14, 28, 38], [0, 15, 25, 27],
         [0, 16, 32, 35], [0, 17, 34, 37], [0, 18, 33, 36],
        ...
    """
    if not v%4 in [0,1]:
        raise ValueError
    if v == 0:
        return []
    if v == 13:
        PBD = v_4_1_BIBD(v, check=False)
    elif v == 28:
        PBD = v_4_1_BIBD(v, check=False)
    elif v == 29:
        TD47 = transversal_design(4,7)
        four_more_sets = [[28]+[i*7+j for j in range(7)] for i in range(4)]
        PBD = TD47 + four_more_sets
    elif v == 41:
        TD59 = transversal_design(5,9)
        PBD = ([[x for x in X if x<41] for X in TD59]
                +[[i*9+j for j in range(9)] for i in range(4)]
                +[[36,37,38,39,40]])
    elif v == 44:
        TD59 = transversal_design(5,9)
        PBD = ([[x for x in X if x<44] for X in TD59]
                +[[i*9+j for j in range(9)] for i in range(4)]
                +[[36,37,38,39,40,41,42,43]])
    elif v == 45:
        TD59 = transversal_design(5,9)
        PBD = (TD59+[[i*9+j for j in range(9)] for i in range(5)])
    elif v == 48:
        TD4_12 = transversal_design(4,12)
        PBD = (TD4_12+[[i*12+j for j in range(12)] for i in range(4)])
    elif v == 49:
        # Lemma 7.16 : A (49,{4,13})-PBD
        TD4_12 = transversal_design(4,12)

        # Replacing the block of size 13 with a BIBD
        BIBD_13_4 = v_4_1_BIBD(13)
        for i in range(4):
            for B in BIBD_13_4:
                TD4_12.append([i*12+x if x != 12 else 48
                               for x in B])

        PBD = TD4_12
    else:
        t,u = _get_t_u(v)
        TD = transversal_design(5,t)
        TD = [[x for x in X if x<4*t+u] for X in TD]
        for B in [range(t*i,t*(i+1)) for i in range(4)]:
            TD.extend(_PBD_4_5_8_9_12_closure([B]))

        if u > 1:
            TD.extend(_PBD_4_5_8_9_12_closure([range(4*t,4*t+u)]))

        PBD = TD

    if check:
        _check_pbd(PBD,v,[4,5,8,9,12])

    return PBD

def _PBD_4_5_8_9_12_closure(B):
    r"""
    Makes sure all blocks of `B` have size in `\{4,5,8,9,12\}`.

    This is a helper function for :meth:`PBD_4_5_8_9_12`. Given that
    `\{4,5,8,9,12\}` is PBD-closed, any block of size not in `\{4,5,8,9,12\}`
    can be decomposed further.

    EXAMPLES::

        sage: designs.BalancedIncompleteBlockDesign(40,4).blocks() # indirect doctest
        [[0, 1, 2, 12], [0, 3, 6, 9], [0, 4, 8, 11], [0, 5, 7, 10],
         [0, 13, 26, 39], [0, 14, 28, 38], [0, 15, 25, 27],
         [0, 16, 32, 35], [0, 17, 34, 37], [0, 18, 33, 36],
        ...
    """
    BB = []
    for X in B:
        if len(X) not in [4,5,8,9,12]:
            PBD = PBD_4_5_8_9_12(len(X), check = False)
            X = [[X[i] for i in XX] for XX in PBD]
            BB.extend(X)
        else:
            BB.append(X)
    return BB

table_7_1 = {
    0:{'t':-4,'u':16,'s':2},
    1:{'t':-4,'u':17,'s':2},
    4:{'t':1,'u':0,'s':1},
    5:{'t':1,'u':1,'s':1},
    8:{'t':1,'u':4,'s':1},
    9:{'t':1,'u':5,'s':1},
    12:{'t':1,'u':8,'s':1},
    13:{'t':1,'u':9,'s':1},
    16:{'t':4,'u':0,'s':0},
    17:{'t':4,'u':1,'s':0},
    20:{'t':5,'u':0,'s':0},
    21:{'t':5,'u':1,'s':0},
    24:{'t':5,'u':4,'s':0},
    25:{'t':5,'u':5,'s':0},
    28:{'t':5,'u':8,'s':1},
    29:{'t':5,'u':9,'s':1},
    32:{'t':8,'u':0,'s':0},
    33:{'t':8,'u':1,'s':0},
    36:{'t':8,'u':4,'s':0},
    37:{'t':8,'u':5,'s':0},
    40:{'t':8,'u':8,'s':0},
    41:{'t':8,'u':9,'s':1},
    44:{'t':8,'u':12,'s':1},
    45:{'t':8,'u':13,'s':1},
    }


def _get_t_u(v):
    r"""
    Returns the parameters of table 7.1 from [Stinson2004]_.

    INPUT:

    - ``v`` (integer)

    EXAMPLE::

        sage: from sage.combinat.designs.bibd import _get_t_u
        sage: _get_t_u(20)
        (5, 0)
    """
    # Table 7.1
    v = int(v)
    global table_7_1
    d = table_7_1[v%48]
    s = v//48
    if s < d['s']:
        raise RuntimeError("This should not have happened.")
    t = 12*s+d['t']
    u = d['u']
    return t,u
