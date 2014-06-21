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
    [[0, 1, 2, 12], [0, 3, 6, 9], [0, 4, 8, 10], [0, 5, 7, 11], [1, 3, 8, 11],
     [1, 4, 7, 9], [1, 5, 6, 10], [2, 3, 7, 10], [2, 4, 6, 11], [2, 5, 8, 9],
     [3, 4, 5, 12], [6, 7, 8, 12], [9, 10, 11, 12]]

`K_4`-decompositions of `K_v`
-----------------------------

Decompositions of `K_v` into `K_4` (i.e. `(v,4,1)`-BIBD) are built following
Douglas Stinson's construction as presented in [Stinson2004]_ page 167. It is
based upon the construction of `(v\{4,5,8,9,12\})`-PBD (see the doc of
:meth:`PBD_4_5_8_9_12`), knowing that a `(v\{4,5,8,9,12\})`-PBD on `v` points
can always be transformed into a `((k-1)v+1,4,1)`-BIBD, which covers all
possible cases of `(v,4,1)`-BIBD.

`K_5`-decompositions of `K_v`
-----------------------------

Decompositions of `K_v` into `K_4` (i.e. `(v,4,1)`-BIBD) are built following
Clayton Smith's construction [ClaytonSmith]_.

.. [ClaytonSmith] On the existence of `(v,5,1)`-BIBD.
  http://www.argilo.net/files/bibd.pdf
  Clayton Smith


Functions
---------
"""

from sage.categories.sets_cat import EmptySetError
from sage.misc.unknown import Unknown
from design_catalog import transversal_design
from block_design import BlockDesign
from sage.rings.arith import binomial
from sage.rings.arith import is_prime_power

def BalancedIncompleteBlockDesign(v,k,existence=False,use_LJCR=False):
    r"""
    Returns a BIBD of parameters `v,k`.

    A Balanced Incomplete Block Design of parameters `v,k` is a collection
    `\mathcal C` of `k`-subsets of `V=\{0,\dots,v-1\}` such that for any two
    distinct elements `x,y\in V` there is a unique element `S\in \mathcal C`
    such that `x,y\in S`.

    More general definitions sometimes involve a `\lambda` parameter, and we
    assume here that `\lambda=1`.

    For more information on BIBD, see the
    :wikipedia:`corresponding Wikipedia entry <Block_design#Definition_of_a_BIBD_.28or_2-design.29>`.

    INPUT:

    - ``v,k`` (integers)

    - ``existence`` (boolean) -- instead of building the design, returns:

        - ``True`` -- meaning that Sage knows how to build the design

        - ``Unknown`` -- meaning that Sage does not know how to build the
          design, but that the design may exist (see :mod:`sage.misc.unknown`).

        - ``False`` -- meaning that the design does not exist.

    - ``use_LJCR`` (boolean) -- whether to query the La Jolla Covering
      Repository for the design when Sage does not know how to build it (see
      :meth:`~sage.combinat.designs.covering_design.best_known_covering_design_www`). This
      requires internet.

    .. SEEALSO::

        * :func:`steiner_triple_system`
        * :func:`v_4_1_BIBD`
        * :func:`v_5_1_BIBD`

    TODO:

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
        sage: designs.BalancedIncompleteBlockDesign(16,6)
        Traceback (most recent call last):
        ...
        NotImplementedError: I don't know how to build a (16,6,1)-BIBD!

    TESTS::

        sage: designs.BalancedIncompleteBlockDesign(85,5,existence=True)
        True
        sage: _ = designs.BalancedIncompleteBlockDesign(85,5)

    A BIBD from a Finite Projective Plane::

        sage: _ = designs.BalancedIncompleteBlockDesign(21,5)

    Some trivial BIBD::

        sage: designs.BalancedIncompleteBlockDesign(10,10)
        Incidence structure with 10 points and 1 blocks
        sage: designs.BalancedIncompleteBlockDesign(1,10)
        Incidence structure with 1 points and 0 blocks

    Existence of BIBD with `k=3,4,5`::

        sage: [v for v in xrange(50) if designs.BalancedIncompleteBlockDesign(v,3,existence=True)]
        [1, 3, 7, 9, 13, 15, 19, 21, 25, 27, 31, 33, 37, 39, 43, 45, 49]
        sage: [v for v in xrange(100) if designs.BalancedIncompleteBlockDesign(v,4,existence=True)]
        [1, 4, 13, 16, 25, 28, 37, 40, 49, 52, 61, 64, 73, 76, 85, 88, 97]
        sage: [v for v in xrange(150) if designs.BalancedIncompleteBlockDesign(v,5,existence=True)]
        [1, 5, 21, 25, 41, 45, 61, 65, 81, 85, 101, 105, 121, 125, 141, 145]

    For `k > 5` there are currently very few constructions::

        sage: [v for v in xrange(150) if designs.BalancedIncompleteBlockDesign(v,6,existence=True) is True]
        [1, 6, 31, 91]
        sage: [v for v in xrange(150) if designs.BalancedIncompleteBlockDesign(v,6,existence=True) is Unknown]
        [16, 21, 36, 46, 51, 61, 66, 76, 81, 96, 106, 111, 121, 126, 136, 141]
    """
    if v == 1:
        if existence:
            return True
        return BlockDesign(v, [], test=False)

    if k == v:
        if existence:
            return True
        return BlockDesign(v, [range(v)], test=False)

    if v < k or k < 2 or (v-1) % (k-1) != 0 or (v*(v-1)) % (k*(k-1)) != 0:
        if existence:
            return False
        raise EmptySetError("No such design exists !")

    if k == 2:
        if existence:
            return True
        from itertools import combinations
        return BlockDesign(v, combinations(range(v),2), test = False)
    if k == 3:
        if existence:
            return v%6 == 1 or v%6 == 3
        return steiner_triple_system(v)
    if k == 4:
        if existence:
            return v%12 == 1 or v%12 == 4
        return BlockDesign(v, v_4_1_BIBD(v), test = False)
    if k == 5:
        if existence:
            return v%20 == 1 or v%20 == 5
        return BlockDesign(v, v_5_1_BIBD(v), test = False)

    from difference_family import difference_family

    if BIBD_from_TD(v,k,existence=True):
        if existence:
            return True
        return BlockDesign(v, BIBD_from_TD(v,k))
    if v == (k-1)**2+k and is_prime_power(k-1):
        if existence:
            return True
        from block_design import projective_plane
        return projective_plane(k-1)
    if difference_family(v,k,existence=True):
        if existence:
            return True
        G,D = difference_family(v,k)
        return BlockDesign(v, BIBD_from_difference_family(G,D,check=False), test=False)
    if use_LJCR:
        from covering_design import best_known_covering_design_www
        B = best_known_covering_design_www(v,k,2)

        # Is it a BIBD or just a good covering ?
        expected_n_of_blocks = binomial(v,2)/binomial(k,2)
        if B.low_bd() > expected_n_of_blocks:
            if existence:
                return False
            raise EmptySetError("No such design exists !")
        B = B.incidence_structure()
        if len(B.blcks) == expected_n_of_blocks:
            if existence:
                return True
            else:
                return B

    if existence:
        return Unknown
    else:
        raise NotImplementedError("I don't know how to build a ({},{},1)-BIBD!".format(v,k))

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
        EmptySetError: Steiner triple systems only exist for n = 1 mod 6 or n = 3 mod 6

    REFERENCE:

    .. [AndHonk97] A short course in Combinatorial Designs,
      Ian Anderson, Iiro Honkala,
      Internet Editions, Spring 1997,
      http://www.utu.fi/~honkala/designs.ps
    """

    name = "Steiner Triple System on "+str(n)+" elements"

    if n%6 == 3:
        t = (n-3) // 6
        Z = range(2*t+1)

        T = lambda (x,y) : x + (2*t+1)*y

        sts = [[(i,0),(i,1),(i,2)] for i in Z] + \
            [[(i,k),(j,k),(((t+1)*(i+j)) % (2*t+1),(k+1)%3)] for k in range(3) for i in Z for j in Z if i != j]

    elif n%6 == 1:

        t = (n-1) // 6
        N = range(2*t)
        T = lambda (x,y) : x+y*t*2 if (x,y) != (-1,-1) else n-1

        L1 = lambda i,j : (i+j) % ((n-1)//3)
        L = lambda i,j : L1(i,j)//2 if L1(i,j)%2 == 0 else t+(L1(i,j)-1)//2

        sts = [[(i,0),(i,1),(i,2)] for i in range(t)] + \
            [[(-1,-1),(i,k),(i-t,(k+1) % 3)] for i in range(t,2*t) for k in [0,1,2]] + \
            [[(i,k),(j,k),(L(i,j),(k+1) % 3)] for k in [0,1,2] for i in N for j in N if i < j]

    else:
        raise EmptySetError("Steiner triple systems only exist for n = 1 mod 6 or n = 3 mod 6")

    from sage.sets.set import Set
    sts = Set(map(lambda x: Set(map(T,x)),sts))

    return BlockDesign(n, sts, name=name)

def BIBD_from_TD(v,k,existence=False):
    r"""
    Returns a BIBD through TD-based constructions.

    INPUT:

    - ``v,k`` (integers) -- computes a `(v,k,1)`-BIBD.

    - ``existence`` (boolean) -- instead of building the design, returns:

        - ``True`` -- meaning that Sage knows how to build the design

        - ``Unknown`` -- meaning that Sage does not know how to build the
          design, but that the design may exist (see :mod:`sage.misc.unknown`).

        - ``False`` -- meaning that the design does not exist.

    This method implements three constructions:

    - If there exists a `TD(k,v)` and a `(v,k,1)`-BIBD then there exists a
      `(kv,k,1)`-BIBD.

      The BIBD is obtained from all blocks of the `TD`, and from the blocks of
      the `(v,k,1)`-BIBDs defined over the `k` groups of the `TD`.

    - If there exists a `TD(k,v)` and a `(v+1,k,1)`-BIBD then there exists a
      `(kv+1,k,1)`-BIBD.

      The BIBD is obtained from all blocks of the `TD`, and from the blocks of
      the `(v+1,k,1)`-BIBDs defined over the sets `V_1\cup \infty,\dots,V_k\cup
      \infty` where the `V_1,\dots,V_k` are the groups of the TD.

    - If there exists a `TD(k,v)` and a `(v+k,k,1)`-BIBD then there exists a
      `(kv+k,k,1)`-BIBD.

      The BIBD is obtained from all blocks of the `TD`, and from the blocks of
      the `(v+k,k,1)`-BIBDs defined over the sets `V_1\cup
      \{\infty_1,\dots,\infty_k\},\dots,V_k\cup \{\infty_1,\dots,\infty_k\}`
      where the `V_1,\dots,V_k` are the groups of the TD. By making sure that
      all copies of the `(v+k,k,1)`-BIBD contain the block
      `\{\infty_1,\dots,\infty_k\}`, the result is also a BIBD.

    These constructions can be found in
    `<http://www.argilo.net/files/bibd.pdf>`_.

    EXAMPLES:

    First construction::

        sage: from sage.combinat.designs.bibd import BIBD_from_TD
        sage: BIBD_from_TD(25,5,existence=True)
        True
        sage: _ = BlockDesign(25,BIBD_from_TD(25,5))

    Second construction::

        sage: from sage.combinat.designs.bibd import BIBD_from_TD
        sage: BIBD_from_TD(21,5,existence=True)
        True
        sage: _ = BlockDesign(21,BIBD_from_TD(21,5))

    Third construction::

        sage: from sage.combinat.designs.bibd import BIBD_from_TD
        sage: BIBD_from_TD(85,5,existence=True)
        True
        sage: _ = BlockDesign(85,BIBD_from_TD(85,5))

    No idea::

        sage: from sage.combinat.designs.bibd import BIBD_from_TD
        sage: BIBD_from_TD(20,5,existence=True)
        Unknown
        sage: BIBD_from_TD(20,5)
        Traceback (most recent call last):
        ...
        NotImplementedError: I do not know how to build a (20,5,1)-BIBD!
    """
    from orthogonal_arrays import transversal_design

    # First construction
    if (v%k == 0 and
        BalancedIncompleteBlockDesign(v//k,k,existence=True) and
        transversal_design(k,v//k,existence=True)):

        if existence:
            return True

        v = v//k
        BIBDvk = BalancedIncompleteBlockDesign(v,k)
        TDkv = transversal_design(k,v,check=False)

        BIBD = TDkv
        for i in range(k):
            BIBD.extend([[x+i*v for x in B] for B in BIBDvk])

    # Second construction
    elif ((v-1)%k == 0 and
        BalancedIncompleteBlockDesign((v-1)//k+1,k,existence=True) and
        transversal_design(k,(v-1)//k,existence=True)):

        if existence:
            return True

        v = (v-1)//k
        BIBDv1k = BalancedIncompleteBlockDesign(v+1,k)
        TDkv = transversal_design(k,v,check=False)

        inf = v*k
        BIBD = TDkv
        for i in range(k):
            BIBD.extend([[inf if x == v else x+i*v for x in B] for B in BIBDv1k])

    # Third construction
    elif ((v-k)%k == 0 and
        BalancedIncompleteBlockDesign((v-k)//k+k,k,existence=True) and
        transversal_design(k,(v-k)//k,existence=True)):

        if existence:
            return True

        v = (v-k)//k
        BIBDvpkk = BalancedIncompleteBlockDesign(v+k,k)
        TDkv = transversal_design(k,v,check=False)
        inf = v*k
        BIBD = TDkv

        # makes sure that [v,...,v+k-1] is a block of BIBDvpkk. Then, we remove it.
        BIBDvpkk = _relabel_bibd(BIBDvpkk,v+k)
        BIBDvpkk = [B for B in BIBDvpkk if min(B) < v]

        for i in range(k):
            BIBD.extend([[(x-v)+inf if x >= v else x+i*v for x in B] for B in BIBDvpkk])

        BIBD.append(range(k*v,v*k+k))

    # No idea ...
    else:
        if existence:
            return Unknown
        else:
            raise NotImplementedError("I do not know how to build a ({},{},1)-BIBD!".format(v,k))

    return BIBD

def BIBD_from_difference_family(G, D, check=True):
    r"""
    Return the BIBD associated to the difference family ``D`` on the group ``G``.

    Let `G` be a finite Abelian group. A *simple `(G,k)`-difference family* (or
    a *`(G,k,1)`-difference family*) is a family `B = \{B_1,B_2,\ldots,B_b\}` of
    `k`-subsets of `G` such that for each element of `G \backslash \{0\}` there
    exists a unique `s \in \{1,\ldots,b\}` and a unique pair of distinct
    elements `x,y \in B_s` such that `x - y = g`.

    If `\{B_1, B_2, \ldots, B_b\}` is a simple `(G,k)`-difference family then
    its set of translates `\{B_i + g; i \in \{1,\ldots,b\}, g \in G\}` is a
    `(v,k,1)`-BIBD where `v` is the cardinality of `G`.

    INPUT::

    - ``G`` - a finite additive Abelian group

    - ``D`` - a difference family on ``G``.

    - ``check`` - whether or not we check the output (default: ``True``)

    EXAMPLES::

        sage: G = Zmod(21)
        sage: D = [[0,1,4,14,16]]
        sage: print sorted(G(x-y) for x in D[0] for y in D[0] if x != y)
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

        sage: from sage.combinat.designs.bibd import BIBD_from_difference_family
        sage: BIBD_from_difference_family(G, D)
        [[0, 1, 4, 14, 16],
         [1, 2, 5, 15, 17],
         [2, 3, 6, 16, 18],
         [3, 4, 7, 17, 19],
         [4, 5, 8, 18, 20],
         [5, 6, 9, 19, 0],
         [6, 7, 10, 20, 1],
         [7, 8, 11, 0, 2],
         [8, 9, 12, 1, 3],
         [9, 10, 13, 2, 4],
         [10, 11, 14, 3, 5],
         [11, 12, 15, 4, 6],
         [12, 13, 16, 5, 7],
         [13, 14, 17, 6, 8],
         [14, 15, 18, 7, 9],
         [15, 16, 19, 8, 10],
         [16, 17, 20, 9, 11],
         [17, 18, 0, 10, 12],
         [18, 19, 1, 11, 13],
         [19, 20, 2, 12, 14],
         [20, 0, 3, 13, 15]]
    """
    r = {e:i for i,e in enumerate(G)}
    bibd = [[r[G(x)+g] for x in d] for d in D for g in r]
    if check:
        assert _check_pbd(bibd, G.cardinality(), [len(D[0])])
    return bibd



################
# (v,4,1)-BIBD #
################

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

    TESTS:

    Check that the `(25,4)` and `(37,4)`-difference family are available::

        sage: assert designs.difference_family(25,4,existence=True)
        sage: _ = designs.difference_family(25,4)
        sage: assert designs.difference_family(37,4,existence=True)
        sage: _ = designs.difference_family(37,4)
    """
    from sage.rings.finite_rings.constructor import FiniteField
    k = 4
    if v == 0:
        return []
    if v <= 12 or v%12 not in [1,4]:
        raise EmptySetError("A K_4-decomposition of K_v exists iif v=2,4 mod 12, v>12 or v==0")

    # Step 1. Base cases.
    if v == 13:
        from block_design import projective_plane
        return projective_plane(3).blocks()
    if v == 16:
        from block_design import AffineGeometryDesign
        return AffineGeometryDesign(2,1,FiniteField(4,'x')).blocks()
    if v == 25 or v == 37:
        from difference_family import difference_family
        G,D = difference_family(v,4)
        return BIBD_from_difference_family(G,D,check=False)
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
    r = (v-1) // (k-1)
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

    - ``B`` -- a list of blocks

    - ``v`` (integer) -- number of points

    - ``S`` -- list of integers `\geq 2`.

    EXAMPLE::

        sage: designs.BalancedIncompleteBlockDesign(40,4).blocks() # indirect doctest
        [[0, 1, 2, 12], [0, 3, 6, 9], [0, 4, 8, 10],
         [0, 5, 7, 11], [0, 13, 26, 39], [0, 14, 25, 28],
         [0, 15, 27, 38], [0, 16, 22, 32], [0, 17, 23, 34],
        ...
        sage: from sage.combinat.designs.bibd import _check_pbd
        sage: _check_pbd([[1],[]],1,[1,0])
        Traceback (most recent call last):
        ...
        RuntimeError: All integers of S must be >=2
    """
    from itertools import combinations
    from sage.graphs.graph import Graph

    if not all(len(X) in S for X in B):
        raise RuntimeError("Some block has wrong size: this is not a nice honest PBD from the good old days !")

    if any(x < 2 for x in S):
        raise RuntimeError("All integers of S must be >=2")

    if v == 0 or v == 1:
        if B:
            raise RuntimeError("This is not a nice honest PBD from the good old days!")
        else:
            return

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

def _relabel_bibd(B,n,p=None):
    r"""
    Relabels the BIBD on `n` points and blocks of size k such that
    `\{0,...,k-2,n-1\},\{k-1,...,2k-3,n-1\},...,\{n-k,...,n-2,n-1\}` are blocks
    of the BIBD.

    INPUT:

    - ``B`` -- a list of blocks.

    - ``n`` (integer) -- number of points.

    - ``p`` (optional) -- the point that will be labeled with n-1.

    EXAMPLE::

        sage: designs.BalancedIncompleteBlockDesign(40,4).blocks() # indirect doctest
        [[0, 1, 2, 12], [0, 3, 6, 9], [0, 4, 8, 10],
         [0, 5, 7, 11], [0, 13, 26, 39], [0, 14, 25, 28],
         [0, 15, 27, 38], [0, 16, 22, 32], [0, 17, 23, 34],
        ...
    """
    if p is None:
        p = n-1
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
    d[p] = n-1
    return [[d[x] for x in X] for X in B]

def PBD_4_5_8_9_12(v, check=True):
    """
    Returns a `(v,\{4,5,8,9,12\})-`PBD on `v` elements.

    A `(v,\{4,5,8,9,12\})`-PBD exists if and only if `v\equiv 0,1 \pmod 4`. The
    construction implemented here appears page 168 in [Stinson2004]_.

    INPUT:

    - ``v`` -- an integer congruent to `0` or `1` modulo `4`.

    - ``check`` (boolean) -- whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to ``True``
      by default.

    EXAMPLES::

        sage: designs.BalancedIncompleteBlockDesign(40,4).blocks() # indirect doctest
        [[0, 1, 2, 12], [0, 3, 6, 9], [0, 4, 8, 10],
         [0, 5, 7, 11], [0, 13, 26, 39], [0, 14, 25, 28],
         [0, 15, 27, 38], [0, 16, 22, 32], [0, 17, 23, 34],
        ...

    Check that :trac:`16476` is fixed::

        sage: from sage.combinat.designs.bibd import PBD_4_5_8_9_12
        sage: for v in (0,1,4,5,8,9,12,13,16,17,20,21,24,25):
        ....:     _ = PBD_4_5_8_9_12(v)
    """
    if not v%4 in [0,1]:
        raise ValueError
    if v <= 1:
        PBD = []
    elif v <= 12:
        PBD = [range(v)]
    elif v == 13 or v == 28:
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
        [[0, 1, 2, 12], [0, 3, 6, 9], [0, 4, 8, 10],
         [0, 5, 7, 11], [0, 13, 26, 39], [0, 14, 25, 28],
         [0, 15, 27, 38], [0, 16, 22, 32], [0, 17, 23, 34],
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

################
# (v,5,1)-BIBD #
################


def v_5_1_BIBD(v, check=True):
    r"""
    Returns a `(v,5,1)`-BIBD.

    This method follows the constuction from [ClaytonSmith]_.

    INPUT:

    - ``v`` (integer)

    .. SEEALSO::

        * :meth:`BalancedIncompleteBlockDesign`

    EXAMPLES::

        sage: from sage.combinat.designs.bibd import v_5_1_BIBD
        sage: i = 0
        sage: while i<200:
        ....:    i += 20
        ....:    _ = v_5_1_BIBD(i+1)
        ....:    _ = v_5_1_BIBD(i+5)

    TESTS:

    Check that the needed difference families are there::

        sage: for v in [21,41,61,81,141,161,281]:
        ....:     assert designs.difference_family(v,5,existence=True)
        ....:     _ = designs.difference_family(v,5)
    """
    v = int(v)

    assert (v > 1)
    assert (v%20 == 5 or v%20 == 1)  # note: equivalent to (v-1)%4 == 0 and (v*(v-1))%20 == 0

    # Lemma 27
    if v%5 == 0 and (v//5)%4 == 1 and is_prime_power(v//5):
        bibd = BIBD_5q_5_for_q_prime_power(v//5)
    # Lemma 28
    elif v in [21,41,61,81,141,161,281]:
        from difference_family import difference_family
        G,D = difference_family(v,5)
        bibd = BIBD_from_difference_family(G, D, check=False)
    # Lemma 29
    elif v == 165:
        bibd = BIBD_from_PBD(v_5_1_BIBD(41,check=False),165,5,check=False)
    elif v == 181:
        bibd = BIBD_from_PBD(v_5_1_BIBD(45,check=False),181,5,check=False)
    elif v in (201,285,301,401,421,425):
        # Call directly the BIBD_from_TD function
        bibd = BIBD_from_TD(v,5)
    # Theorem 31.2
    elif (v-1)//4 in [80, 81, 85, 86, 90, 91, 95, 96, 110, 111, 115, 116, 120, 121, 250, 251, 255, 256, 260, 261, 265, 266, 270, 271]:
        r = (v-1)//4
        if r <= 96:
            k,t,u = 5, 16, r-80
        elif r <= 121:
            k,t,u = 10, 11, r-110
        else:
            k,t,u = 10, 25, r-250
        bibd = BIBD_from_PBD(PBD_from_TD(k,t,u),v,5,check=False)

    else:
        r,s,t,u = _get_r_s_t_u(v)
        bibd = BIBD_from_PBD(PBD_from_TD(5,t,u),v,5,check=False)

    if check:
        _check_pbd(bibd,v,[5])

    return bibd

def _get_r_s_t_u(v):
    r"""
    Implements the table from [ClaytonSmith]_

    Returns the parameters ``r,s,t,u`` associated with an integer ``v``.

    INPUT:

    - ``v`` (integer)

    EXAMPLES::

        sage: from sage.combinat.designs.bibd import _get_r_s_t_u
        sage: _get_r_s_t_u(25)
        (6, 0, 1, 1)
    """
    r = int((v-1)/4)
    s = r//150
    x = r%150

    if   x == 0:   t,u = 30*s-5,  25
    elif x == 1:   t,u = 30*s-5,  26
    elif x <= 21:  t,u = 30*s+1,  x-5
    elif x == 25:  t,u = 30*s+5,  0
    elif x == 26:  t,u = 30*s+5,  1
    elif x == 30:  t,u = 30*s+5,  5
    elif x <= 51:  t,u = 30*s+5,  x-25
    elif x <= 66:  t,u = 30*s+11, x-55
    elif x <= 96:  t,u = 30*s+11, x-55
    elif x <= 121: t,u = 30*s+11, x-55
    elif x <= 146: t,u = 30*s+25, x-125

    return r,s,t,u

def PBD_from_TD(k,t,u):
    r"""
    Returns a `(kt,\{k,t\})`-PBD if `u=0` and a `(kt+u,\{k,k+1,t,u\})`-PBD otherwise.

    This is theorem 23 from [ClaytonSmith]_. The PBD is obtained from the blocks
    a truncated `TD(k+1,t)`, to which are added the blocks corresponding to the
    groups of the TD. When `u=0`, a `TD(k,t)` is used instead.

    INPUT:

    - ``k,t,u`` -- integers such that `0\leq u \leq t`.

    EXAMPLES::

        sage: from sage.combinat.designs.bibd import PBD_from_TD
        sage: from sage.combinat.designs.bibd import _check_pbd
        sage: PBD = PBD_from_TD(2,2,1); PBD
        [[0, 2, 4], [0, 3], [1, 2], [1, 3, 4], [0, 1], [2, 3]]
        sage: _check_pbd(PBD,2*2+1,[2,3])
        [[0, 2, 4], [0, 3], [1, 2], [1, 3, 4], [0, 1], [2, 3]]

    """
    from orthogonal_arrays import transversal_design
    TD = transversal_design(k+bool(u),t, check=False)
    TD = [[x for x in X if x<k*t+u] for X in TD]
    for i in range(k):
        TD.append(range(t*i,t*i+t))
    if u>=2:
        TD.append(range(k*t,k*t+u))
    return TD

def BIBD_5q_5_for_q_prime_power(q):
    r"""
    Returns a `(5q,5,1)`-BIBD with `q\equiv 1\pmod 4` a prime power.

    See Theorem 24 [ClaytonSmith]_.

    INPUT:

    - ``q`` (integer) -- a prime power such that `q\equiv 1\pmod 4`.

    EXAMPLES::

        sage: from sage.combinat.designs.bibd import BIBD_5q_5_for_q_prime_power
        sage: for q in [25, 45, 65, 85, 125, 145, 185, 205, 305, 405, 605]: # long time
        ....:     _ = BIBD_5q_5_for_q_prime_power(q/5)                      # long time
    """
    from sage.rings.arith import is_prime_power
    from sage.rings.finite_rings.constructor import FiniteField

    if q%4 != 1 or not is_prime_power(q):
        raise ValueError("q is not a prime power or q%4!=1.")

    d = (q-1)/4
    B = []
    F = FiniteField(q,'x')
    a = F.primitive_element()
    L = {b:i for i,b in enumerate(F)}
    for b in L:
        B.append([i*q + L[b] for i in range(5)])
        for i in range(5):
            for j in range(d):
                B.append([        i*q + L[b          ],
                          ((i+1)%5)*q + L[ a**j+b    ],
                          ((i+1)%5)*q + L[-a**j+b    ],
                          ((i+4)%5)*q + L[ a**(j+d)+b],
                          ((i+4)%5)*q + L[-a**(j+d)+b],
                          ])

    return B
