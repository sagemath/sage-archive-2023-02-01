r"""
Difference families

This module gathers everything related to difference families. One can build a
difference family (or check that it can be built) with :func:`difference_family`::

    sage: G,F = designs.difference_family(13,4,1)

It defines the following functions:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`is_difference_family` | Check if the input is a (``k``, ``l``)-difference family.
    :func:`singer_difference_set` | Return a difference set associated to hyperplanes in a projective space.
    :func:`difference_family` | Return a (``k``, ``l``)-difference family on an Abelian group of size ``v``.

REFERENCES:

.. [Wi72] R. M. Wilson "Cyclotomy and difference families in elementary Abelian
   groups", J. of Num. Th., 4 (1972), pp. 17-47.

Functions
---------
"""
#*****************************************************************************
#       Copyright (C) 2014 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.sets_cat import EmptySetError
import sage.rings.arith as arith
from sage.misc.unknown import Unknown
from sage.rings.integer import Integer

def group_law(G):
    r"""
    Return a triple ``(identity, operation, inverse)`` that define the
    operations on the group ``G``.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import group_law
        sage: group_law(Zmod(3))
        (0, <built-in function add>, <built-in function neg>)
        sage: group_law(SymmetricGroup(5))
        ((), <built-in function mul>, <built-in function inv>)
        sage: group_law(VectorSpace(QQ,3))
        ((0, 0, 0), <built-in function add>, <built-in function neg>)
    """
    import operator
    from sage.categories.groups import Groups
    from sage.categories.additive_groups import AdditiveGroups

    if G in Groups():            # multiplicative groups
        return (G.one(), operator.mul, operator.inv)
    elif G in AdditiveGroups():  # additive groups
        return (G.zero(), operator.add, operator.neg)
    else:
        raise ValueError("%s does not seem to be a group"%G)

def block_stabilizer(G, B):
    r"""
    Compute the left stabilizer of the block ``B`` under the action of ``G``.

    This function return the list of all `x\in G` such that `x\cdot B=B` (as a
    set).

    INPUT:

    - ``G`` -- a group (additive or multiplicative).

    - ``B`` -- a subset of ``G``.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import block_stabilizer

        sage: Z8 = Zmod(8)
        sage: block_stabilizer(Z8, [Z8(0),Z8(2),Z8(4),Z8(6)])
        [0, 2, 4, 6]
        sage: block_stabilizer(Z8, [Z8(0),Z8(2)])
        [0]

        sage: C = cartesian_product([Zmod(4),Zmod(3)])
        sage: block_stabilizer(C, [C((0,0)),C((2,0)),C((0,1)),C((2,1))])
        [(0, 0), (2, 0)]

        sage: b = map(Zmod(45),[1, 3, 7, 10, 22, 25, 30, 35, 37, 38, 44])
        sage: block_stabilizer(Zmod(45),b)
        [0]
    """
    if not B:
        return list(G)
    identity, op, inv = group_law(G)
    b0 = inv(B[0])
    S = []
    for b in B:
        # fun: if we replace +(-b) with -b it completely fails!!
        bb0 = op(b,b0) # bb0 = b-B[0]
        if all(op(bb0,c) in B for c in B):
            S.append(bb0)
    return S

def is_difference_family(G, D, v=None, k=None, l=None, verbose=False):
    r"""
    Check wether ``D`` forms a difference family in the group ``G``.

    INPUT:

    - ``G`` - group of cardinality ``v``

    - ``D`` - a set of ``k``-subsets of ``G``

    - ``v``, ``k`` and ``l`` - optional parameters of the difference family

    - ``verbose`` - whether to print additional information

    .. SEEALSO::

        :func:`difference_family`

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import is_difference_family
        sage: G = Zmod(21)
        sage: D = [[0,1,4,14,16]]
        sage: is_difference_family(G, D, 21, 5)
        True

        sage: G = Zmod(41)
        sage: D = [[0,1,4,11,29],[0,2,8,17,21]]
        sage: is_difference_family(G, D, verbose=True)
        Too few:
          5 is obtained 0 times in blocks []
          14 is obtained 0 times in blocks []
          27 is obtained 0 times in blocks []
          36 is obtained 0 times in blocks []
        Too much:
          4 is obtained 2 times in blocks [0, 1]
          13 is obtained 2 times in blocks [0, 1]
          28 is obtained 2 times in blocks [0, 1]
          37 is obtained 2 times in blocks [0, 1]
        False
        sage: D = [[0,1,4,11,29],[0,2,8,17,22]]
        sage: is_difference_family(G, D)
        True

        sage: G = Zmod(61)
        sage: D = [[0,1,3,13,34],[0,4,9,23,45],[0,6,17,24,32]]
        sage: is_difference_family(G, D)
        True

        sage: G = AdditiveAbelianGroup([3]*4)
        sage: a,b,c,d = G.gens()
        sage: D = [[d, -a+d, -c+d, a-b-d, b+c+d],
        ....:      [c, a+b-d, -b+c, a-b+d, a+b+c],
        ....:      [-a-b+c+d, a-b-c-d, -a+c-d, b-c+d, a+b],
        ....:      [-b-d, a+b+d, a-b+c-d, a-b+c, -b+c+d]]
        sage: is_difference_family(G, D)
        True

    The following example has a third block with a non-trivial stabilizer::

        sage: G = Zmod(15)
        sage: D = [[0,1,4],[0,2,9],[0,5,10]]
        sage: is_difference_family(G,D,verbose=True)
        It is a (15,3,1)-difference family
        True

    The function also supports multiplicative groups (non necessarily Abelian)::

        sage: G = DihedralGroup(8)
        sage: x,y = G.gens()
        sage: i = G.one()
        sage: D1 = [[i,x,x^4], [i,x^2, y*x], [i,x^5,y], [i,x^6,y*x^2], [i,x^7,y*x^5]]
        sage: is_difference_family(G, D1, 16, 3, 2)
        True
        sage: from sage.combinat.designs.bibd import BIBD_from_difference_family
        sage: bibd = BIBD_from_difference_family(G,D1,lambd=2)
    """
    import operator

    identity, mul, inv = group_law(G)

    Glist = list(G)

    D = [map(G,d) for d in D]

    # Check v (and define it if needed)
    if v is None:
        v = len(Glist)
    else:
        if len(Glist) != v:
            if verbose:
                print "G must have cardinality v (=%d)"%int(v)
            return False

    # Check k (and define it if needed)
    if k is None:
        k = len(D[0])
    else:
        k = int(k)

    for d in D:
        if len(d) != k:
            if verbose:
                print "the block {} does not have length {}".format(d,k)
            return False

    # Check l (and define it if needed)
    #
    # - nb_diff: the number of pairs (with multiplicity) covered by the BIBD
    #            generated by the DF.
    #
    # - stab: the stabilizer of each set.
    nb_diff = 0
    stab = []
    for d in D:
        s = block_stabilizer(G,d)
        stab.append(s)
        nb_diff += k*(k-1) / len(s)
    if l is None:
        if nb_diff % (v-1) != 0:
            if verbose:
                print "the number of differences (={}) must be a multiple of v-1={}".format(nbdiff,v-1)
            return False
        l = nb_diff // (v-1)
    else:
        if nb_diff != l*(v-1):
            if verbose:
                print "the number of differences (={}) is not equal to l*(v-1) = {}".format(nb_diff, l*(v-1))
            return False

    # Check that every x \in G-{0},occurs exactly l times as a difference
    counter = {g: 0 for g in Glist}
    where   = {g: set() for g in Glist}
    del counter[identity]

    for i,d in enumerate(D):
        tmp_counter = {}
        for b in d:
            for c in d:
                if b == c:
                    continue
                gg = mul(b,inv(c)) # = b-c or bc^{-1}
                if gg not in tmp_counter:
                    tmp_counter[gg] = 0
                where[gg].add(i)
                tmp_counter[gg] += 1

        if sum(tmp_counter.itervalues()) != k*(k-1):
            if verbose:
                print "repeated element in the {}-th block {}".format(i,dd)
            return False

        # Normalized number of occurrences added to counter
        stabi = len(stab[i])
        for gg in tmp_counter:
            counter[gg] += tmp_counter[gg]//stabi

    # Check the counter and report any error
    too_few  = []
    too_much = []
    for g in Glist:
        if g == identity:
            continue
        if counter[g] < l:
            if verbose:
                too_few.append(g)
            else:
                return False
        if counter[g] > l:
            if verbose:
                too_much.append(g)
            else:
                return False

    if too_few:
        print "Too few:"
        for g in too_few:
            print "  {} is obtained {} times in blocks {}".format(
                        g,counter[g],sorted(where[g]))
    if too_much:
        print "Too much:"
        for g  in too_much:
            print "  {} is obtained {} times in blocks {}".format(
                        g,counter[g],sorted(where[g]))
    if too_few or too_much:
        return False

    if verbose:
        print "It is a ({},{},{})-difference family".format(v,k,l)
    return True

def singer_difference_set(q,d):
    r"""
    Return a difference set associated to the set of hyperplanes in a projective
    space of dimension `d` over `GF(q)`.

    A Singer difference set has parameters:

    .. MATH::

        v = \frac{q^{d+1}-1}{q-1}, \quad
        k = \frac{q^d-1}{q-1}, \quad
        \lambda = \frac{q^{d-1}-1}{q-1}.

    The idea of the construction is as follows. One consider the finite field
    `GF(q^{d+1})` as a vector space of dimension `d+1` over `GF(q)`. The set of
    `GF(q)`-lines in `GF(q^{d+1})` is a projective plane and its set of
    hyperplanes form a balanced incomplete block design.

    Now, considering a multiplicative generator `z` of `GF(q^{d+1})`, we get a
    transitive action of a cyclic group on our projective plane from which it is
    possible to build a difference set.

    The construction is given in details in [Stinson2004]_, section 3.3.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import singer_difference_set, is_difference_family
        sage: G,D = singer_difference_set(3,2)
        sage: is_difference_family(G,D,verbose=True)
        It is a (13,4,1)-difference family
        True

        sage: G,D = singer_difference_set(4,2)
        sage: is_difference_family(G,D,verbose=True)
        It is a (21,5,1)-difference family
        True

        sage: G,D = singer_difference_set(3,3)
        sage: is_difference_family(G,D,verbose=True)
        It is a (40,13,4)-difference family
        True

        sage: G,D = singer_difference_set(9,3)
        sage: is_difference_family(G,D,verbose=True)
        It is a (820,91,10)-difference family
        True
    """
    q = Integer(q)
    assert q.is_prime_power()
    assert d >= 2

    from sage.rings.finite_rings.constructor import GF
    from sage.rings.finite_rings.conway_polynomials import conway_polynomial
    from sage.rings.finite_rings.integer_mod_ring import Zmod

    # build a polynomial c over GF(q) such that GF(q)[x] / (c(x)) is a
    # GF(q**(d+1)) and such that x is a multiplicative generator.
    p,e = q.factor()[0]
    c = conway_polynomial(p,e*(d+1))
    if e != 1:  # i.e. q is not a prime, so we factorize c over GF(q) and pick
                # one of its factor
        K = GF(q,'z')
        c = c.change_ring(K).factor()[0][0]
    else:
        K = GF(q)
    z = c.parent().gen()

    # Now we consider the GF(q)-subspace V spanned by (1,z,z^2,...,z^(d-1)) inside
    # GF(q^(d+1)). The multiplication by z is an automorphism of the
    # GF(q)-projective space built from GF(q^(d+1)). The difference family is
    # obtained by taking the integers i such that z^i belong to V.
    powers = [0]
    i = 1
    x = z
    k = (q**d-1)//(q-1)
    while len(powers) < k:
        if x.degree() <= (d-1):
            powers.append(i)
        x = (x*z).mod(c)
        i += 1

    return Zmod((q**(d+1)-1)//(q-1)), [powers]

def difference_family(v, k, l=1, existence=False, check=True):
    r"""
    Return a (``k``, ``l``)-difference family on an Abelian group of cardinality ``v``.

    Let `G` be a finite Abelian group. For a given subset `D` of `G`, we define
    `\Delta D` to be the multi-set of differences `\Delta D = \{x - y; x \in D,
    y \in D, x \not= y\}`. A `(G,k,\lambda)`-*difference family* is a collection
    of `k`-subsets of `G`, `D = \{D_1, D_2, \ldots, D_b\}` such that the union
    of the difference sets `\Delta D_i` for `i=1,...b`, seen as a multi-set,
    contains each element of `G \backslash \{0\}` exactly `\lambda`-times.

    When there is only one block, i.e. `\lambda(v - 1) = k(k-1)`, then a
    `(G,k,\lambda)`-difference family is also called a *difference set*.

    See also :wikipedia:`Difference_set`.

    If there is no such difference family, an ``EmptySetError`` is raised and if
    there is no construction at the moment ``NotImplementedError`` is raised.

    EXAMPLES::

        sage: K,D = designs.difference_family(73,4)
        sage: D
        [[0, 1, 8, 64],
         [0, 25, 54, 67],
         [0, 41, 36, 69],
         [0, 3, 24, 46],
         [0, 2, 16, 55],
         [0, 50, 35, 61]]

        sage: K,D = designs.difference_family(337,7)
        sage: D
        [[1, 175, 295, 64, 79, 8, 52],
         [326, 97, 125, 307, 142, 249, 102],
         [121, 281, 310, 330, 123, 294, 226],
         [17, 279, 297, 77, 332, 136, 210],
         [150, 301, 103, 164, 55, 189, 49],
         [35, 59, 215, 218, 69, 280, 135],
         [289, 25, 331, 298, 252, 290, 200],
         [191, 62, 66, 92, 261, 180, 159]]

    For `k=6,7` we look at the set of small prime powers for which a
    construction is available::

        sage: def prime_power_mod(r,m):
        ....:     k = m+r
        ....:     while True:
        ....:         if is_prime_power(k):
        ....:             yield k
        ....:         k += m

        sage: from itertools import islice
        sage: l6 = {True:[], False: [], Unknown: []}
        sage: for q in islice(prime_power_mod(1,30), 60):
        ....:     l6[designs.difference_family(q,6,existence=True)].append(q)
        sage: l6[True]
        [31, 121, 151, 181, 211, ...,  3061, 3121, 3181]
        sage: l6[Unknown]
        [61]
        sage: l6[False]
        []

        sage: l7 = {True: [], False: [], Unknown: []}
        sage: for q in islice(prime_power_mod(1,42), 60):
        ....:     l7[designs.difference_family(q,7,existence=True)].append(q)
        sage: l7[True]
        [337, 421, 463, 883, 1723, 3067, 3319, 3529, 3823, 3907, 4621, 4957, 5167]
        sage: l7[Unknown]
        [43, 127, 169, 211, ..., 4999, 5041, 5209]
        sage: l7[False]
        []

    Other constructions for `\lambda > 1`::

        sage: for v in xrange(2,100):
        ....:     constructions = []
        ....:     for k in xrange(2,10):
        ....:         for l in xrange(2,10):
        ....:             if designs.difference_family(v,k,l,existence=True):
        ....:                 constructions.append((k,l))
        ....:                 _ = designs.difference_family(v,k,l)
        ....:     if constructions:
        ....:         print "%2d: %s"%(v, ', '.join('(%d,%d)'%(k,l) for k,l in constructions))
         2: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
         3: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
         4: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
         5: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
         7: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
         8: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
         9: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        11: (3,2), (4,3), (4,6), (5,2), (5,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        13: (3,2), (4,3), (5,4), (5,5), (6,5), (7,6), (8,7), (9,8)
        15: (4,6), (5,6), (7,3)
        16: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        17: (3,2), (4,3), (5,4), (5,5), (6,5), (7,6), (8,7), (9,8)
        19: (3,2), (4,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,4), (9,5), (9,6), (9,7), (9,8)
        21: (4,3), (6,3), (6,5)
        22: (4,2), (6,5), (7,4), (8,8)
        23: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        25: (3,2), (4,3), (5,4), (6,5), (7,6), (7,7), (8,7), (9,8)
        27: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        28: (3,2), (6,5)
        29: (3,2), (4,3), (5,4), (6,5), (7,3), (7,6), (8,4), (8,6), (8,7), (9,8)
        31: (3,2), (4,2), (4,3), (5,2), (5,4), (6,5), (7,6), (8,7), (9,8)
        32: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        33: (5,5), (6,5)
        34: (4,2)
        35: (5,2), (8,4)
        37: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,2), (9,3), (9,8)
        39: (6,5)
        40: (3,2)
        41: (3,2), (4,3), (5,4), (6,3), (6,5), (7,6), (8,7), (9,8)
        43: (3,2), (4,2), (4,3), (5,4), (6,5), (7,2), (7,3), (7,6), (8,4), (8,7), (9,8)
        46: (4,2), (6,2)
        47: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        49: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,3), (9,8)
        51: (5,2), (6,3)
        53: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        55: (9,4)
        57: (7,3)
        59: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        61: (3,2), (4,3), (5,4), (6,2), (6,3), (6,5), (7,6), (8,7), (9,8)
        64: (3,2), (4,3), (5,4), (6,5), (7,2), (7,6), (8,7), (9,8)
        67: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        71: (3,2), (4,3), (5,2), (5,4), (6,5), (7,3), (7,6), (8,4), (8,7), (9,8)
        73: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        75: (5,2)
        79: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        81: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        83: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        85: (7,2), (7,3), (8,2)
        89: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,8)
        97: (3,2), (4,3), (5,4), (6,5), (7,6), (8,7), (9,3), (9,8)

    TESTS:

    Check more of the Wilson constructions from [Wi72]_::

        sage: Q5 = [241, 281,421,601,641, 661, 701, 821,881]
        sage: Q9 = [73, 1153, 1873, 2017]
        sage: Q15 = [76231]
        sage: Q4 = [13, 73, 97, 109, 181, 229, 241, 277, 337, 409, 421, 457]
        sage: Q8 = [1009, 3137, 3697]
        sage: for Q,k in [(Q4,4),(Q5,5),(Q8,8),(Q9,9),(Q15,15)]:
        ....:     for q in Q:
        ....:         assert designs.difference_family(q,k,1,existence=True) is True
        ....:         _ = designs.difference_family(q,k,1)

    Check Singer difference sets::

        sage: sgp = lambda q,d: ((q**(d+1)-1)//(q-1), (q**d-1)//(q-1), (q**(d-1)-1)//(q-1))

        sage: for q in range(2,10):
        ....:     if is_prime_power(q):
        ....:         for d in [2,3,4]:
        ....:           v,k,l = sgp(q,d)
        ....:           assert designs.difference_family(v,k,l,existence=True) is True
        ....:           _ = designs.difference_family(v,k,l)

    Check twin primes difference sets::

        sage: for p in [3,5,7,9,11]:
        ....:     v = p*(p+2); k = (v-1)/2;  lmbda = (k-1)/2
        ....:     G,D = designs.difference_family(v,k,lmbda)

    Check the database:

        sage: from sage.combinat.designs.database import DF
        sage: for v,k,l in DF:
        ....:     df = designs.difference_family(v,k,l,check=True)

    .. TODO::

        Implement recursive constructions from Buratti "Recursive for difference
        matrices and relative difference families" (1998) and Jungnickel
        "Composition theorems for difference families and regular planes" (1978)
    """
    from block_design import are_hyperplanes_in_projective_geometry_parameters

    from database import DF

    if (v,k,l) in DF:
        if existence:
            return True

        vv, blocks = DF[v,k,l].iteritems().next()

        # Build the group
        from sage.rings.finite_rings.integer_mod_ring import Zmod
        if len(vv) == 1:
            G = Zmod(vv[0])
        else:
            from sage.categories.cartesian_product import cartesian_product
            G = cartesian_product([Zmod(i) for i in vv])

        df = [[G(i) for i in b] for b in blocks]

        if check:
            assert is_difference_family(G, df, v=v, k=k, l=l), "Sage built an invalid ({},{},{})-DF!".format(v,k,l)

        return G,df

    e = k*(k-1)
    t = l*(v-1) // e  # number of blocks

    D = None

    factorization = arith.factor(v)

    if len(factorization) == 1:  # i.e. is v a prime power
        from sage.rings.finite_rings.constructor import GF
        G = K = GF(v,'z')
        x = K.multiplicative_generator()

        if l == (k-1):
            if existence:
                return True
            return K, K.cyclotomic_cosets(x**((v-1)//k))[1:]

        if t == 1:
            # some of the difference set constructions VI.18.48 from the
            # Handbook of combinatorial designs
            # q = 3 mod 4
            if v%4 == 3 and k == (v-1)//2:
                if existence:
                    return True
                D = K.cyclotomic_cosets(x**2, [1])

            # q = 4t^2 + 1, t odd
            elif v%8 == 5 and k == (v-1)//4 and arith.is_square((v-1)//4):
                if existence:
                    return True
                D = K.cyclotomic_cosets(x**4, [1])

            # q = 4t^2 + 9, t odd
            elif v%8 == 5 and k == (v+3)//4 and arith.is_square((v-9)//4):
                if existence:
                    return True
                D = K.cyclotomic_cosets(x**4, [1])
                D[0].insert(0,K.zero())

        if D is None and l == 1:
            one = K.one()

            # Wilson (1972), Theorem 9
            if k%2 == 1:
                m = (k-1) // 2
                xx = x**m
                to_coset = {x**i * xx**j: i for i in xrange(m) for j in xrange((v-1)/m)}
                r = x ** ((v-1) // k)  # primitive k-th root of unity
                if len(set(to_coset[r**j-one] for j in xrange(1,m+1))) == m:
                    if existence:
                        return True
                    B = [r**j for j in xrange(k)]  # = H^((k-1)t) whose difference is
                                                   # H^(mt) (r^i - 1, i=1,..,m)
                    # Now pick representatives a translate of R for by a set of
                    # representatives of H^m / H^(mt)
                    D = [[x**(i*m) * b for b in B] for i in xrange(t)]

            # Wilson (1972), Theorem 10
            else:
                m = k//2
                xx = x**m
                to_coset = {x**i * xx**j: i for i in xrange(m) for j in xrange((v-1)/m)}
                r = x ** ((v-1) // (k-1))  # primitive (k-1)-th root of unity
                if (all(to_coset[r**j-one] != 0 for j in xrange(1,m)) and
                    len(set(to_coset[r**j-one] for j in xrange(1,m))) == m-1):
                    if existence:
                        return True
                    B = [K.zero()] + [r**j for j in xrange(k-1)]
                    D = [[x**(i*m) * b for b in B] for i in xrange(t)]

            # Wilson (1972), Theorem 11
            if D is None and k == 6:
                r = x**((v-1)//3)  # primitive cube root of unity
                r2 = r*r
                xx = x**5
                to_coset = {x**i * xx**j: i for i in xrange(5) for j in xrange((v-1)/5)}
                for c in to_coset:
                    if c == 1 or c == r or c == r2:
                        continue
                    if len(set(to_coset[elt] for elt in (r-1, c*(r-1), c-1, c-r, c-r**2))) == 5:
                        if existence:
                            return True
                        B = [one,r,r**2,c,c*r,c*r**2]
                        D = [[x**(i*5) * b for b in B] for i in xrange(t)]
                        break

    # Twin prime powers construction (see :wikipedia:`Difference_set`)
    #
    # i.e. v = p(p+2) where p and p+2 are prime powers
    #      k = (v-1)/2
    #      lambda = (k-1)/2
    elif (len(factorization) == 2 and
          abs(pow(*factorization[0])-pow(*factorization[1])) == 2 and
          k == (v-1)//2 and
          (l is None or 2*l == (v-1)//2-1)):

        # A difference set can be built from the set of elements
        # (x,y) in GF(p) x GF(p+2) such that:
        #
        # - either y=0
        # - x and y with x and y     squares
        # - x and y with x and y non-squares
        if existence:
            return True

        from sage.rings.finite_rings.constructor import FiniteField
        from sage.categories.cartesian_product import cartesian_product
        from itertools import product
        p,q = pow(*factorization[0]), pow(*factorization[1])
        if p>q:
            p,q=q,p
        Fp = FiniteField(p,'x')
        Fq = FiniteField(q,'x')
        Fpset = set(Fp)
        Fqset = set(Fq)
        Fp_squares = set(x**2 for x in Fpset)
        Fq_squares = set(x**2 for x in Fqset)

        # Pairs of squares, pairs of non-squares
        d = []
        d.extend(product(Fp_squares.difference([0]),Fq_squares.difference([0])))
        d.extend(product(Fpset.difference(Fp_squares),Fqset.difference(Fq_squares)))

        # All (x,0)
        d.extend((x,0) for x in Fpset)

        G = cartesian_product([Fp,Fq])
        D = [d]

    if D is None and are_hyperplanes_in_projective_geometry_parameters(v,k,l):
        _, (q,d) = are_hyperplanes_in_projective_geometry_parameters(v,k,l,True)
        if existence:
            return True
        else:
            G,D = singer_difference_set(q,d)

    if D is None:
        if existence:
            return Unknown
        raise NotImplementedError("No constructions for these parameters")

    if check and not is_difference_family(G,D,verbose=False):
        raise RuntimeError

    return G, D
