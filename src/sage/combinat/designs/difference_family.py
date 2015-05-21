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
    :func:`difference_family` | Return a (``k``, ``l``)-difference family on an Abelian group of size ``v``.
    :func:`singer_difference_set` | Return a difference set associated to hyperplanes in a projective space.
    :func:`df_q_6_1` | Return a difference set with parameter `k=6` on a finite field.
    :func:`radical_difference_set` | Return a radical difference set
    :func:`radical_difference_family` | Return a radical difference family
    :func:`twin_prime_powers_difference_set` | Return a twin prime powers difference family.
    :func:`wilson_1972_difference_family` | Return a radical difference family on a finite field following a construction of Wilson.

REFERENCES:

.. [Bu95] M. Buratti "On simple radical difference families", J. of
   Combinatorial Designs, vol. 3, no. 2 (1995)

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

    TESTS::

        sage: K = GF(3^2,'z')
        sage: z = K.gen()
        sage: D = [[1,z+1,2]]
        sage: _ = is_difference_family(K, D, verbose=True)
        the number of differences (=6) must be a multiple of v-1=8
        sage: _
        False
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
                print "the number of differences (={}) must be a multiple of v-1={}".format(nb_diff,v-1)
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

def df_q_6_1(K, existence=False, check=True):
    r"""
    Return a `(q,6,1)`-difference family over the finite field `K`.

    The construction uses Theorem 11 of [Wi72]_.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import is_difference_family, df_q_6_1
        sage: prime_powers = [v for v in xrange(31,500,30) if is_prime_power(v)]
        sage: parameters = [v for v in prime_powers if df_q_6_1(GF(v,'a'), existence=True)]
        sage: print parameters
        [31, 151, 181, 211, 241, 271, 331, 361, 421]
        sage: for v in parameters:
        ....:     K = GF(v, 'a')
        ....:     df = df_q_6_1(K, check=True)
        ....:     assert is_difference_family(K, df, v, 6, 1)

    .. TODO:

        Do improvements due to Zhen and Wu 1999.
    """
    v = K.cardinality()
    x = K.multiplicative_generator()
    one = K.one()
    if v % 30 != 1:
        if existence:
            return False
        raise EmptySetError("k(k-1)=30 should divide (v-1)")
    t = (v-1) // 30  # number of blocks

    r = x**((v-1)//3)  # primitive cube root of unity
    r2 = r*r           # the other primitive cube root

    # we now compute the cosets of x**i
    xx = x**5
    to_coset = {x**i * xx**j: i for i in xrange(5) for j in xrange((v-1)/5)}

    for c in to_coset: # the loop runs through all nonzero elements of K
        if c == one or c == r or c == r2:
            continue
        if len(set(to_coset[elt] for elt in (r-one, c*(r-one), c-one, c-r, c-r**2))) == 5:
            if existence:
                return True
            B = [one,r,r2,c,c*r,c*r2]
            D = [[xx**i * b for b in B] for i in xrange(t)]
            break
    else:
        if existence:
            return Unknown
        raise NotImplementedError("Wilson construction failed for v={}".format(v))

    if check and not is_difference_family(K, D, v, 6, 1):
        raise RuntimeError("Wilson 1972 construction failed! Please e-mail sage-devel@googlegroups.com")

    return D

def radical_difference_set(K, k, l=1, existence=False, check=True):
    r"""
    Return a difference set made of a cyclotomic coset in the finite field
    ``K`` and with paramters ``k`` and ``l``.

    Most of these difference sets appear in chapter VI.18.48 of the Handbook of
    combinatorial designs.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import radical_difference_set

        sage: D = radical_difference_set(GF(7), 3, 1); D
        [[1, 2, 4]]
        sage: sorted(x-y for x in D[0] for y in D[0] if x != y)
        [1, 2, 3, 4, 5, 6]

        sage: D = radical_difference_set(GF(16,'a'), 6, 2)
        sage: sorted(x-y for x in D[0] for y in D[0] if x != y)
        [1,
         1,
         a,
         a,
         a + 1,
         a + 1,
         a^2,
         a^2,
         ...
         a^3 + a^2 + a + 1,
         a^3 + a^2 + a + 1]

        sage: for (v,k,l) in [(3,2,1), (7,3,1), (7,4,2), (11,5,2), (11,6,3),
        ....:                 (13,4,1), (16,6,2), (19,9,4), (19,10,5)]:
        ....:
        ....:     assert radical_difference_set(GF(v,'a'), k, l, existence=True), "pb with v={} k={} l={}".format(v,k,l)
    """
    v = K.cardinality()
    one = K.one()
    x = K.multiplicative_generator()

    if l*(v-1) != k*(k-1):
        if existence:
            return False
        raise EmptySetError("l*(v-1) is not equal to k*(k-1)")

    # trivial case
    if (v-1) == k:
        if existence:
            return True
        return K.cyclotomic_cosets(x, [one])

    # q = 3 mod 4
    elif v%4 == 3 and k == (v-1)//2:
        if existence:
            return True
        D = K.cyclotomic_cosets(x**2, [one])

    # q = 3 mod 4
    elif v%4 == 3 and k == (v+1)//2:
        if existence:
            return True
        D = K.cyclotomic_cosets(x**2, [one])
        D[0].insert(0, K.zero())

    # q = 4t^2 + 1, t odd
    elif v%8 == 5 and k == (v-1)//4 and arith.is_square((v-1)//4):
        if existence:
            return True
        D = K.cyclotomic_cosets(x**4, [one])

    # q = 4t^2 + 9, t odd
    elif v%8 == 5 and k == (v+3)//4 and arith.is_square((v-9)//4):
        if existence:
            return True
        D = K.cyclotomic_cosets(x**4, [one])
        D[0].insert(0,K.zero())

    # one case with k-1 = (v-1)/3
    elif (v,k,l) == (16,6,2):
        if existence:
            return True
        D = K.cyclotomic_cosets(x**3, [one])
        D[0].insert(0,K.zero())

    # one case with k = (v-1)/8
    elif (v,k,l) == (73,9,1):
        if existence:
            return True
        D = K.cyclotomic_cosets(x**8, [one])

    else:
        if existence:
            return Unknown
        raise NotImplementedError("no radical difference set is "
                "implemented for the parameters (v,k,l) = ({},{},{}".format(v,k,l))

    if check and not is_difference_family(K, D, v, k, l):
        raise RuntimeError("Sage tried to build a cyclotomic coset with "
                "parameters ({},{},{}) but it seems that it failed! Please "
                "e-mail sage-devel@googlegroups.com".format(v,k,l))

    return D

def wilson_1972_difference_family(K, k, existence=False, check=True):
    r"""
    Wilson construction of difference families on finite field.
    
    The construction appears in [Wi72]_.

    INPUT:

    - ``K`` -- a finite field

    - ``k`` -- an integer such that `k(k-1)` divides `v-1` where `v` is the
      cardinality of the finite field `K`

    - ``existence`` -- if ``True`` then return either ``True`` if the
      construction is possible or ``False`` if it can not.

    - ``check`` -- boolean (default: ``True``) whether to check that the output
      is valid. This should not be needed, but it guarantee that the output is
      correct. It might be faster to set it to ``False``.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import wilson_1972_difference_family

        sage: wilson_1972_difference_family(GF(13),4)
        [[0, 1, 3, 9]]

        sage: wilson_1972_difference_family(GF(337), 7, existence=True)
        True

        sage: from itertools import ifilter, count
        sage: ppap = lambda k,r: ifilter(is_prime_power, (k*i+r for i in count(1)))


        sage: it4 = ppap(4*3,1)
        sage: for _ in range(7):
        ....:     v = next(it4)
        ....:     existence = wilson_1972_difference_family(GF(v,'a'), 4, existence=True)
        ....:     print "v = {}: {}".format(v, existence)
        v = 13: True
        v = 25: True
        v = 37: False
        v = 49: True
        v = 61: False
        v = 73: True
        v = 97: True

        sage: it5 = ppap(5*4,1)
        sage: for _ in range(7):
        ....:     v = next(it5)
        ....:     existence = wilson_1972_difference_family(GF(v,'a'), 5, existence=True)
        ....:     print "v = {:3}: {}".format(v, existence)
        v =  41: True
        v =  61: True
        v =  81: False
        v = 101: False
        v = 121: False
        v = 181: False
        v = 241: True

        sage: it6 = ppap(6*5,1)
        sage: for _ in range(7):
        ....:     v = next(it6)
        ....:     existence = wilson_1972_difference_family(GF(v,'a'), 6, existence=True)
        ....:     print "v = {:3}: {}".format(v, existence)
        v =  31: False
        v =  61: False
        v = 121: False
        v = 151: False
        v = 181: True
        v = 211: True
        v = 241: True

        sage: it7 = ppap(7*6,1)
        sage: for _ in range(7):
        ....:     v = next(it7)
        ....:     existence = wilson_1972_difference_family(GF(v,'a'), 7, existence=True)
        ....:     print "v = {:3}: {}".format(v, existence)
        v =  43: False
        v = 127: False
        v = 169: False
        v = 211: False
        v = 337: True
        v = 379: False
        v = 421: True
    """
    v = K.cardinality()
    zero = K.zero()
    one = K.one()
    x = K.multiplicative_generator()
    e = k*(k-1)

    if (v-1) % e != 0:
        if existence:
            return False
        raise EmptySetError("In the wilson construction, k(k-1) must divide (v-1) but k={} and v={}".format(k,v))

    t = (v-1) // (k*(k-1))

    if k%2 == 1:
        m = (k-1) // 2
        r = x ** ((v-1) // k)     # k-th root of unity
        xx = x ** m
        A = set(r**i - one for i in range(1,m+1))
    else:
        m = k // 2
        r = x ** ((v-1) // (k-1)) # (k-1)-th root of unity
        xx = x ** m
        A = set(r**i - one for i in range(1,m))
        A.add(one)

    # now, we check whether the elements of A belong to distinct cosets modulo
    # H^m where H = K \ {0}
    AA = set(y/x for x in A for y in A if x != y)
    xxi = one
    for i in range((v-1)/m):
        if xxi in AA:
                if existence:
                    return False
                raise EmptySetError("In the Wilson construction with v={} "
                "and k={}, the roots of unity fail to belong to distinct "
                "cosets modulo H^m")
        xxi *= xx
    
    if existence:
        return True

    D = K.cyclotomic_cosets(r, [xx**i for i in xrange(t)])
    if k%2 == 0:
        for d in D:
            d.insert(0, zero)

    if check and not is_difference_family(K,D,v,k,1):
        raise RuntimeError("Wilson construction of difference family "
                           "failed with parameters v={} and k={}. "
                           "Please contact sage-devel@googlegroups.com".format(v,k))

    return D


def radical_difference_family(K, k, l=1, existence=False, check=True):
    r"""
    Return a ``(v,k,l)``-radical difference family.

    Let `K` be a finite field. A *radical difference family* is a difference
    family on `K` whose blocks are made of either cyclotomic cosets or
    cyclotomic cosets together with `0`. The terminology comes from
    M. Buratti article [Bu95]_ but the constructions go back to R. Wilson
    [Wi72]_.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import radical_difference_family

        sage: radical_difference_family(GF(73),9)
        [[1, 2, 4, 8, 16, 32, 37, 55, 64]]

        sage: radical_difference_family(GF(281),5)
        [[1, 86, 90, 153, 232],
         [4, 50, 63, 79, 85],
         [5, 36, 149, 169, 203],
         [7, 40, 68, 219, 228],
         [9, 121, 212, 248, 253],
         [29, 81, 222, 246, 265],
         [31, 137, 167, 247, 261],
         [32, 70, 118, 119, 223],
         [39, 56, 66, 138, 263],
         [43, 45, 116, 141, 217],
         [98, 101, 109, 256, 279],
         [106, 124, 145, 201, 267],
         [111, 123, 155, 181, 273],
         [156, 209, 224, 264, 271]]

    .. TODO:

        Implement the more general Buratti construction from [Bu95]_
    """
    v = K.cardinality()
    x = K.multiplicative_generator()
    one = K.one()
    e = k*(k-1)
    if (l*(v-1)) % e:
        raise ValueError("k (k-1) = {} should be a multiple of l (v-1) ={}".format(
                         k*(k-1), l*(v-1)))
    t = l*(v-1) // e  # number of blocks

    if t == 1:
        return radical_difference_set(K, k, l, existence=existence, check=check)

    elif l == (k-1):
        if existence:
            return True
        else:
            return K.cyclotomic_cosets(x**((v-1)//k))[1:]

    # all the other cases below concern the case l == 1
    elif l != 1:
        if existence:
            return Unknown
        raise NotImplementedError("No radical families implemented for l > 2")

    # Wilson (1972), Theorem 9 and 10
    elif wilson_1972_difference_family(K,k,existence=True):
        if existence:
            return True
        D = wilson_1972_difference_family(K,k,check=False)

    elif existence:
      return Unknown

    else:
      raise NotImplementedError("Sage does not know how to build a radical "
               "difference family with parameters v={} and k={}".format(v,k))

    if check and not is_difference_family(K, D, v, k, l):
        raise RuntimeError("radical_difference_family produced a wrong "
                           "difference family with parameters v={}, "
                           "k={}, l={}. Please contact "
                           "sage-devel@googlegroups.com".format(v,k,l))

    return D

def twin_prime_powers_difference_set(p, check=True):
    r"""
    Return a difference set on `GF(p) \times GF(p+2)`.

    The difference set is built from the following element of the cartesian
    product of finite fields `GF(p) \times GF(p+2)`:

    - `(x,0)` with any `x`
    - `(x,y)` with `x` and `y` squares
    - `(x,y)` with `x` and `y` non-squares

    For more information see :wikipedia:`Difference_set`.

    INPUT:

    - ``check`` -- boolean (default: ``True``). If ``True`` then the result of
      the computation is checked before being returned. This should not be
      needed but ensures that the output is correct.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import twin_prime_powers_difference_set
        sage: G,D = twin_prime_powers_difference_set(3)
        sage: G
        The cartesian product of (Finite Field of size 3, Finite Field of size 5)
        sage: D
        [[(1, 1), (1, 4), (2, 2), (2, 3), (0, 0), (1, 0), (2, 0)]]
    """
    from sage.rings.finite_rings.constructor import FiniteField
    from sage.categories.cartesian_product import cartesian_product
    from itertools import product
    Fp = FiniteField(p,'x')
    Fq = FiniteField(p+2,'x')
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

    if check and not is_difference_family(G, [d]):
        raise RuntimeError("twin_prime_powers_difference_set produced a wrong "
                           "difference set with p={}. Please contact "
                           "sage-devel@googlegroups.com".format(p))

    return G, [d]

def difference_family(v, k, l=1, existence=False, explain_construction=False, check=True):
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

    INPUT:

    - ``v,k,l`` -- parameters of the difference family. If ``l`` is not provided
      it is assumed to be ``1``.

    - ``existence`` -- if ``True``, then return either ``True`` if Sage knows
      how to build such design, ``Unknown`` if it does not and ``False`` if it
      knows that the design does not exist..

    - ``explain_construction`` -- instead of returning a difference family,
      returns a string that explains the construction used.

    - ``check`` -- boolean (default: ``True``). If ``True`` then the result of
      the computation is checked before being returned. This should not be
      needed but ensures that the output is correct.

    OUTPUT:

    A pair ``(G,D)`` made of a group `G` and a difference family `D` on that
    group. Or, if ``existence`` is ``True`` a troolean or if
    ``explain_construction`` is ``True`` a string.

    EXAMPLES::

        sage: G,D = designs.difference_family(73,4)
        sage: G
        Finite Field of size 73
        sage: D
        [[0, 1, 8, 64],
         [0, 2, 16, 55],
         [0, 3, 24, 46],
         [0, 25, 54, 67],
         [0, 35, 50, 61],
         [0, 36, 41, 69]]
        sage: print designs.difference_family(73, 4, explain_construction=True)
        Radical difference family on a finite field

        sage: G,D = designs.difference_family(15,7,3)
        sage: G
        The cartesian product of (Finite Field of size 3, Finite Field of size 5)
        sage: D
        [[(1, 1), (1, 4), (2, 2), (2, 3), (0, 0), (1, 0), (2, 0)]]
        sage: print designs.difference_family(15,7,3,explain_construction=True)
        Twin prime powers difference family

        sage: print designs.difference_family(91,10,1,explain_construction=True)
        Singer difference set

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

    List available constructions::

        sage: for v in xrange(2,100):
        ....:     constructions = []
        ....:     for k in xrange(2,10):
        ....:         for l in xrange(1,10):
        ....:             if designs.difference_family(v,k,l,existence=True):
        ....:                 constructions.append((k,l))
        ....:                 _ = designs.difference_family(v,k,l)
        ....:     if constructions:
        ....:         print "%2d: %s"%(v, ', '.join('(%d,%d)'%(k,l) for k,l in constructions))
         3: (2,1)
         4: (3,2)
         5: (2,1), (4,3)
         6: (5,4)
         7: (2,1), (3,1), (3,2), (4,2), (6,5)
         8: (7,6)
         9: (2,1), (4,3), (8,7)
        10: (9,8)
        11: (2,1), (4,6), (5,2), (5,4), (6,3)
        13: (2,1), (3,1), (3,2), (4,1), (4,3), (5,5), (6,5)
        15: (3,1), (4,6), (5,6), (7,3)
        16: (3,2), (5,4), (6,2)
        17: (2,1), (4,3), (5,5), (8,7)
        19: (2,1), (3,1), (3,2), (4,2), (6,5), (9,4), (9,8)
        21: (3,1), (4,3), (5,1), (6,3), (6,5)
        22: (4,2), (6,5), (7,4), (8,8)
        23: (2,1)
        25: (2,1), (3,1), (3,2), (4,1), (4,3), (6,5), (7,7), (8,7)
        27: (2,1), (3,1)
        28: (3,2), (6,5)
        29: (2,1), (4,3), (7,3), (7,6), (8,4), (8,6)
        31: (2,1), (3,1), (3,2), (4,2), (5,2), (5,4), (6,1), (6,5)
        33: (3,1), (5,5), (6,5)
        34: (4,2)
        35: (5,2)
        37: (2,1), (3,1), (3,2), (4,1), (4,3), (6,5), (9,2), (9,8)
        39: (3,1), (6,5)
        40: (3,2), (4,1)
        41: (2,1), (4,3), (5,1), (5,4), (6,3), (8,7)
        43: (2,1), (3,1), (3,2), (4,2), (6,5), (7,2), (7,3), (7,6), (8,4)
        45: (3,1), (5,1)
        46: (4,2), (6,2)
        47: (2,1)
        49: (2,1), (3,1), (3,2), (4,1), (4,3), (6,5), (8,7), (9,3)
        51: (3,1), (5,2), (6,3)
        52: (4,1)
        53: (2,1), (4,3)
        55: (3,1), (9,4)
        57: (3,1), (7,3), (8,1)
        59: (2,1)
        61: (2,1), (3,1), (3,2), (4,3), (5,1), (5,4), (6,2), (6,3), (6,5)
        63: (3,1)
        64: (3,2), (4,1), (7,2), (7,6), (9,8)
        65: (5,1)
        67: (2,1), (3,1), (3,2), (6,5)
        69: (3,1)
        71: (2,1), (5,2), (5,4), (7,3), (7,6), (8,4)
        73: (2,1), (3,1), (3,2), (4,1), (4,3), (6,5), (8,7), (9,1), (9,8)
        75: (3,1), (5,2)
        76: (4,1)
        79: (2,1), (3,1), (3,2), (6,5)
        81: (2,1), (3,1), (4,3), (5,1), (5,4), (8,7)
        83: (2,1)
        85: (4,1), (7,2), (7,3), (8,2)
        89: (2,1), (4,3), (8,7)
        91: (6,1)
        97: (2,1), (3,1), (3,2), (4,1), (4,3), (6,5), (8,7), (9,3)

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

    Check a failing construction (:trac:`17528`):

        sage: designs.difference_family(9,3)
        Traceback (most recent call last):
        ...
        NotImplementedError: No construction available for (9,3,1)-difference family

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
        elif explain_construction:
            return "The database contains a ({},{},{})-difference family".format(v,k,l)

        vv, blocks = next(DF[v,k,l].iteritems())

        # Build the group
        from sage.rings.finite_rings.integer_mod_ring import Zmod
        if len(vv) == 1:
            G = Zmod(vv[0])
        else:
            from sage.categories.cartesian_product import cartesian_product
            G = cartesian_product([Zmod(i) for i in vv])

        df = [[G(i) for i in b] for b in blocks]

        if check and not is_difference_family(G, df, v=v, k=k, l=l):
            raise RuntimeError("There is an invalid ({},{},{})-difference "
                    "family in the database... Please contact "
                    "sage-devel@googlegroups.com".format(v,k,l))

        return G,df

    e = k*(k-1)
    if (l*(v-1)) % e:
        if existence:
            return Unknown
        raise NotImplementedError("No construction available for ({},{},{})-difference family".format(v,k,l))
    t = l*(v-1) // e  # number of blocks

    # trivial construction
    if k == (v-1) and l == (v-2):
        from sage.rings.finite_rings.integer_mod_ring import Zmod
        G = Zmod(v)
        return G, [range(1,v)]

    factorization = arith.factor(v)
    D = None

    if len(factorization) == 1:  # i.e. is v a prime power
        from sage.rings.finite_rings.constructor import GF
        G = K = GF(v,'z')

        if radical_difference_family(K, k, l, existence=True):
            if existence:
                return True
            elif explain_construction:
                return "Radical difference family on a finite field"
            else:
                D = radical_difference_family(K,k,l)

        elif l == 1 and k == 6 and df_q_6_1(K,existence=True):
            if existence:
                return True
            elif explain_construction:
                return "Wilson 1972 difference family made from the union of two cyclotomic cosets"
            else:
                D = df_q_6_1(K)

    # Twin prime powers construction
    # i.e. v = p(p+2) where p and p+2 are prime powers
    #      k = (v-1)/2
    #      lambda = (k-1)/2 (ie 2l+1 = k)
    elif (k == (v-1)//2 and
          l == (k-1)//2 and
          len(factorization) == 2 and
          abs(pow(*factorization[0]) - pow(*factorization[1])) == 2):
        if existence:
            return True
        elif explain_construction:
            return "Twin prime powers difference family"
        else:
            p = pow(*factorization[0])
            q = pow(*factorization[1])
            if p > q:
                p,q = q,p
            G,D = twin_prime_powers_difference_set(p,check=False)

    if D is None and are_hyperplanes_in_projective_geometry_parameters(v,k,l):
        _, (q,d) = are_hyperplanes_in_projective_geometry_parameters(v,k,l,True)
        if existence:
            return True
        elif explain_construction:
            return "Singer difference set"
        else:
            G,D = singer_difference_set(q,d)

    if D is None:
        if existence:
            return Unknown
        raise NotImplementedError("No constructions for these parameters")

    if check and not is_difference_family(G,D,v=v,k=k,l=l,verbose=False):
        raise RuntimeError("There is a problem. Sage built the following "
                "difference family on G='{}' with parameters ({},{},{}):\n "
                "{}\nwhich seems to not be a difference family... "
                "Please contact sage-devel@googlegroups.com".format(G,v,k,l,D))

    return G, D
