r"""
Difference families

.. TODO::

    There is a slightly more general version where the stabilizers of the blocks
    are taken into account. A block is *short* if the stabilizer is not trivial.
    In the case is not trivial we consider slight less differences. It is still
    possible to construct BIBD from this more general version (see the
    Handbook of combinatorial designs).

REFERENCES:

.. [We72] R. M. Wilson "Cyclotomy and difference families in elementary Abelian
   groups", J. of Num. Th., 4 (1972), pp. 17-47.
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

def zmod_or_finite_field(q):
    from sage.rings.finite_rings.integer_mod_ring import Zmod
    from sage.rings.finite_rings.constructor import GF

    if arith.is_prime(q):
        return Zmod(q)
    return GF(q,'z')

def group_law(G):
    import operator
    from sage.categories.groups import Groups
    from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
    from sage.categories.modules import Modules

    if G in Groups():  # multiplicative groups
        return (G.one(), operator.mul, operator.inv)
    elif G in CommutativeAdditiveGroups() or G in Modules(ZZ):  # additive groups
        return (G.zero(), operator.add, operator.neg)
    else:
        raise ValueError("%s does not seem to be a group"%G)


def cyclotomic_cosets(K, e, x=None):
    r"""
    Return the `e`-th cyclotomic cosets on `K`.

    Let `e` be an integer that divides the cardinality of `K`. The `e`-th
    cyclotomic cosets are the cosets modulo the group of `e`-th power in `K^*`.
    They form a ``(q,f,f-1)`` difference family in `K`.

    INPUT:

    - ``K`` -- a finite field or an integer which is a prime power

    - ``e`` -- an integer that divides the cardinality of ``K`` minus one

    - ``x`` -- an optional generator for the multiplicative group `K - \{0\}`.

    OUTPUT:

    The output is a pair ``(K,H)`` where ``K`` is a finite field and ``H`` is
    the list of cosets.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import cyclotomic_cosets, is_difference_family
        sage: K,H = cyclotomic_cosets(7,2)
        sage: H
        [[1, 2, 4], [3, 6, 5]]
        sage: is_difference_family(K,H,7,3,2)
        True

        sage: K,H = cyclotomic_cosets(16,5)
        sage: H
        [[1, z^2 + z, z^2 + z + 1],
         [z, z^3 + z^2, z^3 + z^2 + z],
         [z^2, z^3 + z + 1, z^3 + z^2 + z + 1],
         [z^3, z^2 + 1, z^3 + z^2 + 1],
         [z + 1, z^3 + z, z^3 + 1]]
        sage: is_difference_family(K,H,16,3,2)
        True
    """
    if isinstance(K, (int,Integer)):
        q = K
        K = zmod_or_finite_field(K)
    else:
        q = K.cardinality()
    assert q%e == 1
    f = (q-1) // e
    if x is None:
        x = K.multiplicative_generator()
    return K, [[x**(e*s+i) for s in xrange(f)] for i in xrange(e)]

def quadratic_residues(K):
    r"""
    Return the set of quadratic residues in the finite field `K`.

    If `q` is congruent to `3` mod `4` they form a `(q, (q-1)//2, (q-1)//4)`-difference family.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import quadratic_residues, is_difference_family
        sage: K = GF(19)
        sage: quad_res = quadratic_residues(K)
        sage: is_difference_family(K,[quad_res],19,9,4)
        True

        sage: K = GF(23)
        sage: quad_res = quadratic_residues(K)
        sage: is_difference_family(K, [quad_res],23,11,5)
        True
    """
    quad_res = set(x**2 for x in K)
    quad_res.remove(K.zero())
    return sorted(quad_res)

def quartic_residues(K):
    r"""
    Return the set of quartic residues in the finite field `K`.

    If `q = 4t^2 + 1` with `t` odd, then it is a `(q, (q-1)//4, (q-5)//16)`-difference family.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import quartic_residues, is_difference_family

        sage: K = GF(37)   # 4*3^2 + 1
        sage: B = quartic_residues(K)
        sage: is_difference_family(K,[B],37,9,2)
        True

        sage: K = GF(101)   # 4*5^2 + 1
        sage: B = quartic_residues(K)
        sage: is_difference_family(K,[B],101,25,6)
        True
    """
    quart_res = set(x**4 for x in K)
    quart_res.remove(K.zero())
    return sorted(quart_res)

def is_difference_family(G, D, v=None, k=None, l=None, verbose=False):
    r"""
    Check wether ``D`` forms a difference family in ``G``.

    INPUT:

    - ``G`` - Abelian group of cardinality ``v``

    - ``D`` - a set of ``k``-subsets of ``G``

    - ``v``, ``k`` and ``l`` - optional parameters of the difference family

    - ``verbose`` - whether to print additional information

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import is_difference_family
        sage: G = Zmod(21)
        sage: D = [[0,1,4,14,16]]
        sage: is_difference_family(G, D, 21, 5)
        True

        sage: G = Zmod(41)
        sage: D = [[0,1,4,11,29],[0,2,8,17,21]]
        sage: is_difference_family(G, D, verbose=True)
        the element 28 in G is obtained more than 1 times
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

    Examples with a (non Abelian) multiplicative group::

        sage: G = DihedralGroup(8)
        sage: x,y = G.gens()
        sage: D1 = [[1,x,x^4], [1,x^2, y*x], [1,x^5,y], [1,x^6,y*x^2], [1,x^7,y*x^5]]
        sage: is_difference_family(G, D1, 16, 3, 2)
        True
    """
    if v is None:
        v = G.cardinality()
    else:
        v = int(v)
        if G.cardinality() != v:
            if verbose:
                print "G must have cardinality v (=%d)"%v
            return False

    if k is None:
        k = len(D[0])
    else:
        k = int(k)

    b = len(D)

    if any(len(d) != k for d in D):
        if verbose:
            print "each element of D must have cardinality k (=%d)"%k
        return False

    if l is None:
        if (b*k*(k-1)) % (v-1) != 0:
            if verbose:
                print "bk(k-1) is not 0 mod (v-1)"
            return False
        l = b*k*(k-1) // (v-1)
    else:
        l = int(l)

    if b*k*(k-1) != l*(v-1):
        if verbose:
            print "the relation bk(k-1) == l(v-1) is not satisfied with b=%d, k=%d, l=%d, v=%d"%(b,k,l,v)
        return False

    identity, op, inv = group_law(G)

    # now we check that every non-identity element of G occurs exactly l-time
    # as a difference
    counter = {g: 0 for g in G}
    del counter[identity]

    for d in D:
        dd = map(G,d)
        for i in xrange(k):
            for j in xrange(k):
                if i == j:
                    continue
                g = op(dd[i], inv(dd[j]))
                counter[g] += 1
                if counter[g] > l:
                    if verbose:
                        print "the element %s in G is obtained more than %s times"%(g,l)
                    return False
    if verbose:
        print "It is a ({},{},{})-difference family".format(v,k,l)
    return True


def _have_distinct_images(elts, f, images=None):
    r"""
    Check whether ``elts`` are non zero and mapped to distinct image through
    ``f``. If ``images`` is provided, the image must also be distinct from it.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_family import _have_distinct_images
        sage: f = {1:0, 2:1, 3:1, 4:2}
        sage: [(i,j) for i in srange(5) for j in srange(5) if _have_distinct_images([i,j], f)]
        [(1, 2), (1, 3), (1, 4), (2, 1), (2, 4), (3, 1), (3, 4), (4, 1), (4, 2), (4, 3)]

        sage: _have_distinct_images([2,4], f, set([0]))
        True
        sage: _have_distinct_images([2,4], f, set([2]))
        False
    """
    if images is None:
        images = set()
    for elt in elts:
        if elt.is_zero():
            return False
        i = f[elt]
        if i in images:
            return False
        images.add(i)
    return True

def difference_family(v, k, l=1, existence=False, check=True):
    r"""
    Return a `(k,l)`-difference family on a group of size `v`.

    Let `G` be a finite Abelian group. A *`(G,k,\lambda)`-difference family* is
    a collection of `k`-subsets of `G`, `D = \{D_1, D_2, \ldots, D_b\}` such
    that the set of differences seen as a multi-set contains each element of `G
    \backslash \{0\}` exactly `\lambda`-times.

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
        [31, 151, 181, 211, ...,  3061, 3121, 3181]
        sage: l6[Unknown]
        [61, 121]
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

    Other constructions for `l > 1`::

        sage: _ = designs.difference_family(23,11,5)
        sage: _ = designs.difference_family(101,25,6)
    """
    from database import DF_constructions

    if (l*(v-1)) % (k*(k-1)) != 0:
        if existence:
            return False
        raise EmptySetError("A (v,%d,%d)-difference family may exist only if %d*(v-1) = mod %d"%(k,l,l,k*(k-1)))

    D = None
    e = k*(k-1)
    t = (v-1) // e  # number of blocks

    if v in DF_constructions and (k,l) in DF_constructions[v]:
        return DF_constructions[v][(k,l)]()

    if arith.is_prime_power(v):
        K = zmod_or_finite_field(v)

        if l == (k-1):
            if existence:
                return True
            return K, cyclotomic_cosets(K, (v-1)//k)

        if v%4 == 3 and k == (v-1)//2 and l == (v-1)//4:
            if existence:
                return True
            return K, [quadratic_residues(K)]

        if v%8 == 5 and k == (v-1)//4 and l == (v-5)//16 and arith.is_square((v-1)//4):
            if existence:
                return True
            return K, [quartic_residues(K)]

        if l != 1:
            raise NotImplementedError

        one = K.one()
        x = K.multiplicative_generator()

        # Wilson (1972), Theorem 9
        if k%2 == 1:
            m = (k-1) // 2
            xx = x**m
            to_coset = {x**i * xx**j: i for i in xrange(m) for j in xrange((v-1)/m)}
            r = x ** ((v-1) // k)  # primitive k-th root of unity
            cosets = set()
            if _have_distinct_images((r**j-one for j in xrange(1,m+1)), to_coset):
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
            cosets = set([0])
            if _have_distinct_images((r**j-one for j in xrange(1,m)), to_coset, set([0])):
                if existence:
                    return True
                B = [K.zero()] + [r**j for j in xrange(k-1)]
                D = [[x**(i*m) * b for b in B] for i in xrange(t)]

        # Wilson (1972), Theorem 11
        if D is None and k == 6:
            r = x**((v-1)//3)  # primitive cube root of unity
            xx = x**5
            to_coset = {x**i * xx**j: i for i in xrange(5) for j in xrange((v-1)/5)}
            for c in to_coset:
                if _have_distinct_images((r-1, c*(r-1), c-1, c-r, c-r**2), to_coset):
                    if existence:
                        return True
                    B = [one,r,r**2,c,c*r,c*r**2]
                    D = [[x**(i*5) * b for b in B] for i in xrange(t)]
                    break

    if D is None:
        if existence:
            return Unknown
        raise NotImplementedError("No constructions for these parameters")

    if check and not is_difference_family(K,D,verbose=False):
        raise RuntimeError

    return K, D


