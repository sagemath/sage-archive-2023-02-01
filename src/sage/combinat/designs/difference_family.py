r"""
Difference families

REFERENCES:

.. [Wi72] R. M. Wilson "Cyclotomy and difference families in elementary Abelian
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

def is_difference_family(G, D, v=None, k=None, l=None, verbose=False):
    r"""
    Check wether ``D`` forms a difference family in ``G``.

    INPUT:

    - ``G`` - Abelian group of cardinality ``v``

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

    The function also supports multiplicative groups (non necessarily Abelian)::

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
                if g == identity:
                    if verbose:
                        print "two identical elements in the same block"
                    return False
                counter[g] += 1
                if counter[g] > l:
                    if verbose:
                        print "the element %s in G is obtained more than %s times"%(g,l)
                    return False
    if verbose:
        print "It is a ({},{},{})-difference family".format(v,k,l)
    return True


def difference_family(v, k, l=1, existence=False, check=True):
    r"""
    Return a (``k``, ``l``)-difference family on a group of size ``v``.

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
         4: (3,2)
         5: (4,3)
         7: (3,2), (6,5)
         8: (7,6)
         9: (4,3), (8,7)
        11: (5,2), (5,4)
        13: (3,2), (4,3), (6,5)
        16: (3,2), (5,4)
        17: (4,3), (8,7)
        19: (3,2), (6,5), (9,4), (9,8)
        25: (3,2), (4,3), (6,5), (8,7)
        29: (4,3), (7,6)
        31: (3,2), (5,4), (6,5)
        37: (3,2), (4,3), (6,5), (9,2), (9,8)
        41: (4,3), (5,4), (8,7)
        43: (3,2), (6,5), (7,6)
        49: (3,2), (4,3), (6,5), (8,7)
        53: (4,3)
        61: (3,2), (4,3), (5,4), (6,5)
        64: (3,2), (7,6), (9,8)
        67: (3,2), (6,5)
        71: (5,4), (7,6)
        73: (3,2), (4,3), (6,5), (8,7), (9,8)
        79: (3,2), (6,5)
        81: (4,3), (5,4), (8,7)
        89: (4,3), (8,7)
        97: (3,2), (4,3), (6,5), (8,7)

    TESTS:

    Check more of the Wilson constructions from [Wi72]_::

        sage: Q5 = [241, 281,421,601,641, 661, 701, 821,881]
        sage: Q9 = [73, 1153, 1873, 2017]
        sage: Q15 = [76231]
        sage: Q4 = [13, 73, 97, 109, 181, 229, 241, 277, 337, 409, 421, 457]
        sage: Q8 = [1009, 3137, 3697]
        sage: for Q,k in [(Q4,4),(Q5,5),(Q8,8),(Q9,9),(Q15,15)]:
        ....:     for q in Q:
        ....:         assert designs.difference_family(q,k,1,existence=True)
        ....:         _ = designs.difference_family(q,k,1)

    .. TODO::

        There is a slightly more general version of difference families where
        the stabilizers of the blocks are taken into account. A block is *short*
        if the stabilizer is not trivial. The more general version is called a
        *partial difference family*. It is still possible to construct BIBD from
        this more general version (see the chapter 16 in the Handbook
        [DesignHandbook]_).

        Implement recursive constructions from Buratti "Recursive for difference
        matrices and relative difference families" (1998) and Jungnickel
        "Composition theorems for difference families and regular planes" (1978)
    """
    if (l*(v-1)) % (k*(k-1)) != 0:
        if existence:
            return False
        raise EmptySetError("A (v,%d,%d)-difference family may exist only if %d*(v-1) = mod %d"%(k,l,l,k*(k-1)))

    from database import DF_constructions
    if (v,k,l) in DF_constructions:
        return DF_constructions[(v,k,l)]()

    e = k*(k-1)
    t = l*(v-1) // e  # number of blocks

    D = None

    if arith.is_prime_power(v):
        from sage.rings.finite_rings.constructor import GF
        K = GF(v,'z')
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

    if D is None:
        if existence:
            return Unknown
        raise NotImplementedError("No constructions for these parameters")

    if check and not is_difference_family(K,D,verbose=False):
        raise RuntimeError

    return K, D


