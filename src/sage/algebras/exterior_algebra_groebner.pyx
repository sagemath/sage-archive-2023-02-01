"""
Exterior algebras Gröbner bases

This contains the backend implementations in Cython for the Gröbner bases
of exterior algebra.

AUTHORS:

- Trevor Karn, Travis Scrimshaw (July 2022): Initial implementation
"""

#*****************************************************************************
#       Copyright (C) 2022 Trevor Karn <karnx018 at umn.edu>
#                 (C) 2022 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.gmp.mpz cimport mpz_sizeinbase, mpz_setbit, mpz_tstbit, mpz_cmp_si, mpz_sgn
from sage.data_structures.bitset_base cimport (bitset_t, bitset_init, bitset_first,
                                               bitset_next, bitset_set_to, bitset_len)
from sage.structure.parent cimport Parent

cdef inline long degree(FrozenBitset X):
    """
    Compute the degree of ``X``.
    """
    return bitset_len(X._bitset)


cdef inline CliffordAlgebraElement build_monomial(Parent E, FrozenBitset supp):
    """
    Helper function for the fastest way to build a monomial.
    """
    return <CliffordAlgebraElement> E.element_class(E, {supp: E._base.one()})

cdef class GroebnerStrategy:
    """
    A strategy for computing a Gröbner basis.

    INPUT:

    - ``I`` -- the ideal
    """
    def __init__(self, I):
        """
        Initialize ``self``.
        """
        self.ideal = I
        self.groebner_basis = (None,)
        self.E = <Parent> I.ring()
        self.homogeneous = all(x.is_super_homogeneous() for x in I.gens())
        if self.homogeneous or I.side() == "left":
            self.side = 0
        elif I.side() == "right":
            self.side = 1
        else:
            self.side = 2

    cdef inline FrozenBitset leading_supp(self, CliffordAlgebraElement f):
        """
        Return the leading support of the exterior algebra element ``f``.
        """
        cdef dict mc = <dict> f._monomial_coefficients
        return self.int_to_bitset(max(self.bitset_to_int(k) for k in mc))

    cdef inline partial_S_poly_left(self, CliffordAlgebraElement f, CliffordAlgebraElement g):
        """
        Compute one half of the `S`-polynomial for ``f`` and ``g``.

        This computes:

        .. MATH::

            LCM(LM(f), LM(g)) / LT(f) \cdot f.
        """
        cdef FrozenBitset lmf = self.leading_supp(f)
        cdef FrozenBitset lmg = self.leading_supp(g)
        cdef FrozenBitset D = <FrozenBitset> lmg.difference(lmf)
        ret = build_monomial(self.E, D) * f
        return (~ret[lmf._union(lmg)]) * ret

    cdef inline partial_S_poly_right(self, CliffordAlgebraElement f, CliffordAlgebraElement g):
        """
        Compute one half of the `S`-polynomial for ``f`` and ``g``.

        This computes:

        .. MATH::

            f \cdot LCM(LM(f), LM(g)) / LT(f).
        """
        cdef FrozenBitset lmf = self.leading_supp(f)
        cdef FrozenBitset lmg = self.leading_supp(g)
        cdef FrozenBitset D = <FrozenBitset> lmg.difference(lmf)
        ret = f * build_monomial(self.E, D)
        return ret * (~ret[lmf._union(lmg)])

    cdef inline bint build_S_poly(self, CliffordAlgebraElement f, CliffordAlgebraElement g):
        r"""
        Check to see if we should build the `S`-polynomial.

        For homogeneous ideals, we throw out all those pairs `(f, g)` such that

        .. MATH::

            LCM(LM(f), LM(g)) == LM(f) \cdot LM(g).
        """
        if not self.homogeneous:
            return True

        return (<FrozenBitset> self.leading_supp(f).intersection(self.leading_supp(g))).isempty()

    cdef inline set preprocessing(self, list P, list G):
        """
        Perform the preprocessing step.
        """
        cdef CliffordAlgebraElement f, g, f0, f1

        cdef set L = set()
        if self.side == 0:
            for f0, f1 in P:
                if self.build_S_poly(f0, f1):
                    L.add(self.partial_S_poly_left(f0, f1))
                    L.add(self.partial_S_poly_left(f1, f0))
        elif self.side == 1:
            for f0, f1 in P:
                if self.build_S_poly(f0, f1):
                    L.add(self.partial_S_poly_right(f0, f1) for f0,f1 in P)
                    L.add(self.partial_S_poly_right(f1, f0) for f0,f1 in P)  
        if self.side == 2:
            for f0, f1 in P:
                if self.build_S_poly(f0, f1):
                    L.add(self.partial_S_poly_left(f0, f1))
                    L.add(self.partial_S_poly_left(f1, f0))
                    L.add(self.partial_S_poly_right(f0, f1))
                    L.add(self.partial_S_poly_right(f1, f0))

        cdef set done = set(self.leading_supp(f) for f in L)
        cdef set monL = set()
        for f in L:
            monL.update(f._monomial_coefficients)
        monL.difference_update(done)

        while monL:
            m = self.int_to_bitset(max(self.bitset_to_int(k) for k in monL))
            done.add(m)
            monL.remove(m)
            for g in G:
                lm = self.leading_supp(g)
                if lm <= m:
                    f = <CliffordAlgebraElement> build_monomial(self.E, <FrozenBitset> m.difference(lm)) * g
                    if f in L:
                        break
                    monL.update(set(f._monomial_coefficients) - done)
                    L.add(f)
                    break
        return L

    cdef inline list reduction(self, list P, list G):
        """
        Perform the reduction of ``P`` mod ``G`` in ``E``.
        """
        cdef set L = self.preprocessing(P, G)
        cdef Py_ssize_t i
        from sage.matrix.constructor import matrix
        M = matrix({(i, self.bitset_to_int(<FrozenBitset> m)): c
                    for i,f in enumerate(L)
                    for m,c in (<CliffordAlgebraElement> f)._monomial_coefficients.items()},
                   sparse=True)
        M.echelonize()  # Do this in place
        lead_supports = set(self.leading_supp(<CliffordAlgebraElement> f) for f in L)
        return [self.E.element_class(self.E, {self.int_to_bitset(Integer(j)): c for j,c in M[i].iteritems()})
                for i,p in enumerate(M.pivots())
                if self.int_to_bitset(Integer(p)) not in lead_supports]

    def compute_groebner(self):
        """
        Compute the reduced ``side`` Gröbner basis for the ideal ``I``.
        """
        cdef FrozenBitset p0, p1
        cdef long deg
        cdef Py_ssize_t i, j, k
        cdef list G = [f for f in self.ideal.gens() if f]  # Remove 0s

        cdef Py_ssize_t n = len(G)
        cdef dict P = {}
        cdef list Gp

        for i in range(n):
            f0 = G[i]
            p0 = self.leading_supp(f0)
            for j in range(i+1, n):
                f1 = G[j]
                p1 = self.leading_supp(f1)
                deg = degree(<FrozenBitset> (p0._union(p1)))
                if deg in P:
                    P[deg].append((f0, f1))
                else:
                    P[deg] = [(f0, f1)]

        while P:
            Pp = P.pop(min(P))  # The selection: lowest lcm degree
            Gp = self.reduction(Pp, G)
            G.extend(Gp)
            for j in range(n, len(G)):
                f1 = G[j]
                p1 = self.leading_supp(f1)
                for i in range(j):
                    f0 = G[i]
                    p0 = self.leading_supp(f0)
                    deg = degree(<FrozenBitset> (p0._union(p1)))
                    if deg in P:
                        P[deg].append((f0, f1))
                    else:
                        P[deg] = [(f0, f1)]
            n = len(G)

        # Now that we have a Gröbner basis, we make this into a reduced Gröbner basis
        cdef set pairs = set((i, j) for i in range(n) for j in range(n) if i != j)
        cdef tuple supp
        cdef bint did_reduction
        cdef FrozenBitset lm, s
        while pairs:
            i,j = pairs.pop()
            # We perform the classical reduction algorithm here on each pair
            # TODO: Make this faster by using the previous technique?
            f = self.reduce_single(G[i], G[j])
            if G[i] != f:
                G[i] = f
                if not f:
                    pairs.difference_update((k, i) for k in range(n))
                else:
                    pairs.update((k, i) for k in range(n) if k != i)

        self.groebner_basis = tuple([~f[self.leading_supp(f)] * f for f in G if f])

    cpdef CliffordAlgebraElement reduce(self, CliffordAlgebraElement f):
        """
        Reduce ``f`` modulo the ideal with Gröbner basis ``G``.
        """
        for g in self.groebner_basis:
            f = self.reduce_single(f, <CliffordAlgebraElement> g)
        return f

    cdef CliffordAlgebraElement reduce_single(self, CliffordAlgebraElement f, CliffordAlgebraElement g):
        """
        Reduce ``f`` by ``g``.

        .. TODO::

            Optimize this by doing it in-place and changing the underlying dict of ``f``.
        """
        cdef FrozenBitset lm, s
        cdef tuple supp

        lm = self.leading_supp(g)
        did_reduction = True
        while did_reduction:
            supp = tuple(f._monomial_coefficients)
            did_reduction = False
            for s in supp:
                if lm <= s:
                    did_reduction = True
                    mon = self.E.monomial(s - lm)
                    if self.side == 0:
                        gp = mon * g
                        f -= f[s] / gp[s] * gp
                    else:
                        gp = g * mon
                        f -= f[s] / gp[s] * gp
                    break
        return f


    cdef Integer bitset_to_int(self, FrozenBitset X):
        raise NotImplementedError

    cdef FrozenBitset int_to_bitset(self, Integer n):
        raise NotImplementedError

    def sorted_monomials(self, as_dict=False):
        """
        Helper function to display the monomials in their term order
        used by ``self``.

        EXAMPLES::

            sage: from sage.algebras.exterior_algebra_groebner import *
            sage: E.<x,y,z> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([x, y])
            sage: GroebnerStrategyNegLex(I).sorted_monomials()
            [1, x, y, x*y, z, x*z, y*z, x*y*z]
            sage: GroebnerStrategyDegLex(I).sorted_monomials()
            [1, x, y, z, x*y, x*z, y*z, x*y*z]
            sage: GroebnerStrategyDegRevLex(I).sorted_monomials()
            [1, z, y, x, y*z, x*z, x*y, x*y*z]

            sage: E.<a,b,c,d> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([a, b])
            sage: GroebnerStrategyNegLex(I).sorted_monomials()
            [1,
             a,
             b, a*b,
             c, a*c, b*c, a*b*c,
             d, a*d, b*d, a*b*d, c*d, a*c*d, b*c*d, a*b*c*d]
            sage: GroebnerStrategyDegLex(I).sorted_monomials()
            [1,
             a, b, c, d,
             a*b, a*c, a*d, b*c, b*d, c*d,
             a*b*c, a*b*d, a*c*d, b*c*d,
             a*b*c*d]
            sage: GroebnerStrategyDegRevLex(I).sorted_monomials()
            [1,
             d, c, b, a,
             c*d, b*d, a*d, b*c, a*c, a*b,
             b*c*d, a*c*d, a*b*d, a*b*c,
             a*b*c*d]
        """
        cdef FrozenBitset X
        cdef Integer i
        cdef list D = [self.bitset_to_int(X) for X in self.E._indices]
        D.sort()
        if as_dict:
            return {i: build_monomial(self.E, self.int_to_bitset(i)) for i in D}
        return [build_monomial(self.E, self.int_to_bitset(i)) for i in D]

    def monomial_to_int(self):
        """
        Helper function to display the monomials in their term order
        used by ``self``.

        EXAMPLES::

            sage: from sage.algebras.exterior_algebra_groebner import *
            sage: E.<a,b,c,d> = ExteriorAlgebra(QQ)
            sage: I = E.ideal([a, b])
            sage: GroebnerStrategyDegLex(I).monomial_to_int()
            {1: 0,
             a: 1, b: 2, c: 3, d: 4,
             a*b: 5, a*c: 6, a*d: 7, b*c: 8, b*d: 9, c*d: 10,
             a*b*c: 11, a*b*d: 12, a*c*d: 13, b*c*d: 14,
             a*b*c*d: 15}
            sage: GroebnerStrategyDegRevLex(I).monomial_to_int()
            {1: 0,
             a: 4, b: 3, c: 2, d: 1,
             a*b: 10, a*c: 9, a*d: 7, b*c: 8, b*d: 6, c*d: 5,
             a*b*c: 14, a*b*d: 13, a*c*d: 12, b*c*d: 11,
             a*b*c*d: 15}
        """
        B = self.E.basis()
        return {B[X]: self.bitset_to_int(X) for X in self.E._indices}

cdef class GroebnerStrategyNegLex(GroebnerStrategy):
    """
    Gröbner basis strategy implementing neglex ordering.
    """
    cdef inline Integer bitset_to_int(self, FrozenBitset X):
        """
        Convert ``X`` to an :class:`Integer`.
        """
        cdef Integer ret = Integer(0)
        cdef long elt = bitset_first(X._bitset)
        while elt >= 0:
            mpz_setbit(ret.value, elt)
            elt = bitset_next(X._bitset, elt + 1)
        return ret

    cdef inline FrozenBitset int_to_bitset(self, Integer n):
        """
        Convert a nonnegative integer ``n`` to a :class:`FrozenBitset`.
        """
        cdef size_t i

        if mpz_sgn(n.value) == 0:
            return FrozenBitset()

        cdef FrozenBitset ret = <FrozenBitset> FrozenBitset()
        cdef size_t s = mpz_sizeinbase(n.value, 2)
        bitset_init(ret._bitset, s)
        for i in range(s):
            bitset_set_to(ret._bitset, i, mpz_tstbit(n.value, i))
        return ret

cdef class GroebnerStrategyDegRevLex(GroebnerStrategy):
    """
    Gröbner basis strategy implementing degree revlex ordering.
    """
    def __init__(self, I):
        """
        Initialize ``self``.
        """
        GroebnerStrategy.__init__(self, I)
        self.rank = Integer(self.E.ngens())

    cdef inline Integer bitset_to_int(self, FrozenBitset X):
        """
        Convert ``X`` to an :class:`Integer`.
        """
        if X.isempty():
            return Integer(0)

        cdef Integer n = self.rank
        cdef long i, deg = degree(X)
        cdef long r = 1
        cdef Integer t = Integer(0)

        cdef long elt = bitset_first(X._bitset)
        while elt >= 0:
            t += Integer(elt).binomial(r)
            r += 1
            elt = bitset_next(X._bitset, elt + 1)
        return Integer(sum(n.binomial(i) for i in range(deg+1)) - t - 1)

    cdef inline FrozenBitset int_to_bitset(self, Integer n):
        """
        Convert a nonnegative integer ``n`` to a :class:`FrozenBitset`.
        """
        cdef size_t i

        if mpz_sgn(n.value) == 0:
            return FrozenBitset()

        cdef Py_ssize_t deg = 0
        cdef Integer binom = Integer(1)
        while n >= binom:
            n -= binom
            deg += 1
            binom = self.rank.binomial(deg)

        # TODO: Cythonize the from_rank
        from sage.combinat.combination import from_rank
        return FrozenBitset([self.rank - val - 1 for val in from_rank(n, self.rank, deg)])

cdef class GroebnerStrategyDegLex(GroebnerStrategy):
    """
    Gröbner basis strategy implementing degree lex ordering.
    """
    def __init__(self, I):
        """
        Initialize ``self``.
        """
        GroebnerStrategy.__init__(self, I)
        self.rank = Integer(self.E.ngens())

    cdef inline Integer bitset_to_int(self, FrozenBitset X):
        """
        Convert ``X`` to an :class:`Integer`.
        """
        if X.isempty():
            return Integer(0)

        cdef Integer n = self.rank
        cdef long i, deg = degree(X)
        cdef long r = deg
        cdef Integer t = Integer(0)

        cdef long elt = bitset_first(X._bitset)
        while elt >= 0:
            t += Integer(n - 1 - elt).binomial(r)
            r -= 1
            elt = bitset_next(X._bitset, elt + 1)
        return Integer(sum(n.binomial(i) for i in range(deg+1)) - t - 1)

    cdef inline FrozenBitset int_to_bitset(self, Integer n):
        """
        Convert a nonnegative integer ``n`` to a :class:`FrozenBitset`.
        """
        cdef size_t i

        if mpz_sgn(n.value) == 0:
            return FrozenBitset()

        cdef Py_ssize_t deg = 0
        cdef Integer binom = Integer(1)
        while n >= binom:
            n -= binom
            deg += 1
            binom = self.rank.binomial(deg)

        from sage.combinat.combination import from_rank
        return FrozenBitset(from_rank(n, self.rank, deg))

