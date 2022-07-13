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
        self.E = <Parent> I.ring()
        if I.side() == "left":
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
        return self.int_to_bitset(min(self.bitset_to_int(k) for k in mc))

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

    cdef inline set preprocessing(self, list P, list G):
        """
        Perform the preprocessing step.
        """
        #print("Start preprocessing:", P)
        cdef CliffordAlgebraElement f, g

        cdef set L = set()
        if self.side == 0:
            L.update(self.partial_S_poly_left(f0, f1) for f0,f1 in P)
            L.update(self.partial_S_poly_left(f1, f0) for f0,f1 in P)
        elif self.side == 1:
            L.update(self.partial_S_poly_right(f0, f1) for f0,f1 in P)    
            L.update(self.partial_S_poly_right(f1, f0) for f0,f1 in P)    
        if self.side == 2:
            L.update(self.partial_S_poly_left(f0, f1) for f0,f1 in P)
            L.update(self.partial_S_poly_left(f1, f0) for f0,f1 in P)
            L.update(self.partial_S_poly_right(f0, f1) for f0,f1 in P)
            L.update(self.partial_S_poly_right(f1, f0) for f0,f1 in P)

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
        #print("preprocessing:", L)
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
        cdef list G = [f for f in self.ideal.gens() if f]  # Remove 0s TODO: We should make this unnecessary here

        cdef Py_ssize_t n = len(G)
        cdef dict P = {}
        cdef list Gp

        # for ideals generated by homogeneous (wrt Z_2-grading) polynomials, we can just consider it as a left ideal
        # TODO: We can reduce the number of S-poly computations for Z_2-graded homogeneous
        #   ideals by throwing out those such that LCM(LM(f), LM(g)) == LM(f) * LM(g).
        if all(f.is_super_homogeneous() for f in G):
            side = 0

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
            #print("Cur G:", G)
            Pp = P.pop(min(P))  # The selection: lowest lcm degree
            Gp = self.reduction(Pp, G)
            #print("Reduction yielded:", Gp)
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

        #print(G)

        # Now that we have a Gröbner basis, we make this into a reduced Gröbner basis
        cdef set pairs = set((i, j) for i in range(n) for j in range(n) if i != j)
        cdef list supp
        cdef bint did_reduction
        cdef FrozenBitset lm, s
        while pairs:
            i,j = pairs.pop()
            # We perform the classical reduction algorithm here on each pair
            # TODO: Make this faster by using the previous technique
            f = G[i]
            g = G[j]
            lm = self.leading_supp(g)
            did_reduction = True
            while did_reduction:
                supp = f.support()
                did_reduction = False
                for s in supp:
                    if lm <= s:
                        did_reduction = True
                        mon = self.E.monomial(s - lm)
                        if side == 0:
                            gp = mon * g
                            f = f - f[s] / gp[s] * gp
                        else:
                            gp = g * mon
                            f = f - f[s] / gp[s] * gp
                        break
            if G[i] != f:
                G[i] = f
                #print("reduction:", G)
                if not f:
                    pairs.difference_update((k, i) for k in range(n))
                else:
                    pairs.update((k, i) for k in range(n) if k != i)

        return tuple([f for f in G if f])


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
        cdef Py_ssize_t i

        if mpz_sgn(n.value) == 0:
            return FrozenBitset()

        cdef FrozenBitset ret = <FrozenBitset> FrozenBitset()
        cdef size_t s = mpz_sizeinbase(n.value, 2)
        bitset_init(ret._bitset, s)
        for i in range(s):
            bitset_set_to(ret._bitset, i, mpz_tstbit(n.value, i))
        return ret


