"""
Exterior algebras backend

This contains the backend implementations in Cython for the exterior algebra.
"""

from sage.libs.gmp.mpz cimport mpz_sizeinbase, mpz_setbit, mpz_tstbit, mpz_cmp_si, mpz_sgn
from sage.data_structures.bitset_base cimport bitset_t, bitset_init, bitset_first, bitset_next, bitset_set_to
from sage.structure.parent cimport Parent

cdef inline Integer bitset_to_int(FrozenBitset X):
    """
    Convert ``X`` to an :class:`Integer`.
    """
    cdef Integer ret = Integer(0)
    cdef long elt = bitset_first(X._bitset)
    while elt >= 0:
        mpz_setbit(ret.value, elt)
        elt = bitset_next(X._bitset, elt + 1)
    return ret

cdef inline FrozenBitset int_to_bitset(Integer n):
    """
    Convert a nonnegative integer ``n`` to a :class:`FrozenBitset`.
    """
    cdef unsigned long i

    if mpz_sgn(n.value) == 0:
        return FrozenBitset()

    cdef FrozenBitset ret = <FrozenBitset> FrozenBitset()
    cdef size_t s = mpz_sizeinbase(n.value, 2)
    bitset_init(ret._bitset, s)
    for i in range(s):
        bitset_set_to(ret._bitset, i, mpz_tstbit(n.value, i))
    return ret


cdef inline unsigned long degree(FrozenBitset X):
    """
    Compute the degree of ``X``.
    """
    cdef unsigned long ret = 0
    cdef long elt = bitset_first(X._bitset)
    while elt >= 0:
        ret += 1
        elt = bitset_next(X._bitset, elt + 1)
    return ret


#TODO: Bring the exterior algebra elements as a cdef class and make f know its type!
cdef inline FrozenBitset leading_supp(f):
    """
    Return the leading support of the exterior algebra element ``f``.
    """
    cdef dict mc = <dict> f._monomial_coefficients
    return int_to_bitset(min(bitset_to_int(k) for k in mc))

def leading_support(f):
    return leading_supp(f)

cpdef tuple get_leading_supports(tuple I):
    """
    Return the leading supports of the elements in ``I``.

    Used for testing mostly

    INPUT:

    - ``I`` -- a tuple of elements of an exterior algebra
    """
    # We filter out any elements that are 0
    return tuple(set([leading_supp(f) for f in I if f._monomial_coefficients]))


cdef inline build_monomial(E, supp):
    """
    Helper function for the fastest way to build a monomial.
    """
    return E.element_class(E, {supp: (<Parent> E)._base.one()})

cdef inline partial_S_poly(f, g, E, int side):
    """
    Compute one half of the `S`-polynomial for ``f`` and ``g``.

    This computes:

    .. MATH::

        LCM(LM(f), LM(g)) / LT(f) \cdot f.
    """
    cdef FrozenBitset lmf = leading_supp(f)
    cdef FrozenBitset lmg = leading_supp(g)
    cdef FrozenBitset D = <FrozenBitset> lmg.difference(lmf)
    if side == 0:
        ret = build_monomial(E, D) * f
        return (~ret[lmf._union(lmg)]) * ret
    ret = f * build_monomial(E, D)
    return ret * (~ret[lmf._union(lmg)])

cdef inline set preprocessing(list P, list G, E, int side):
    """
    Perform the preprocessing step.
    """
    #print("Start preprocessing:", P)
    cdef set L = set(partial_S_poly(f0, f1, E, side) for f0,f1 in P)
    L.update(partial_S_poly(f1, f0, E, side) for f0,f1 in P)
    if side == 2:
        # in partial_S_poly, side == 2 gets treated like right (== 1)
        L.update(partial_S_poly(f0, f1, E, 0) for f0,f1 in P)
        L.update(partial_S_poly(f1, f0, E, 0) for f0,f1 in P)

    cdef set done = set(leading_supp(f) for f in L)
    cdef set monL = set()
    for f in L:
        monL.update(f.support())
    monL.difference_update(done)

    while monL:
        m = int_to_bitset(max(bitset_to_int(k) for k in monL))
        done.add(m)
        monL.remove(m)
        for g in G:
            lm = leading_supp(g)
            if lm <= m:
                f = build_monomial(E, m.difference(lm)) * g
                if f in L:
                    break
                monL.update(set(f.support()) - done)
                L.add(f)
                break
    #print("preprocessing:", L)
    return L

cdef inline list reduction(list P, list G, E, int side):
    """
    Perform the reduction of ``P`` mod ``G`` in ``E``.
    """
    cdef set L = preprocessing(P, G, E, side)
    cdef Py_ssize_t i
    from sage.matrix.constructor import matrix
    M = matrix({(i, bitset_to_int(<FrozenBitset> m)): c for i,f in enumerate(L) for m,c in f._monomial_coefficients.items()},
               sparse=True)
    M.echelonize()  # Do this in place
    lead_supports = set(leading_supp(f) for f in L)
    return [E.element_class(E, {int_to_bitset(Integer(j)): c for j,c in M[i].iteritems()})
            for i,p in enumerate(M.pivots())
            if int_to_bitset(Integer(p)) not in lead_supports]

def compute_groebner(I, side):
    """
    Compute the reduced ``side`` Gröbner basis for the ideal ``I``.

    INPUT:

    - ``I`` -- the ideal
    - ``side`` -- integer; the side of the ideal: ``0`` for left, ``1`` for
      right, and ``2`` for two-sided
    """
    E = I.ring()
    cdef FrozenBitset p0, p1
    cdef unsigned long deg
    cdef Py_ssize_t i, j, k

    cdef list G = [f for f in I.gens() if f]  # Remove 0s TODO: We should make this unnecessary here
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
        p0 = leading_supp(f0)
        for j in range(i+1, n):
            f1 = G[j]
            p1 = leading_supp(f1)
            deg = degree(<FrozenBitset> (p0._union(p1)))
            if deg in P:
                P[deg].append((f0, f1))
            else:
                P[deg] = [(f0, f1)]

    while P:
        #print("Cur G:", G)
        Pp = P.pop(min(P))  # The selection: lowest lcm degree
        Gp = reduction(Pp, G, E, side)
        #print("Reduction yielded:", Gp)
        G.extend(Gp)
        for j in range(n, len(G)):
            f1 = G[j]
            p1 = leading_supp(f1)
            for i in range(j):
                f0 = G[i]
                p0 = leading_supp(f0)
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
        lm = leading_supp(g)
        did_reduction = True
        while did_reduction:
            supp = f.support()
            did_reduction = False
            for s in supp:
                if lm <= s:
                    did_reduction = True
                    mon = E.monomial(s - lm)
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

