from sage.rings.ideal import Ideal_generic
from sage.misc.cachefunc import cached_method

from sage.structure.richcmp import op_EQ, op_NE, op_LT, op_GT, op_LE, op_GE

from sage.structure.element cimport Element
from sage.rings.tate_algebra_element cimport TateAlgebraTerm
from sage.rings.tate_algebra_element cimport TateAlgebraElement
from heapq import heappush, heappop


cdef _groebner_basis_buchberger(I, prec):
    r"""
    Compute a Gröbner basis of the Tate series ideal ``I`` using Buchberger's algorithm

    INPUT:
    
    - ``I`` - an ideal in a Tate series algebra

    - ``prec`` - the required precision of the calculation

    EXAMPLES::

    TODO
    
    """
    cdef list gb, rgb, indices, ts, S = [ ]
    cdef int i, j, l
    cdef TateAlgebraTerm ti, tj, t
    cdef TateAlgebraElement f, g, r, s
    cdef bint reduce = True

    gb = [ f.add_bigoh(prec) for f in I.gens() ]
    l = len(gb); indices = range(l)

    # We minimize the family of generators
    rgb = gb[:]
    i = 0
    while i < len(rgb):
        ti = (<TateAlgebraElement>rgb[i])._terms_c()[0]
        for j in range(len(rgb)):
            tj = (<TateAlgebraElement>rgb[j])._terms_c()[0]
            if j != i and tj._divides_c(ti, True):
                del rgb[i]
                del indices[i]
                break
        else:
            i += 1

    # We compute the initial S-polynomials
    for i in range(l):
        ti = (<TateAlgebraElement>gb[i])._terms_c()[0]
        for j in range(i+1, l):
            tj = (<TateAlgebraElement>gb[j])._terms_c()[0]
            if not ti.is_coprime_with(tj):
                s = (<TateAlgebraElement>gb[i])._Spoly_c(<TateAlgebraElement>gb[j])
                if not s.is_zero():
                    t = s._terms_c()[0]
                    heappush(S, (t._valuation_c(), t._exponent, i, j, s))

    # Main loop of Buchberger algorithm
    while S:
        # We reduce the Grobner basis if needed
        if reduce:
            reduce = False
            for i in range(len(rgb)-1, -1, -1):
                g = rgb[i]
                rgb[i] = g._positive_lshift_c(1)
                rgb[i] = g._mod_c(rgb)
                gb[indices[i]] = rgb[i]

        # We pop a new S-polynomial
        _, _, i, j, f = heappop(S)
        if gb[i] is None or gb[j] is None:
            continue
        r = f._mod_c(rgb)
        if r.is_zero():
            continue

        # We add it to our Grobner basis
        tj = r._terms_c()[0]
        j = len(gb)
        for i in range(j):
            g = gb[i]
            if g is None: continue
            ti = g._terms_c()[0]
            if not ti.is_coprime_with(tj):  # first Buchberger criterium
                s = g._Spoly_c(r)
                if not s.is_zero():
                    t = s._terms_c()[0]
                    heappush(S, (t._valuation_c(), t._exponent, i, j, s))
        gb.append(r)

        # We minimize the Grobner basis
        i = 0
        while i < len(rgb):
            ti = (<TateAlgebraElement>rgb[i])._terms_c()[0]
            if tj._divides_c(ti, True):
                if indices[i] >= l:
                    gb[indices[i]] = None
                del rgb[i]
                del indices[i]
            else:
                i += 1
        rgb.append(r)
        indices.append(j)
        # and reduce it
        reduce = True

    return rgb


class TateAlgebraIdeal(Ideal_generic):
    r"""
    Initialize a class for ideals in a Tate series algebra

    EXAMPLES::

    
    
    """
    
    @cached_method
    def groebner_basis(self, prec=None, algorithm=None):
        if prec is None:
            prec = self.ring().precision_cap()
        if algorithm is None:
            algorithm = "buchberger"
        if algorithm == "buchberger":
            return _groebner_basis_buchberger(self, prec)
        else:
            raise NotImplementedError("Only Buchberger algorithm is implemented so far")

    def _contains_(self, x):
        rgb = self.groebner_basis()
        return (x % rgb).is_zero()

    def _contains_ideal(self, I):
        rgb = self.groebner_basis()
        for f in I.gens():
            if not f in self:
                return False
        return True

    def _richcmp_(self, other, op):
        if not isinstance(other, TateAlgebraIdeal):
            raise NotImplementedError
        if op == op_GT:
            return self._contains_ideal(other) and not other._contains_ideal(self)
        elif op == op_GE:
            return self._contains_ideal(other)
        elif op == op_EQ:
            return self._contains_ideal(other) and other._contains_ideal(self)
        elif op == op_NE:
            return not(self._contains_ideal(other) and other._contains_ideal(self))
        elif op == op_LE:
            return other._contains_ideal(self)
        elif op == op_LT:
            return other._contains_ideal(self) and not self._contains_ideal(other)

    def is_saturated(self):
        if self.ring().base_ring().is_field():
            return True
        gb = self.groebner_basis()
        for g in gb:
            if g.valuation() > 0:
                return False
        return True

    def saturate(self):
        gb = self.groebner_basis()
        gens = [ g.monic() for g in gb ]
        return self.ring().ideal(gens)
