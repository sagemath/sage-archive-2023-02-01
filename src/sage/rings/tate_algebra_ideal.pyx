from sage.rings.ideal import Ideal_generic
from sage.misc.cachefunc import cached_method

from sage.structure.richcmp import op_EQ, op_NE, op_LT, op_GT, op_LE, op_GE

from sage.structure.element cimport Element
from sage.rings.tate_algebra_element cimport TateAlgebraTerm
from sage.rings.tate_algebra_element cimport TateAlgebraElement
from heapq import heappush, heappop


cdef _groebner_basis_buchberger(I, prec):
    r"""
    Compute a Gröbner basis of the Tate algebra ideal I using Buchberger's algorithm

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

        sage: R = Zp(3,prec=10,print_mode="digits"); R
        3-adic Ring with capped relative precision 10
        sage: A.<x,y> = TateAlgebra(R); A
        Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
        sage: f = 3*x^2 + 5*x*y^2; f
        (...0000000012)*x*y^2 + (...00000000010)*x^2
        sage: g = 5*x^2*y + 3; g
        (...0000000012)*x^2*y + (...00000000010)
        sage: I = A.ideal([f,g]); I
        Ideal ((...0000000012)*x*y^2 + (...00000000010)*x^2, (...0000000012)*x^2*y + (...00000000010)) of Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
    
    """
    
    @cached_method
    def groebner_basis(self, prec=None, algorithm=None):
        r"""
        Compute a Gröbner basis of the ideal

        INPUT:

        - ``prec`` - (default: None) the precision required in the
          computations. If ``None``, defaults to the algebra precision cap

        - ``algorithm`` - (default: None) the algorithm to use in the
          calculations. Currently, only Buchberger's algorithm is implemented,
          and it is also the default.

        NOTE::

        The result of this method is cached.
        
        EXAMPLES::

            sage: R = Zp(3,prec=10,print_mode="digits"); R
            3-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
            sage: f = 3*x^2 + 5*x*y^2; f
            (...0000000012)*x*y^2 + (...00000000010)*x^2
            sage: g = 5*x^2*y + 3; g
            (...0000000012)*x^2*y + (...00000000010)
            sage: I = A.ideal([f,g]); I
            Ideal ((...0000000012)*x*y^2 + (...00000000010)*x^2, (...0000000012)*x^2*y + (...00000000010)) of Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
            sage: I.groebner_basis()
            [(...0000000012)*x*y^2 + (...0000000010)*x^2 + O(3^10),
             (...0000000012)*x^2*y + (...0000000010) + O(3^10),
             (...0000000010)*x^3 + (...2222222220)*y + O(3^10),
             (...0000000010)*y^2 + (...2101210200)*x + O(3^10)]
            sage: I.groebner_basis(algorithm="F4")
            Traceback (most recent call last):
            ...
            NotImplementedError: Only Buchberger algorithm is implemented so far


        """
        if prec is None:
            prec = self.ring().precision_cap()
        if algorithm is None:
            algorithm = "buchberger"
        if algorithm == "buchberger":
            return _groebner_basis_buchberger(self, prec)
        else:
            raise NotImplementedError("Only Buchberger algorithm is implemented so far")

    def _contains_(self, x):
        r"""
        Test whether `x` lies in the ideal

        INPUT:

        - ``x`` - a Tate series

        EXAMPLES::

            sage: R = Zp(3,prec=10,print_mode="digits"); R
            3-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
            sage: f = 3*x^2 + 5*x*y^2; f
            (...0000000012)*x*y^2 + (...00000000010)*x^2
            sage: g = 5*x^2*y + 3; g
            (...0000000012)*x^2*y + (...00000000010)
            sage: I = A.ideal([f,g]); I
            Ideal ((...0000000012)*x*y^2 + (...00000000010)*x^2, (...0000000012)*x^2*y + (...00000000010)) of Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
            sage: f in I # indirect doctest
            True
            sage: (f+g) in I # indirect doctest
            True
            sage: (f+1) in I # indirect doctest
            False

        
        """
        rgb = self.groebner_basis()
        return (x % rgb).is_zero()

    def _contains_ideal(self, I):
        r"""
        Test whether the ideal `I` is contained in the ideal

        INPUT:

        - ``I`` - an ideal in a Tate series algebra

        EXAMPLES::

            sage: R = Zp(3,prec=10,print_mode="digits"); R
            3-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
            sage: f = 3*x^2 + 5*x*y^2; f
            (...0000000012)*x*y^2 + (...00000000010)*x^2
            sage: g = 5*x^2*y + 3; g
            (...0000000012)*x^2*y + (...00000000010)
            sage: I = A.ideal([f,g]); I
            Ideal ((...0000000012)*x*y^2 + (...00000000010)*x^2, (...0000000012)*x^2*y + (...00000000010)) of Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
            sage: A.ideal([f]) < I # indirect doctest
            True
            sage: I < A.ideal([f]) # indirect doctest
            False
            sage: A.ideal([1]) < I # indirect doctest
            False
            sage: I < A.ideal([1]) # indirect doctest
            True
        
        """
        rgb = self.groebner_basis()
        for f in I.gens():
            if not f in self:
                return False
        return True

    def _richcmp_(self, other, op):
        r"""
        Compare the ideal with another according to the inclusion partial order

        INPUT:

        - ``other`` - an ideal in a Tate series algebra

        - ``op`` - a comparison operator

        EXAMPLES::

            sage: R = Zp(3,prec=10,print_mode="digits"); R
            3-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
            sage: f = 3*x^2 + 5*x*y^2; f
            (...0000000012)*x*y^2 + (...00000000010)*x^2
            sage: g = 5*x^2*y + 3; g
            (...0000000012)*x^2*y + (...00000000010)
            sage: I = A.ideal([f,g]); I
            Ideal ((...0000000012)*x*y^2 + (...00000000010)*x^2, (...0000000012)*x^2*y + (...00000000010)) of Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
            sage: A.ideal([f]) < I # indirect doctest
            True
            sage: I < A.ideal([f]) # indirect doctest
            False
            sage: A.ideal([1]) < I # indirect doctest
            False
            sage: I < A.ideal([1]) # indirect doctest
            True
            sage: I <= A.ideal([f,g])
            True
            sage: I == A.ideal([f,g])
            True
            sage: I <= A.ideal([f])
            False
            sage: A.ideal([f]) <= I
            True
            sage: A.ideal([f]) == I
            False
        
        """
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
        r"""
        Test whether the Tate series ideal is saturated

        The ideal `I` is saturated if and only it is generated by elements with
        0 valuation, or equivalently, if the quotient ideal `(I:\pi^\infty) =
        \{f \in A : \exists n \in \mathbb{N}, \pi^n f \in I\}` is exactly `I`,
        where `\pi` is the uniformizer of the discrete valuation ring.

        If the base ring is a field, all ideals are saturated.

        EXAMPLES::

        By default, Tate algebras are defined over a field, so this test is
        always true:
        
            sage: R = Zp(3,prec=10,print_mode="digits"); R
            3-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
            sage: f = 3*x^2 + 5*x*y^2; f
            (...0000000012)*x*y^2 + (...00000000010)*x^2
            sage: g = 5*x^2*y + 3; g
            (...0000000012)*x^2*y + (...00000000010)
            sage: A.ideal([f,g]).is_saturated()
            True
            sage: A.ideal([3*f,3*g]).is_saturated()
            True
            sage: A.base_ring()
            3-adic Field with capped relative precision 10

        The test is only relevant if the base ring is not a field:
        
            sage: AA = A.integer_ring()
            sage: II = AA.ideal([f,g])
            sage: II.is_saturated()
            False
            sage: II.groebner_basis()
            [(...0000000012)*x*y^2 + (...0000000010)*x^2 + O(3^10),
             (...0000000012)*x^2*y + (...0000000010) + O(3^10),
             (...0000000010)*x^3 + (...2222222220)*y + O(3^10),
             (...0000000010)*y^2 + (...2101210200)*x + O(3^10)]
            sage: AA.ideal([f]).is_saturated()
            True

        
        """
        if self.ring().base_ring().is_field():
            return True
        gb = self.groebner_basis()
        for g in gb:
            if g.valuation() > 0:
                return False
        return True

    def saturate(self):
        r"""
        Saturate the ideal by the uniformizer of the discrete valuation ring

        The result is the quotient ideal `(I:\pi^\infty) = \{f \in A : \exists n
        \in \mathbb{N}, \pi^n f \in I\}`, where `\pi` is the uniformizer of the
        discrete valuation ring.

        EXAMPLES::

            sage: R = Zp(3,prec=10,print_mode="digits"); R
            3-adic Ring with capped relative precision 10
            sage: A.<x,y> = TateAlgebra(R); A
            Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
            sage: f = 3*x^2 + 5*x*y^2; f
            (...0000000012)*x*y^2 + (...00000000010)*x^2
            sage: g = 5*x^2*y + 3; g
            (...0000000012)*x^2*y + (...00000000010)
            sage: I = A.ideal([f,g]); I
            Ideal ((...0000000012)*x*y^2 + (...00000000010)*x^2, (...0000000012)*x^2*y + (...00000000010)) of Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
            sage: I.saturate()
            Ideal ((...0000000001)*x*y^2 + (...1210121020)*x^2 + O(3^10), (...0000000001)*x^2*y + (...1210121020) + O(3^10), (...000000001)*x^3 + (...222222222)*y + O(3^9), (...000000001)*y^2 + (...210121020)*x + O(3^9)) of Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10

        Over a ring:
        
            sage: AA = A.integer_ring()
            sage: II = AA.ideal([f,g]); II
            Ideal ((...0000000012)*x*y^2 + (...00000000010)*x^2, (...0000000012)*x^2*y + (...00000000010)) of Integer ring of the Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
            sage: IIs = II.saturate(); IIs
            Ideal ((...0000000001)*x*y^2 + (...1210121020)*x^2 + O(3^10), (...0000000001)*x^2*y + (...1210121020) + O(3^10), (...000000001)*x^3 + (...222222222)*y + O(3^9), (...000000001)*y^2 + (...210121020)*x + O(3^9)) of Integer ring of the Tate Algebra in x (val >= 0), y (val >= 0) over 3-adic Field with capped relative precision 10
            sage: IIs.is_saturated()
            True

        
        
        """
        # TODO: Should we skip the GB computation + generation of another ideal
        # if the base ring is a field?
        gb = self.groebner_basis()
        gens = [ g.monic() for g in gb ]
        return self.ring().ideal(gens)
