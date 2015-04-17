"""
Augmented valuations polynomial rings

Implements inductive valuations as defined in [ML1936].

REFERENCES:

.. [ML1936] Mac Lane, S. (1936). A construction for prime ideals as absolute
values of an algebraic field. Duke Mathematical Journal, 2(3), 492-510.

.. [ML1936'] MacLane, S. (1936). A construction for absolute values in
polynomial rings. Transactions of the American Mathematical Society, 40(3),
363-395.

AUTHORS:

- Julian Rueth (15-04-2013): initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from developing_valuation import DevelopingValuation, _lift_to_maximal_precision

from sage.misc.cachefunc import cached_method
from sage.rings.all import infinity

class AugmentedValuation(DevelopingValuation):
    """
    An augmented valuation is a discrete valuation on a polynomial ring. It
    extends another discrete valuation `v` by setting the valuation of a
    polynomial `f` to the minumum of `v(f_i)i\mu` when writing `f=\sum_i
    f_i\mu^i`.

    INPUT:

    - ``v`` -- a discrete valuation on a polynomial ring

    - ``phi`` -- a key polynomial over ``v``

    - ``mu`` -- a rational number such that ``mu > v(phi)`` or ``infinity``

    - ``check`` -- a boolean (default: ``True``), whether to check that the
      parameters define an augemented value of ``v``

    EXAMPLES::

            sage: K.<u> = Qq(4,5)
            sage: R.<x> = K[]
            sage: v = GaussValuation(R)
            sage: v.extension(x, 1/2) # indirect doctest
            [ Gauss valuation induced by 2-adic valuation, v((1 + O(2^5))*x) = 1/2 ]

    """
    def __init__(self, v, phi, mu, check=True):
        """
        Initialization.

        TESTS::

            sage: from sage.rings.padics.augmented_valuation import AugmentedValuation
            sage: K.<u> = Qq(4,5)
            sage: R.<x> = K[]
            sage: v = GaussValuation(R)
            sage: w = AugmentedValuation(v, x, 1/2)
            sage: type(w)
            <class 'sage.rings.padics.augmented_valuation.AugmentedValuation'>

        """
        if phi.parent() is not v.domain():
            raise ValueError("phi must be in the domain of v")

        if check:
            is_key, reason = v.is_key(phi, explain=True)
            if not is_key:
                raise ValueError(reason)
            if mu <= v(phi):
                raise ValueError("the value of the key polynomial must strictly increase but `%s` does not exceed `%s`."%(phi, v(phi)))

        DevelopingValuation.__init__(self, v.domain(), phi)

        self._base_valuation = v
        self._mu = mu

        if check and not self.residue_field().has_coerce_map_from(v.residue_field()):
            raise ValueError("the residue field `%s` does not embed into `%s`"%(v.residue_field(), self.residue_field()))

    def __hash__(self):
        return hash(self._base_valuation, self._mu, self._phi)

    def _cache_key(self):
        return self._base_valuation, self._mu, self._phi

    def _call_(self, f):
        """
        Evaluate this valuation at ``f``.

        INPUT::

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: R = Qp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 3

            sage: v = v.extension( x^2 + x + 1, 1)
            sage: v(f)
            0
            sage: v(f * v.phi()^3 )
            3
            sage: v(S.zero())
            +Infinity

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation %s but is in %s"%(self.domain(),f.parent()))

        from sage.rings.all import infinity
        if f.is_zero():
            return infinity

        if f.degree() < self.phi().degree():
            return self._base_valuation(f)

        # We can slightly optimize the approach of DevelopingValuation._call_
        # We know that self(f) >= self._base_valuation(f)
        # as soon as we find a coefficient of f with self._base_valuation(c) ==
        # self._base_valuation(f) we know that this is the valuation of f

        # this optimization does only pay off for polynomials of large degree:
        if f.degree() // self.phi().degree() <= 3:
            return DevelopingValuation._call_(self, f)

        ret = infinity

        lower_bound = self._base_valuation(f)

        for v in self.valuations(f):
            ret = min(ret, v)
            if ret == lower_bound:
                return ret
        return ret

    @cached_method
    def Q(self):
        if self._mu is infinity or self.tau().is_zero():
            raise NotImplementedError

        return self.equivalence_unit(self._mu * self.tau())

    @cached_method
    def Q_(self):
        try:
            ret = self.equivalence_reciprocal(self.Q())

            assert self.is_equivalence_unit(ret)
            # esentially this checks that the reduction of Q'*phi^tau is the
            # generator of the residue field
            assert self._base_valuation.reduce(self.Q()*ret)(self.residue_field_generator()).is_one()

        except ValueError:
            print "CHEATING - HARD CODED RECIPROCAL"
            Q = self.Q()
            pi = Q.parent().base().constant_base_field().uniformizer()
            ret = Q/(pi**(self(Q)*2))

        assert self.is_equivalence_unit(ret)
        # esentially this checks that the reduction of Q'*phi^tau is the
        # generator of the residue field
        assert self._base_valuation.reduce(self.Q()*ret)(self.residue_field_generator()).is_one()

        return ret

    def equivalence_unit(self, s):
        """
        Return an equivalence unit of minimal degree and valuation ``s``.

        INPUT:

        - ``s`` -- a rational number

        OUTPUT:

        An element of the domain of this valuation

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.extension(x^2 + x + u, 1)
            sage: w.equivalence_unit(0)
            1 + O(2^5)
            sage: w.equivalence_unit(-4)
            2^-4 + O(2)

        Since an equivalence unit is of effective degree zero, `\phi` must not
        divide it. Therefore, its valuation is in the value group of the base
        valuation::

            sage: w = v.extension(x, 1/2)
            sage: w.equivalence_unit(3/2)
            Traceback (most recent call last):
            ...
            ValueError: v must be an integer
            sage: w.equivalence_unit(1)
            2 + O(2^6)

        An equivalence unit might not be integral, even if ``s >= 0``::

            sage: v = v.extension(x, 3/4)
            sage: w = v.extension(x^4 + 8, 5)
            sage: w.equivalence_unit(1/2)
            (2^-1 + O(2^4))*x^2

        .. SEEALSO::

            :meth:`effective_degree`, :meth:`is_equivalence_unit`

        """
        if s < 0 and not self.domain().base_ring().is_field():
            raise NotImplementedError("only implemented for polynomial rings over fields")

        ret = self._base_valuation.element_with_valuation(s)

        assert self.is_equivalence_unit(ret)
        assert self(ret) == s
        return ret

    def element_with_valuation(self, s):
        """
        Create an element of minimal degree and of valuation ``s``.

        INPUT:

        - ``s`` -- a rational number in the value group of this valuation

        OUTPUT:

        An element in the domain of this valuation

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.extension(x^2 + x + u, 1/2)
            sage: w.element_with_valuation(0)
            1 + O(2^5)
            sage: w.element_with_valuation(1/2)
            (1 + O(2^5))*x^2 + (1 + O(2^5))*x + u + O(2^5)
            sage: w.element_with_valuation(1)
            2 + O(2^6)
            sage: c = w.element_with_valuation(-1/2); c
            (2^-1 + O(2^4))*x^2 + (2^-1 + O(2^4))*x + u*2^-1 + O(2^4)
            sage: w(c)
            -1/2
            sage: w.element_with_valuation(1/3)
            Traceback (most recent call last):
            ...
            ValueError: s must be in the value group of the valuation

        """
        if s not in self.value_group():
            raise ValueError("s must be in the value group of the valuation")
        if s < 0 and not self.domain().base_ring().is_field():
            raise NotImplementedError("only implemented for polynomial rings over fields")

        ret = self.domain().one()
        while s not in self._base_valuation.value_group():
            ret *= self._phi
            s -= self._mu
        return ret * self._base_valuation.element_with_valuation(s)

    def _latex_(self):
        vals = [self]
        v = self
        while isinstance(v, AugmentedValuation):
            v = v._base_valuation
            vals.append(v)
        vals.reverse()
        from sage.misc.latex import latex
        vals = [ "v_%s(%s) = %s"%(i,latex(v._phi), latex(v._mu)) if isinstance(v, AugmentedValuation) else latex(v) for i,v in enumerate(vals) ]
        return "[ %s ]"%", ".join(vals)

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.extension(x^2 + x + u, 1/2)
            sage: w # indirect doctest
            [ Gauss valuation induced by 2-adic valuation, v((1 + O(2^5))*x^2 + (1 + O(2^5))*x + u + O(2^5)) = 1/2 ]

        """
        vals = [self]
        v = self
        while isinstance(v, AugmentedValuation):
            v = v._base_valuation
            vals.append(v)
        vals.reverse()
        vals = [ "v(%s) = %s"%(v._phi, v._mu) if isinstance(v, AugmentedValuation) else str(v) for v in vals ]
        return "[ %s ]"%", ".join(vals)

    @cached_method
    def constant_valuation(self):
        """
        Return the restriction of this valuation to the constants of the
        polynomial ring.

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.extension(x^2 + x + u, 1/2)
            sage: w.constant_valuation()
            2-adic valuation

        """
        return self._base_valuation.constant_valuation()

    @cached_method
    def value_group(self):
        """
        Return the value group of this valuation.

        OUTPUT:

        Currently, there is no support for additive subgroups of `\QQ`.
        Therefore we create this value group as a fractional ideal of `\QQ`.
        However, `\QQ` does not support fractional ideals, so we use fractional
        ideals of the trivial extensions `\QQ[x]/(x)`

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.extension(x, infinity)
            sage: w.value_group()
            Fractional ideal (1)

            sage: w = v.extension(x^2 + x + u, 1/2)
            sage: w.value_group()
            Fractional ideal (1/2)
            sage: w = w.extension((x^2 + x + u)^2 + 2, 5/3)
            sage: w.value_group()
            Fractional ideal (1/6)

        """
        base = self._base_valuation.value_group()
        if self._mu is infinity:
            return base
        return base + base.number_field().ideal(self._mu)

    def valuations(self, f):
        """
        Return the valuations of the `f_i\phi^i` in the expansion `f=\sum_i
        f_i\phi^i`.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        OUTPUT:

        An iterator over rational numbers, `[v(f_0), v(f_1\phi), \dots]`

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.extension(x^2 + x + u, 1/2)
            sage: list(w.valuations( x^2 + 1 ))
            [0, 1/2]
            sage: w = w.extension((x^2 + x + u)^2 + 2, 5/3)
            sage: list(w.valuations( ((x^2 + x + u)^2 + 2)^3 ))
            [+Infinity, +Infinity, +Infinity, 5]

        TESTS::

            sage: w = v.extension(x, infinity)
            sage: list(w.valuations( x^2 + 1 ))
            [0, +Infinity, +Infinity]

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of this valuation")

        if self._mu is infinity:
            num_infty_coefficients = f.degree() // self.phi().degree()
            yield self._base_valuation(self.coefficients(f).next())
            for i in range(num_infty_coefficients):
                yield infinity
        else:
            for i,c in enumerate(self.coefficients(f)):
                yield self._base_valuation(c) + i*self._mu

    def reduce(self, f):
        r"""
        Reduce ``f`` module this valuation.

        INPUT:

        - ``f`` -- an element of the domain of this valuation of non-negative
          valuation

        OUTPUT:

        an element of the :meth:`residue_ring` of this valuation, the reduction
        modulo the ideal of elements of positive valuation

        ALGORITHM:

        We follow the algorithm given in the proof of Theorem 12.1 of [ML1936]:
        If ``f`` has positive valuation, the reduction is simply zero.
        Otherwise, let `f=\sum f_i\phi^i` be the expansion of `f`, as computed
        by :meth:`coefficients`. Since the valuation is zero, the exponents `i`
        must all be multiples of :meth:`tau`.
        Hence, there is an :meth:`equivalence_unit` `Q` with the same valuation
        as `\phi^\tau`. Let `Q'` be its :meth:`reciprocal_inverse`.
        Now, rewrite each term `f_i\phi^{i\tau}=(f_iQ^i)(\phi^\tauQ^{-1})^i`;
        it turns out that the second factor in this expression is a lift of the
        generator of the :meth:`residue_field`. The reduction of the first
        factor can be computed recursively.

        EXAMPLES::

            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.reduce(x)
            x
            sage: v.reduce(S(u))
            u0

            sage: w = v.extension(x^2 + x + u, 1/2)
            sage: w.reduce(S.one())
            1
            sage: w.reduce(S(2))
            0
            sage: w.reduce(S(u))
            u0
            sage: w.reduce(x) # this gives the generator of the residue field extension of w over v
            u1
            sage: f = (x^2 + x + u)^2 / 2
            sage: w.reduce(f)
            x
            sage: w.reduce(f + x + 1)
            x + u1 + 1

            sage: w = w.extension((x^2 + x + u)^2 + 2, 5/3)
            sage: g = ((x^2 + x + u)^2 + 2)^3 / 2^5
            sage: w.reduce(g)
            x
            sage: w.reduce(f)
            1
            sage: w(f-1) > 0 # checking the above
            True
            sage: w.reduce(f * g)
            x
            sage: w.reduce(f + g)
            x + 1

        REFERENCES:

        .. [ML1936] Mac Lane, S. (1936). A construction for prime ideals as absolute
        values of an algebraic field. Duke Mathematical Journal, 2(3), 492-510.

        .. SEEALSO::

            :meth:`lift`

        """
        if f.parent() is not self.domain():
            raise ValueError("f must be in the domain of the valuation")
        if not self.domain().base_ring().is_field():
            raise NotImplementedError("only implemented for polynomial rings over fields")

        if self(f) < 0:
            assert self(f) < 0
            print self(f)
            raise ValueError("f must have non-negative valuation")
        elif self(f) > 0:
            return self.residue_ring().zero()

        # if this extends a trivial valuation, then this is very easy: just
        # return the constant coefficient in the phi-adic expansion; everything
        # else must have positive valuation
        if self._base_valuation.value_group() == 0:
            assert self.valuations(f).next() == 0
            if self.value_group() == 0:
                raise NotImplementedError
            return self.residue_ring()(self.coefficients(f).next())(self.residue_field_generator())

        # if this is an infinite valuation, then we can simply drop all but the
        # constant term
        if self._mu is infinity:
            return self.residue_ring()(self._base_valuation.reduce(self.coefficients(f).next())(self.residue_field().gen()))

        CV = zip(self.coefficients(f), self.valuations(f))
        # rewrite as sum of f_i phi^{i tau}, i.e., drop most coefficients
        assert not any([v==0 for i,(c,v) in enumerate(CV) if i % self.tau() != 0])
        CV = CV[::self.tau()]

        # replace f_i by f_i Q^{i tau}
        vQ = self._mu * self.tau()
        CV = [(c*self.Q()**i, v - vQ*i) for i,(c,v) in enumerate(CV)]
        assert all([self._base_valuation(c)>=0 for c,v in CV])

        # recursively reduce the f_i Q^{i tau}
        C = [self._base_valuation.reduce(c)(self.residue_field_generator()) for c,v in CV]

        # reduce the Q'^i phi^i
        return self.residue_ring()(C)

    def lift(self, F):
        """
        Return a polynomial whose :meth:`reduction` is ``F``.

        INPUT:

        - ``F`` -- an element of the :meth:`residue_ring`

        OUTPUT:

        a polynomial in the domain of the valuation with reduction ``F``, monic
        if ``F`` is monic

        ALGORITHM:

        Since this is the inverse of :meth:`reduce`, we only have to go backwards
        through the algorithm described there.

        EXAMPLES::

            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: y = v.residue_ring().gen()
            sage: u0 = v.residue_field().gen()
            sage: v.lift(y)
            (1 + O(2^10))*x

            sage: w = v.extension(x^2 + x + u, 1/2)
            sage: r = w.residue_ring()
            sage: y = r.gen()
            sage: u1 = w.residue_field().gen()
            sage: w.lift(r.one())
            1 + O(2^10)
            sage: w.lift(r.zero())
            0
            sage: w.lift(r(u0))
            u + O(2^10)
            sage: w.lift(r(u1))
            (1 + O(2^10))*x
            sage: w.reduce(w.lift(y)) == y
            True
            sage: w.reduce(w.lift(y + u1 + 1)) == y + u1 + 1
            True

            sage: w = w.extension((x^2 + x + u)^2 + 2, 5/3)
            sage: r = w.residue_ring()
            sage: y = r.gen()
            sage: u2 = w.residue_field().gen()
            sage: w.reduce(w.lift(y)) == y
            True
            sage: w.reduce(w.lift(r.one())) == 1
            True
            sage: w.reduce(w.lift(y + 1)) == y +  1
            True

        A more complicated example::

            sage: v = GaussValuation(S)
            sage: v = v.extension(x^2 + x + u, 1)
            sage: v = v.extension((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)

            sage: u = v.residue_field().gen()
            sage: F = v.residue_ring()(u); F
            u2
            sage: f = v.lift(F); f
            (2^-1 + O(2^9))*x^2 + (2^-1 + O(2^9))*x + u*2^-1 + O(2^9)
            sage: F == v.reduce(f)
            True

        """
        if F.parent() is self.residue_ring().base():
            F = self.residue_ring()(F)
        if F.parent() is not self.residue_ring():
            raise ValueError("F must be an element of the residue ring of the valuation")
        if not self.domain().base_ring().is_field():
            raise NotImplementedError("only implemented for polynomial rings over fields")
        if F.is_constant():
            if F.is_zero():
                return self.domain().zero()
            if F.is_one():
                return self.domain().one()

        if self.tau().is_zero():
            if not self._mu > 0:
                raise NotImplementedError
            if not F.is_constant():
                raise ValueError("any reduction is constant in this valuation")
            if self.phi() != self.domain().gen():
                raise NotImplementedError
            constant = F[0].lift()
            assert constant.is_constant()
            return self.domain()(constant[0])

        R0 = self._base_valuation.residue_ring()

        # in the last step of reduce, the f_iQ^i are reduced, and evaluated at
        # the generator of the residue field
        # here, we undo this:
        coeffs = [ R0(c if self.psi().degree()==1 else list(c._vector_())) for c in F.coefficients(sparse=False) ]
        coeffs = [ self._base_valuation.lift(c) for c in coeffs ]
        # now the coefficients correspond to the expansion with (f_iQ^i)(Q^{-1} phi)^i

        # now we undo the factors of Q^i (the if else is necessary to handle the case when mu is infinity, i.e., when Q_() is undefined)
        coeffs = [ (c if i == 0 else c*self.Q_()**i).map_coefficients(lambda d:_lift_to_maximal_precision(d)) for i,c in enumerate(coeffs) ]

        RR = self.domain().change_ring(self.domain())

        if self._mu is infinity:
            assert len(coeffs) <= 1
            ret = RR(coeffs)[0]
        else:
            ret = RR(coeffs)(self.phi()**self.tau())
        ret = ret.map_coefficients(lambda c:_lift_to_maximal_precision(c))
        return ret

    def lift_to_key(self, F):
        """
        Lift the irreducible polynomial ``F`` to a key polynomial.

        INPUT:

        - ``F`` -- an irreducible non-constant polynomial in the
          :meth:`residue_ring` of this valuation

        OUTPUT:

        A polynomial `f` in the domain of this valuation which is a key
        polynomial for this valuation and which, for a suitable equivalence
        unit `R`, satifies that the reduction of `Rf` is ``F``

        ALGORITHM:

        We follow the algorithm described in Theorem 13.1 [ML1936] which, after
        a :meth:`lift` of ``F``, essentially shifts the valuations of all terms
        in the `\phi`-adic expansion up and then kills the leading coefficient.

        EXAMPLES::

            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)

            sage: w = v.extension(x^2 + x + u, 1/2)
            sage: y = w.residue_ring().gen()
            sage: f = w.lift_to_key(y + 1); f
            (1 + O(2^10))*x^4 + (2 + O(2^11))*x^3 + (1 + u*2 + O(2^10))*x^2 + (u*2 + O(2^11))*x + (u + 1) + u*2 + u*2^2 + u*2^3 + u*2^4 + u*2^5 + u*2^6 + u*2^7 + u*2^8 + u*2^9 + O(2^10)
            sage: w.is_key(f)
            True

        A more complicated example::

            sage: v = GaussValuation(S)
            sage: v = v.extension(x^2 + x + u, 1)
            sage: v = v.extension((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)

            sage: u = v.residue_field().gen()
            sage: y = v.residue_ring().gen()
            sage: f = v.lift_to_key(y^3+y+u)
            sage: f.degree()
            12
            sage: v.is_key(f)
            True

        """
        if F.parent() is not self.residue_ring():
            raise ValueError("F must be an element of the residue ring of the valuation")
        if not self.domain().base_ring().is_field():
            raise NotImplementedError("only implemented for polynomial rings over fields")
        if F.is_constant():
            raise ValueError("F must not be constant")
        if not F.is_monic():
            raise ValueError("F must be monic")
        if not F.is_irreducible():
            raise ValueError("F must be irreducible")
        if F == F.parent().gen():
            return self.phi()

        f = self.lift(F)
        assert self(f) == 0
        assert self.reduce(f) == F

        f *= self.Q()**F.degree()
        CV = zip(self.coefficients(f), self.valuations(f))
        vf = self(f)
        CV = [(c,v) if v==vf else (c.parent().zero(),infinity) for c,v in CV]
        while CV[-1][1] is infinity:
            CV.pop()

        CV[-1] = (CV[-1][0].parent().one(), vf)
        ret = self.domain().change_ring(self.domain())([c for c,v in CV])(self.phi())
        ret = ret.map_coefficients(lambda c:_lift_to_maximal_precision(c))
        assert (ret == self.phi()) == (F == F.parent().gen())
        assert self.is_key(ret)
        return ret

    def tau(self):
        """
        Return the factor which is needed to turn an element of the value group
        of this valuation into an element of the value group of its base
        valuation.

        OUTPUT:

        a positive integer

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.extension(x^2 + x + u, 1/2)
            sage: w.tau()
            2
            sage: w = w.extension((x^2 + x + u)^2 + 2, 5/3)
            sage: w.tau()
            3

        """
        from sage.rings.all import ZZ

        if self._base_valuation.value_group() == 0:
            return ZZ.zero()

        assert self.value_group().numerator() == 1
        assert self._base_valuation.value_group().numerator() == 1
        return ZZ(self.value_group().denominator().gen(0)) // ZZ(self._base_valuation.value_group().denominator().gen(0))

    @cached_method
    def psi(self):
        """
        Return the minimal polynomial of the residue field extension of this valuation.

        OUTPUT:

        a polynomial in the residue ring of the base valuation

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.extension(x^2 + x + u, 1/2)
            sage: w.psi()
            x^2 + x + u0
            sage: w = w.extension((x^2 + x + u)^2 + 2, 5/3)
            sage: w.psi()
            x + 1

        """
        R = self._base_valuation.equivalence_unit(-self._base_valuation(self._phi))
        F = self._base_valuation.reduce(self._phi*R)
        assert F.is_irreducible()
        return F

    @cached_method
    def residue_field(self, generator=None):
        """
        Return the residue field of this valuation.

        INPUT:

        - ``generator`` -- a string or ``None`` (default: ``None``), the name
          of the generator of the residue field extension, if ``None`` a name
          is automatically chosen

        OUTPUT:

        a relative extension of the residue field of the base valuation

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.extension(x^2 + x + u, 1/2)
            sage: w.residue_field()
            Univariate Quotient Polynomial Ring in u1 over Finite Field in u0 of size 2^2 with modulus u1^2 + u1 + u0
            sage: w = w.extension((x^2 + x + u)^2 + 2, 5/3)
            sage: w.residue_field()
            Univariate Quotient Polynomial Ring in u1 over Finite Field in u0 of size 2^2 with modulus u1^2 + u1 + u0

        """
        if generator is None:
            level = 0
            v = self
            while isinstance(v, AugmentedValuation):
                level += 1
                v = v._base_valuation
            generator = 'u%s'%level
        if self.psi().degree() == 1:
            ret = self._base_valuation.residue_field()
            #ret.gen = "use augmented_valuation.residue_field_generator() instead" # temporary hack to find bugs when .gen() is used - .residue_field_generator() instead
            return ret
        else:
            return self._base_valuation.residue_field().extension(self.psi(), names=generator)

    @cached_method
    def residue_field_generator(self):
        if self.psi().degree() == 1:
            ret = -self.psi()[0]
        else:
            ret = self.residue_field().gen()

        assert ret.parent() is self.residue_field()
        assert self.psi()(ret).is_zero()
        return ret

    def is_commensurable_inductive(self):
        """
        Return whether this valuation is a commensurable inductive valuation
        over the discrete valuation of the base ring of the polynomial ring, as
        defined in section 4 of [ML1936].

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.extension(x^2 + x + u, 1)
            sage: w.is_commensurable_inductive()
            True

        REFERENCES:

        .. [ML1936] Mac Lane, S. (1936). A construction for prime ideals as absolute
        values of an algebraic field. Duke Mathematical Journal, 2(3), 492-510.

        """
        return self._base_valuation.is_commensurable_inductive()

    def E(self):
        """
        Return the ramification index of this valuation over its underlying
        Gauss valuation.

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.extension(x^2 + x + u, 1)
            sage: w.E()
            1
            sage: w = v.extension(x, 1/2)
            sage: w.E()
            2

        """
        return self.tau() * self._base_valuation.E()

    def F(self):
        """
        Return the degree of the residue field extension of this valuation
        over the Gauss valuation.

        EXAMPLES::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: w = v.extension(x^2 + x + u, 1)
            sage: w.F()
            2
            sage: w = v.extension(x, 1/2)
            sage: w.F()
            1

        """
        return self.psi().degree() * self._base_valuation.F()

    def change_ring(self, base_ring):
        return AugmentedValuation(self._base_valuation.change_ring(base_ring), self.phi().change_ring(base_ring), self._mu)

    def uniformizer(self):
        return self.element_with_valuation(1/self.E())
