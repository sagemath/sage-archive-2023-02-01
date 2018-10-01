"""
Tate algebra element

A class for series in Tate algebras.

AUTHOR:

- Xavier Caruso, Thibaut Verron (2018-09)

"""

# ***************************************************************************
#    Copyright (C) 2018 Xavier Caruso <xavier.caruso@normalesup.org>
#                       Thibaut Verron <thibaut.verron@gmail.com>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ***************************************************************************

from sage.structure.element import coerce_binop
from sage.structure.richcmp import op_EQ, op_NE, op_LT, op_GT, op_LE, op_GE

from sage.structure.element cimport Element
from sage.structure.element cimport MonoidElement
from sage.structure.element cimport CommutativeAlgebraElement
from sage.structure.sequence import Sequence

from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage.rings.polynomial.polydict cimport PolyDict
from sage.rings.polynomial.polydict cimport ETuple
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement



def _pushout_family(elements, initial=ZZ):
    """
    Return a parent in which ``initial`` and all elements 
    in ``elements`` coerce.

    INPUT:

    - ``elements`` -- a list of elements

    - ``initial`` -- a parent

    EXAMPLES:

        sage: from sage.rings.tate_algebra_element import _pushout_family

        sage: R = Zp(2)
        sage: A.<x,y> = TateAlgebra(R)
        sage: S.<a> = Zq(4)

        sage: _pushout_family([a, x, 3])
        Tate Algebra in x (val >= 0), y (val >= 0) over 2-adic Unramified Extension Field in a defined by x^2 + x + 1

    """
    from sage.structure.coerce_exceptions import CoercionException
    from sage.categories.pushout import pushout
    A = initial
    for elt in elements:
        if not isinstance(elt, Element):
            raise TypeError("cannot coerce all the elements to the same parent")
        try:
            A = pushout(A, elt.parent())
        except CoercionException:
            raise TypeError("cannot coerce all the elements to the same parent")
    return A


cdef class TateAlgebraTerm(MonoidElement):
    r"""
    A class for Tate algebra terms.

    A term in `K\{X_1,\dots,X_n\}` is the product of a coefficient in `K` and a
    monomial in the variables `X_1,\dots,X_n`.

    Those terms form a partially ordered monoid, with term multiplication and the
    term order of the parent Tate algebra.

    INPUT:

    - ``coeff`` -- an element in the base field

    - ``exponent`` - a tuple of length ``n``
    
    EXAMPLES::

        sage: R = Zp(2, print_mode='digits', prec=10)
        sage: A.<x,y> = TateAlgebra(R)
        sage: T = A.monoid_of_terms(); T
        Monoid of terms in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10

        sage: T(2*x*y)
        (...00000000010)*x*y
        sage: T(0)
        Traceback (most recent call last):
        ...
        ValueError: a term cannot be zero

    """
    def __init__(self, parent, coeff, exponent=None):
        """
        Initialize a Tate algebra term

        INPUT:

        - ``coeff`` -- an element in the base field

        - ``exponent`` -- a tuple

        TESTS::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()

            sage: t = T(x)
            sage: TestSuite(t).run()

        """
        MonoidElement.__init__(self, parent)
        field = parent.base_ring().fraction_field()
        if isinstance(coeff, TateAlgebraElement):
            coeff = coeff.leading_term()
        if isinstance(coeff, TateAlgebraTerm):
            if coeff.parent().variable_names() != self._parent.variable_names():
                raise ValueError("the variable names do not match")
            self._coeff = field((<TateAlgebraTerm>coeff)._coeff)
            self._exponent = (<TateAlgebraTerm>coeff)._exponent
        else:
            if coeff.is_zero():
                raise ValueError("a term cannot be zero")
            self._coeff = field(coeff)
            self._exponent = ETuple([0] * parent.ngens())
        if exponent is not None:
            if not isinstance(exponent, ETuple):
                exponent = ETuple(exponent)
            self._exponent = self._exponent.eadd(exponent)
        if len(self._exponent) != parent.ngens():
            raise ValueError("the length of the exponent does not match the number of variables")
        for i in self._exponent.nonzero_positions():
            if self._exponent[i] < 0:
                raise ValueError("only nonnegative exponents are allowed")
        if not parent.base_ring().is_field() and self.valuation() < 0:
            raise ValueError("this term is not in the ring of integers")

    cdef TateAlgebraTerm _new_c(self):
        r"""
        Fast creation of a Tate algebra term.

        TESTS::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R);
            sage: T = A.monoid_of_terms(); T
            Monoid of terms in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: t = T(x*y); t
            (...0000000001)*x*y
            sage: 2*t  # indirect doctest
            (...00000000010)*x*y

        """
        cdef TateAlgebraTerm ans = TateAlgebraTerm.__new__(TateAlgebraTerm)
        ans._parent = self._parent
        return ans

    def __reduce__(self):
        """
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: A.<x,y> = TateAlgebra(Zp(2))
            sage: t = x.leading_term()
            sage: loads(dumps(t)) == t  # indirect doctest
            True
        """
        return TateAlgebraTerm, (self.parent(), self._coeff, self._exponent)

    def _repr_(self):
        r"""
        Return a string representation of this Tate algebra term.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: T(2*x*y)  # indirect doctest
            (...00000000010)*x*y
        
        """
        parent = self._parent
        s = "(%s)" % self._coeff
        for i in range(parent._ngens):
            if self._exponent[i] == 1:
                s += "*%s" % parent._names[i]
            elif self._exponent[i] > 1:
                s += "*%s^%s" % (parent._names[i], self._exponent[i])
        return s

    def _latex_(self):
        r"""
        Return a LaTeX representation of this Tate algebra term.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: T(2*x*y)
            (...00000000010)*x*y
            sage: T(2*x*y)._latex_()
            '(...00000000010)xy'

        """
        from sage.misc.latex import latex
        parent = self._parent
        s = "(%s)" % latex(self._coeff)
        for i in range(parent._ngens):
            if self._exponent[i] == 1:
                s += "%s" % parent._latex_names[i]
            elif self._exponent[i] > 1:
                s += "%s^{%s}" % (parent._latex_names[i], self._exponent[i])
        return s

    def coefficient(self):
        r"""
        Return the coefficient of this Tate algebra term.

        EXAMPLES::

            sage: R = Zp(2,prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: t = T(2*x*y); t
            (2 + O(2^11))*x*y
            sage: t.coefficient() 
            2 + O(2^11)

        """
        return self._coeff

    def exponent(self):
        r"""
        Return the exponents of this Tate algebra term.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms(); T
            Monoid of terms in x (val >= 0), y (val >= 0) over 2-adic Field with capped relative precision 10
            sage: t = T(2,(1,1))
            sage: t.exponent()
            (1, 1)

        """
        return self._exponent

    cpdef _mul_(self, other):
        r"""
        Return the product of this Tate algebra term with ``other``.

        INPUT:

        - ``other`` - a Tate algebra term

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(2*x*y); s
            (...00000000010)*x*y
            sage: t = T(3*x^2*y); t
            (...0000000011)*x^2*y
            sage: s*t  # indirect doctest
            (...00000000110)*x^3*y^2

        """
        cdef TateAlgebraTerm ans = self._new_c()
        ans._exponent = self._exponent.eadd((<TateAlgebraTerm>other)._exponent)
        ans._coeff = self._coeff * (<TateAlgebraTerm>other)._coeff
        return ans

    cdef int _cmp_c(self, TateAlgebraTerm other):
        r"""
        Compare the Tate algebra term with ``other``.

        INPUT:

        - ``other`` -- a term

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(x^2*y^2)
            sage: t = T(x*y^3)
            sage: s > t  # indirect doctest
            True

        """
        cdef int c = other._valuation_c() - self._valuation_c()
        if not c:
            skey = self._parent.term_order().sortkey
            ks = skey(self._exponent)
            ko = skey(other._exponent)
            c = (ks > ko) - (ks < ko)
        return c

    cpdef _richcmp_(self, other, int op):
        r"""
        Compare the Tate algebra term with ``other`` according to
        the rich comparison operator ``op``.

        The term `a*X^A` is smaller or equal to `b*X^B` if its valuation is
        greater than the valuation of `b*X^B`, or when equality occurs, when 
        `A` is smaller or equal to `B` for the Tate algebra monomial order.
 
        INPUT:

        - ``other`` -- a Tate algebra term

        - ``op`` -- the comparison operator

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: T.term_order()
            Degree reverse lexicographic term order
            sage: s = T(x^2*y^2); s
            (...0000000001)*x^2*y^2
            sage: t = T(x*y^3); t
            (...0000000001)*x*y^3
            sage: s < t  # indirect doctest
            False
            sage: s > t  # indirect doctest
            True
            sage: s <= t  # indirect doctest
            False
            sage: s >= t  # indirect doctest
            True

        Elements with the same valuation and monomial are equivalent
        for the preorder::
        
            sage: ss = T(3*x^2*y^2); ss
            (...0000000011)*x^2*y^2
            sage: s < ss  # indirect doctest
            False
            sage: s > ss  # indirect doctest
            False
            sage: s <= ss  # indirect doctest
            True
            sage: s >= ss  # indirect doctest
            True
            sage: s == ss
            False
        
        """
        if op == op_EQ:
            return ((<TateAlgebraTerm>self)._coeff == (<TateAlgebraTerm>other)._coeff 
                and (<TateAlgebraTerm>self)._exponent == (<TateAlgebraTerm>other)._exponent)
        if op == op_NE:
            return ((<TateAlgebraTerm>self)._coeff != (<TateAlgebraTerm>other)._coeff 
                 or (<TateAlgebraTerm>self)._exponent != (<TateAlgebraTerm>other)._exponent)
        c = (<TateAlgebraTerm>self)._cmp_c(<TateAlgebraTerm>other)
        if op == op_LT:
            return c < 0
        if op == op_LE:
            return c <= 0
        if op == op_GT:
            return c > 0
        if op == op_GE:
            return c >= 0

    cpdef TateAlgebraTerm monomial(self):
        r"""
        Return this term divided by its coefficient.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(3*x^2*y^2); s
            (...0000000011)*x^2*y^2
            sage: s.monomial()
            (...0000000001)*x^2*y^2

        """
        cdef TateAlgebraTerm ans = self._new_c()
        ans._coeff = self._parent._field.fraction_field()(1)
        ans._exponent = self._exponent
        return ans

    cpdef TateAlgebraTerm monic(self):
        r"""
        Return this term normalized so that it has valuation 0 
        and its coefficient is a power of the uniformizer.

        EXAMPLES:

        When the log radii of convergence are all zero, the
        coefficient of the returned term is `1`. In this case,
        this method does the same thing as :meth:`monomial`::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(3*x^2*y^2); s
            (...0000000011)*x^2*y^2
            sage: s.monic()
            (...0000000001)*x^2*y^2
            sage: s.monomial()
            (...0000000001)*x^2*y^2

        However, when log radii do not vanish, behaviors might
        be different::

            sage: A.<x,y> = TateAlgebra(R, log_radii=1)
            sage: T = A.monoid_of_terms()
            sage: s = T(3*x^2*y^2); s
            (...0000000011)*x^2*y^2
            sage: s.monic()
            (...00000000010000)*x^2*y^2
            sage: s.monomial()
            (...0000000001)*x^2*y^2

        We compare the valuations::

            sage: s.monic().valuation()
            0
            sage: s.monomial().valuation()
            -4

        """
        cdef TateAlgebraTerm ans = self._new_c()
        cdef long v = self._exponent.dotprod(self._parent._log_radii)
        ans._coeff = self._parent._field.uniformizer_pow(v)
        ans._exponent = self._exponent
        return ans

    def valuation(self):
        r"""
        Return the valuation of this term.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: t = T(4*x^2*y^2); t
            (...000000000100)*x^2*y^2
            sage: t.valuation()
            2

        In case of nonzero log radii, the valuations of the variables
        contribute::

            sage: A.<x,y> = TateAlgebra(R, log_radii=1)
            sage: T = A.monoid_of_terms()
            sage: t = T(4*x^2*y^2); t
            (...000000000100)*x^2*y^2
            sage: t.valuation()
            -2

        """
        return ZZ(self._valuation_c())

    cdef long _valuation_c(self):
        r"""
        Return the valuation of this term.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: t = T(4*x^2*y^2); t
            (...000000000100)*x^2*y^2
            sage: t.valuation()  # indirect doctest
            2
        
        """
        return (<pAdicGenericElement>self._coeff).valuation_c() - <long>self._exponent.dotprod(self._parent._log_radii)

    cdef Element _call_c(self, list arg):
        """
        Return this term evaluated at ``args``.

        INPUT:

        - ``args`` -- elements

        TESTS::

            sage: R = Zp(2)
            sage: A.<x,y> = TateAlgebra(R)
            sage: M = A.monoid_of_terms()
            sage: t = M(x*y)
            sage: t(0, 1)
            0

        """
        cdef Element ans = self._coeff
        cdef ETuple exponent = self._exponent
        cdef size_t ind
        for ind in range(exponent._nonzero):
            ans *= arg[exponent._data[2*ind]] ** exponent._data[2*ind+1]
        return ans

    def __call__(self, *args):
        """
        Return this term evaluated at ``args``.

        INPUT:

        - ``args`` -- elements

        EXAMPLES::

            sage: R = Zp(2)
            sage: A.<x,y> = TateAlgebra(R)
            sage: M = A.monoid_of_terms()
            sage: t = M(x*y)
            sage: t(1, 2)
            2 + O(2^21)

        An error is raised if we ask for the evaluation at one
        point which is outside the domain of convergence::

            sage: t(1/2, 1)
            Traceback (most recent call last):
            ...
            ValueError: the given values are outside the domain of convergence

        TESTS::

            sage: t(1/2, 1, 0)
            Traceback (most recent call last):
            ...
            TypeError: wrong number of arguments

            sage: t(1/2, GF(3)(2))
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce all the elements to the same parent

        """
        parent = self._parent
        if len(args) != parent._ngens:
            raise TypeError("wrong number of arguments")
        A = _pushout_family(args, parent._field)
        args = [ A(arg) for arg in args ]
        ratio = A.absolute_e() // parent._base.absolute_e()
        for i in range(parent._ngens):
            if args[i].valuation() < -ratio * parent._log_radii[i]:
                raise ValueError("the given values are outside the domain of convergence")
        res = self._call_c(args)
        if parent._integral:
            try:
                res = res.parent().integer_ring()(res)
            except AttributeError:
                pass
        return res

    @coerce_binop
    def is_coprime_with(self, other):
        r"""
        Return ``True`` if this term is coprime with ``other``.

        INPUT:

        - ``other`` - a Tate term

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: t = T(4*x^2*y^2); t
            (...000000000100)*x^2*y^2
            sage: s = T(y^3); s
            (...0000000001)*y^3
            sage: s.is_coprime_with(t)
            False
            sage: t.is_coprime_with(s)
            False

            sage: tt = T(3*x^2); tt
            (...0000000011)*x^2
            sage: s.is_coprime_with(tt)
            True
            sage: tt.is_coprime_with(s)
            True

        When working over a rational Tate algebra, only the 
        monomial part of terms are compared::
        
            sage: t = T(2*x^2); t
            (...00000000010)*x^2
            sage: s = T(4*y^3); s
            (...000000000100)*y^3
            sage: s.is_coprime_with(t)
            True

        But coefficients play a role when we are working over
        the ring of integers of the Tate Algebra::

            sage: AA = A.integer_ring()
            sage: TT = AA.monoid_of_terms()
            sage: TT(s).is_coprime_with(TT(t))
            False

        """
        for i in range(self._parent.ngens()):
            if self._exponent[i] > 0 and other.exponent()[i] > 0:
                return False
        if self._parent.base_ring().is_field():
            return True
        else:
            return self.valuation() == 0 or other.valuation() == 0

    @coerce_binop
    def gcd(self, other):
        r"""
        Return the greatest common divisor of this term and ``other``.

        The result is normalized so that its valuation is equal to
        the smallest valuation of this term and ``other``.

        INPUT:

        - ``other`` - a Tate term

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(8*x^2*y^2); s
            (...0000000001000)*x^2*y^2
            sage: t = T(4*x*y^3); t
            (...000000000100)*x*y^3
            sage: s.gcd(t)
            (...000000000100)*x*y^2

        """
        return self._gcd_c(other)

    cdef TateAlgebraTerm _gcd_c(self, TateAlgebraTerm other):
        r"""
        Return the greatest common divisor of this term and ``other``.

        The result is normalized so that its valuation is equal to
        the smallest valuation of this term and ``other``.

        INPUT:

        - ``other`` - a Tate term

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(8*x^2*y^2); s
            (...0000000001000)*x^2*y^2
            sage: t = T(4*x*y^3); t
            (...000000000100)*x*y^3
            sage: s.gcd(t)  # indirect doctest
            (...000000000100)*x*y^2

        """
        cdef TateAlgebraTerm ans = self._new_c()
        ans._exponent = self._exponent.emin(other._exponent)
        if self._coeff.valuation() < other._coeff.valuation():
            ans._coeff = self._coeff
        else:
            ans._coeff = other._coeff
        return ans

    @coerce_binop
    def lcm(self, other):
        r"""
        Return the least common multiple of two Tate terms.

        The result is normalized so that its valuation is equal to
        the largest valuation of this term and ``other``.
        
        INPUT:

        - ``other`` - a Tate term

        EXAMPLES::

        In a Tate algebra over a field:

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(8*x^2*y^2); s
            (...0000000001000)*x^2*y^2
            sage: t = T(4*x*y^3); t
            (...000000000100)*x*y^3
            sage: s.lcm(t)
            (...0000000001000)*x^2*y^3
        
        """        
        return self._lcm_c(other)

    cdef TateAlgebraTerm _lcm_c(self, TateAlgebraTerm other):
        r"""
        Return the least common multiple of two Tate terms.

        The result is normalized so that its valuation is equal to
        the largest valuation of this term and ``other``.
        
        INPUT:

        - ``other`` - a Tate term

        EXAMPLES::

        In a Tate algebra over a field:

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(8*x^2*y^2); s
            (...0000000001000)*x^2*y^2
            sage: t = T(4*x*y^3); t
            (...000000000100)*x*y^3
            sage: s.lcm(t)  # indirect doctest
            (...0000000001000)*x^2*y^3
        
        """        
        cdef TateAlgebraTerm ans = self._new_c()
        ans._exponent = self._exponent.emax(other._exponent)
        if self._coeff.valuation() < other._coeff.valuation():
            ans._coeff = other._coeff
        else:
            ans._coeff = self._coeff
        return ans

    @coerce_binop
    def is_divisible_by(self, other, integral=False):
        r"""
        Return ``True`` if this term is divisible by ``other``.

        INPUT:

        - ``other`` - a Tate term
        
        - ``integral`` - (default: ``False``); if ``True``, test 
          for divisibility in the ring of integers of the Tate algebra

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(4*x^2*y^2); s
            (...000000000100)*x^2*y^2
            sage: t = T(4*x*y^3); t
            (...000000000100)*x*y^3
            sage: s.is_divisible_by(t)
            False

            sage: t = T(4*x*y^2); t
            (...000000000100)*x*y^2
            sage: s.is_divisible_by(t)
            True

            sage: t = T(16); t
            (...00000000010000)
            sage: s.is_divisible_by(t)
            True
            sage: s.is_divisible_by(t, integral=True)
            False

        If you are working over the ring of integers of the Tate algebra,
        divisibility is always checked in the ring of integers (even if 
        ``integral`` is set to ``False``)::

            sage: AA = A.integer_ring()
            sage: TT = AA.monoid_of_terms()
            sage: ss = TT(s)
            sage: tt = TT(t)
            sage: ss.is_divisible_by(tt)
            False
            sage: ss.is_divisible_by(tt, integral=False)
            False
            
        Be careful that coercion between the Tate algebra and its ring of 
        integers can be done silently::

            sage: s.is_divisible_by(tt)
            True

        """
        return (<TateAlgebraTerm?>other)._divides_c(self, integral)
    
    @coerce_binop
    def divides(self, other, integral=False):
        r"""
        Return ``True`` if this term divides ``other``.

        INPUT:

        - ``other`` - a Tate term
        
        - ``integral`` - (default: ``False``); if ``True``, test for 
          divisibility in the ring of integers of the Tate algebra

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(4*x^2*y^2); s
            (...000000000100)*x^2*y^2
            sage: t = T(4*x*y^3); t
            (...000000000100)*x*y^3
            sage: t.divides(s)
            False

            sage: t = T(4*x*y^2); t
            (...000000000100)*x*y^2
            sage: t.divides(s)
            True

            sage: t = T(16); t
            (...00000000010000)
            sage: t.divides(s)
            True
            sage: t.divides(s, integral=True)
            False

        If you are working over the ring of integers of the Tate algebra,
        divisibility is always checked in the ring of integers (even if 
        ``integral`` is set to ``False``)::

            sage: AA = A.integer_ring()
            sage: TT = AA.monoid_of_terms()
            sage: ss = TT(s)
            sage: tt = TT(t)
            sage: tt.divides(ss)
            False
            sage: tt.divides(ss, integral=False)
            False
            
        Be careful that coercion between the Tate algebra and its ring of 
        integers can be done silently::

            sage: tt.divides(s)
            True

        """
        return self._divides_c(other, integral)

    cdef bint _divides_c(self, TateAlgebraTerm other, bint integral):
        r"""
        Return ``True`` if this term divides ``other``.

        INPUT:

        - ``other`` - a Tate term
        
        - ``integral`` - (default: ``False``) if ``True``, test for 
          divisibility in the ring of integers of the Tate algebra

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: s = T(4*x^2*y^2); s
            (...000000000100)*x^2*y^2
            sage: t = T(4*x*y^3); t
            (...000000000100)*x*y^3
            sage: t.divides(s)  # indirect doctest
            False

        """
        parent = self._parent
        if (integral or not parent.base_ring().is_field()) and self.valuation() > other.valuation():
            return False
        for i in range(parent._ngens):
            if self._exponent[i] > other._exponent[i]:
                return False
        return True

    @coerce_binop
    def __floordiv__(self, other):
        r"""
        Return the result of the exact division of this term by ``other``.

        INPUT:

        - ``other`` - a Tate term

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: t = T(2*x^2*y^3); t
            (...00000000010)*x^2*y^3
            sage: s = T(6*x*y^2); s
            (...00000000110)*x*y^2
            sage: t // s
            (...1010101011)*x*y

        If the Tate terms are not divisible, an error is raised::

            sage: s // t
            Traceback (most recent call last):
            ...
            ValueError: the division is not exact

        """
        if not self.is_divisible_by(other):
            raise ValueError("the division is not exact")
        return (<TateAlgebraTerm>self)._floordiv_c(<TateAlgebraTerm>other)

    cdef TateAlgebraTerm _floordiv_c(self, TateAlgebraTerm other):
        r"""
        Return the result of the exact division of this term by ``other``.

        INPUT:

        - ``other`` - a Tate term

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: T = A.monoid_of_terms()
            sage: t = T(2*x^2*y^3); t
            (...00000000010)*x^2*y^3
            sage: s = T(6*x*y^2); s
            (...00000000110)*x*y^2
            sage: t // s
            (...1010101011)*x*y

        If the Tate terms are not divisible, an error is raised::

            sage: s // t
            Traceback (most recent call last):
            ...
            ValueError: the division is not exact

        """
        cdef TateAlgebraTerm ans = self._new_c()
        ans._exponent = self.exponent().esub(other.exponent())
        ans._coeff = self._coeff / other._coeff
        return ans


cdef class TateAlgebraElement(CommutativeAlgebraElement):
    r"""
    A class for Tate series, elements of Tate algebras.

    EXAMPLES::

        sage: R = Zp(2,prec=10,print_mode='digits')
        sage: A.<x,y> = TateAlgebra(R)
        sage: A(2*x+1)
        (...0000000001) + (...00000000010)*x
        sage: A(2*x+1, prec=5)
        (...00001) + (...00010)*x + O(2^5)
        sage: A(2*x+1, prec=20)
        (...0000000001) + (...00000000010)*x + O(2^20)

    """
    def __init__(self, parent, x, prec=None, reduce=True):
        r"""
        Initialise a Tate algebra series.

        TESTS::

            sage: R = Zp(2,prec=10,print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: TestSuite(x).run()

        """
        cdef TateAlgebraElement xc
        CommutativeAlgebraElement.__init__(self, parent)
        self._prec = Infinity
        if isinstance(x, TateAlgebraElement):
            xc = <TateAlgebraElement>x
            xparent = x.parent()
            if xparent is parent:
                self._poly = PolyDict(xc._poly.__repn)
                self._prec = xc._prec
            elif xparent.variable_names() == parent.variable_names():
                ratio = parent._base.absolute_e() / xparent.base_ring().absolute_e()
                for i in range(parent.ngens()):
                    if parent.log_radii()[i] > xparent.log_radii()[i] * ratio:
                        raise ValueError("cannot restrict to a bigger domain")
                self._poly = PolyDict({ e: parent._field(v) for (e,v) in xc._poly.__repn.iteritems() })
                if xc._prec is not Infinity:
                    self._prec = (xc._prec * ratio).ceil()
            else:
                raise TypeError("variable names do not match")
        elif isinstance(x, TateAlgebraTerm):
            xparent = x.parent()
            if xparent.variable_names() == parent.variable_names():
                ratio = parent._base.absolute_e() / xparent.base_ring().absolute_e()
                for i in range(parent.ngens()):
                    if parent.log_radii()[i] > xparent.log_radii()[i] * ratio:
                        raise ValueError("cannot restrict to a bigger domain")
                self._poly = PolyDict({(<TateAlgebraTerm>x)._exponent: parent._field((<TateAlgebraTerm>x)._coeff)})
        else:
            poly = parent._polynomial_ring(x)
            self._poly = PolyDict(poly.dict())
        if prec is not None:
            self._prec = min(self._prec, prec)
        self._normalize()
        self._terms = None
        if not parent.base_ring().is_field() and self.valuation() < 0:
            raise ValueError("this series is not in the ring of integers")

    cdef TateAlgebraElement _new_c(self):
        """
        Fast creation of a new Tate series.

        EXAMPLES::

            sage: R = Zp(2,prec=10,print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: x + y  # indirect doctest
            (...0000000001)*x + (...0000000001)*y

        """
        cdef TateAlgebraElement ans = TateAlgebraElement.__new__(TateAlgebraElement)
        ans._parent = self._parent
        ans._is_normalized = False
        return ans

    cdef _normalize(self):
        """
        Normalize this series.

            sage: R = Zp(2,prec=10,print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: A(78612, prec=3)  # indirect doctest
            (...100) + O(2^3)

        """
        self._is_normalized = True
        if self._prec is Infinity: return
        cdef int v
        for (e,c) in self._poly.__repn.iteritems():
            v = (<ETuple>self._parent._log_radii).dotprod(<ETuple>e)
            self._poly.__repn[e] = self._poly.__repn[e].add_bigoh((self._prec - v).ceil())

    def __reduce__(self):
        """
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: A.<x,y> = TateAlgebra(Zp(2))
            sage: loads(dumps(x)) == x  # indirect doctest
            True
        """
        return TateAlgebraElement, (self.parent(), self._poly, self._prec)

    def _repr_(self):
        r"""
        Return a string representation of this series.

        The terms are ordered with decreasing term order 
        (increasing valuation of the coefficients, then 
        the monomial order of the parent algebra).

        EXAMPLES::
        
            sage: R = Zp(2, 10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: x + 2*x^2 + x^3
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2

            sage: A(x + 2*x^2 + x^3, prec=5)
            (...00001)*x^3 + (...00001)*x + (...00010)*x^2 + O(2^5)

        """        
        base = self._parent.base_ring()
        nvars = self._parent.ngens()
        vars = self._parent.variable_names()
        s = ""
        for t in self.terms():
            if t.valuation() >= self._prec:
                continue
            s += " + %s" % t
        if self._prec is not Infinity:
            s += " + %s" % base.change(show_prec="bigoh")(0, self._prec)
        if s == "":
            return "0"
        return s[3:]

    def _latex_(self):
        r"""
        Return a LaTeX representation of this series.

        The terms are ordered with decreasing term order 
        (increasing valuation of the coefficients, then 
        the monomial order of the parent algebra).

        EXAMPLES::
        
            sage: R = Zp(2, 10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: f._latex_()
            '(...0000000001)x^{3} + (...0000000001)x + (...00000000010)x^{2}'

        """
        from sage.misc.latex import latex
        base = self._parent.base_ring()
        nvars = self._parent.ngens()
        vars = self._parent.variable_names()
        s = ""
        for t in self.terms():
            if t.valuation() >= self._prec:
                continue
            s += " + %s" % latex(t)
        if self._prec is not Infinity:
            s += " + %s" % latex(base.change(show_prec="bigoh")(0, self._prec))
        if s == "":
            return "0"
        return s[3:]

    cpdef _add_(self, other):
        r"""
        Return the sum of this series and ``other``.

        The precision of the output is adjusted to the minimum of both precisions.

        INPUT:

        - ``other`` - a Tate series

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: g = A(x + 2*y^2, prec=5); g
            (...00001)*x + (...00010)*y^2 + O(2^5)
            sage: h = f + g; h  # indirect doctest
            (...00001)*x^3 + (...00010)*x^2 + (...00010)*y^2 + (...00010)*x + O(2^5)

            sage: f.precision_absolute()
            +Infinity
            sage: g.precision_absolute()
            5
            sage: h.precision_absolute()
            5

        """
        # TODO: g = x + 2*y^2 + O(2^5) fails coercing O(2^5) into R, is that normal? 
        # If we force conversion with R(O(2^5)) the result is a Tate series with precision infinity... it's understandable but confusing, no? 
        cdef TateAlgebraElement ans = self._new_c()
        ans._poly = self._poly + (<TateAlgebraElement>other)._poly
        ans._prec = min(self._prec, (<TateAlgebraElement>other)._prec)
        if self._prec != (<TateAlgebraElement>other)._prec:
            ans._normalize()
        return ans

    cpdef _neg_(self):
        r"""
        Return the opposite of this series.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: -f  # indirect doctest
            (...1111111111)*x^3 + (...1111111111)*x + (...11111111110)*x^2

        """
        cdef TateAlgebraElement ans = self._new_c()
        cdef Element s = self._parent.base_ring()(-1)
        ans._poly = self._poly.scalar_lmult(s)
        ans._prec = self._prec
        return ans

    cpdef _sub_(self, other):
        r"""
        Return the difference of this series and ``other``.

        The precision of the output is adjusted to the minimum of both precisions.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: g = A(x + 2*y^2, prec=5); g
            (...00001)*x + (...00010)*y^2 + O(2^5)
            sage: h = f - g; h # indirect doctest
            (...00001)*x^3 + (...00010)*x^2 + (...11110)*y^2 + O(2^5)

        ::

            sage: f.precision_absolute()
            +Infinity
            sage: g.precision_absolute()
            5
            sage: h.precision_absolute()
            5

        """
        cdef TateAlgebraElement ans = self._new_c()
        ans._poly = self._poly - (<TateAlgebraElement>other)._poly
        ans._prec = min(self._prec, (<TateAlgebraElement>other)._prec)
        if self._prec != (<TateAlgebraElement>other)._prec:
            ans._normalize()
        return ans

    cpdef _mul_(self, other):
        r"""
        Return the product of this series with ``other``.

        The precision is adjusted to match the best precision on the output.

        EXAMPLES::

            sage: R = Zp(2, 10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x + 4*x^2 + 2*x^3; f
            (...00000000010)*x^3 + (...00000000010)*x + (...000000000100)*x^2
            sage: g = A(x + 2*x^2, prec=5); g
            (...00001)*x + (...00010)*x^2 + O(2^5)
            sage: h = f * g; h # indirect doctest
            (...001010)*x^4 + (...000010)*x^2 + (...000100)*x^5 + (...001000)*x^3 + O(2^6)

        ::

            sage: f.precision_absolute()
            +Infinity
            sage: g.precision_absolute()
            5
            sage: h.precision_absolute()
            6

        """
        cdef TateAlgebraElement ans = self._new_c()
        a = self._prec + (<TateAlgebraElement>other).valuation()
        b = self.valuation() + (<TateAlgebraElement>other)._prec
        ans._poly = self._poly * (<TateAlgebraElement>other)._poly
        ans._prec = min(a, b)
        #ans._normalize()
        return ans

    cpdef _lmul_(self, Element right):
        r"""
        Return the product of this series by ``right``.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: 2*f # indirect doctest
            (...00000000010)*x^3 + (...00000000010)*x + (...000000000100)*x^2

            sage: 6*f # indirect doctest
            (...00000000110)*x^3 + (...00000000110)*x + (...000000001100)*x^2

        """
        cdef TateAlgebraElement ans = self._new_c()
        ans._poly = self._poly.scalar_lmult(right)
        ans._prec = self._prec + (<pAdicGenericElement>self._parent._base(right)).valuation_c()
        return ans

    cpdef _richcmp_(self, other, int op):
        r"""
        Compare this series with ``other`` according to
        the rich comparison operator ``op``.

        INPUT:

        - ``other`` -- a Tate algebra term

        - ``op`` -- the comparison operator

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: A.term_order()
            Degree reverse lexicographic term order

        For terms, we first compare the valuation and ties are broken 
        using the term ordering on `A`::

            sage: 2*x^2*y^2 > x*y^3
            False
            sage: x^2*y^2 > x*y^3
            True

        For general series, leading terms are first compared. In case
        of tie, second leading terms are compared, and so on::
        
            sage: x^2*y^2 > x*y^3 + y^4
            True
            sage: x^2*y^2 > x^2*y^2 + x*y^3
            False

        TESTS::

            sage: f = x^4 + 4*x*y + 1; f
            (...0000000001)*x^4 + (...0000000001) + (...000000000100)*x*y
            sage: g = f + 2
            sage: f == g
            False
            sage: f == g - 2
            True

        """
        diff = self - other
        c = None
        if diff.is_zero():
            c = 0
        if op == op_EQ:
            return c is not None
        if op == op_NE:
            return c is None
        if c is None:
            ts = self.terms()
            to = other.terms()
            for i in range(min(len(ts), len(to))):
                c = (<TateAlgebraTerm>ts[i])._cmp_c(<TateAlgebraTerm>to[i])
                if c: break
            else:
                c = len(ts) - len(to)
        if op == op_LT:
            return c < 0
        if op == op_LE:
            return c <= 0
        if op == op_GT:
            return c > 0
        if op == op_GE:
            return c >= 0

    def __call__(self, *args):
        """
        Return this term evaluated at ``args``.

        INPUT:

        - ``args`` -- elements

        EXAMPLES::

            sage: R = Zp(2)
            sage: A.<u,v> = TateAlgebra(R, log_radii=[0,-1])
            sage: A
            Tate Algebra in u (val >= 0), v (val >= 1) over 2-adic Field with capped relative precision 20

            sage: f = u^2 + v^2
            sage: f(1, 0)
            1 + O(2^20)

        An error is raised if we ask for the evaluation at one
        point which is outside the domain of convergence::

            sage: f(1, 1)
            Traceback (most recent call last):
            ...
            ValueError: the given values are outside the domain of convergence

        Evaluation at points in extensions is allowed::

            sage: S.<a> = Zq(2^3)
            sage: f(a, 2)
            a^2 + 2^2 + O(2^20)

            sage: T.<pi> = S.extension(x^2 - 2)
            sage: f(pi, 2)
            pi^2 + pi^4 + O(pi^42)

            sage: f(pi, pi)
            Traceback (most recent call last):
            ...
            ValueError: the given values are outside the domain of convergence

        This method can also be used to compose Tate series::

            sage: f(u + v, 2*u)
            (1 + 2^2 + O(2^20))*u^2 + (2 + O(2^20))*u*v + (1 + O(2^20))*v^2

        or for partial evaluation::

            sage: f(pi, v)
            (pi^2 + O(pi^42)) + (1 + O(pi^40))*v^2

        """
        cdef TateAlgebraTerm t
        parent = self._parent
        if len(args) != parent._ngens:
            raise TypeError("wrong number of arguments")
        A = _pushout_family(args, parent._field)
        args = [ A(arg) for arg in args ]
        ratio = A.absolute_e() // parent._base.absolute_e()
        for i in range(parent._ngens):
            if args[i].valuation() < -ratio * parent._log_radii[i]:
                raise ValueError("the given values are outside the domain of convergence")
        res = A(0, self._prec)
        for t in self._terms_c():
            if t._valuation_c() >= res.precision_absolute():
                break
            res += t._call_c(args)
        if parent._integral:
            try:
                res = res.parent().integer_ring()(res)
            except AttributeError:
                pass
        return res

    cdef TateAlgebraElement _term_mul_c(self, TateAlgebraTerm term):
        r"""
        Return the product of this series by the term ``term``.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: t = A.monoid_of_terms()(3*x^2); t
            (...0000000011)*x^2
            sage: f = x^4 + 4*x*y + 1; f
            (...0000000001)*x^4 + (...0000000001) + (...000000000100)*x*y
            sage: t*f # indirect doctest
            (...0000000011)*x^6 + (...0000000011)*x^2 + (...000000001100)*x^3*y
        
        """
        cdef TateAlgebraElement ans = self._new_c()
        ans._poly = self._poly.term_lmult(term._exponent, term._coeff)
        ans._prec = self._prec + term._valuation_c()
        return ans

    cdef TateAlgebraElement _positive_lshift_c(self, n):
        r"""
        Return the product of this series by the ``n``-th power 
        of the uniformizer.
        
        INPUT:

        - ``n`` -- a non-negative integer

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: f << 2  # indirect doctest
            (...000000000100)*x^3 + (...000000000100)*x + (...0000000001000)*x^2

        """
        cdef dict coeffs = { }
        cdef ETuple e
        cdef Element c
        cdef TateAlgebraElement ans = self._new_c()
        for (e,c) in self._poly.__repn.iteritems():
            coeffs[e] = c << n
        ans._poly = PolyDict(coeffs, self._poly.__zero)
        ans._prec = self._prec + n
        return ans

    cdef TateAlgebraElement _lshift_c(self, n):
        r"""
        Return the product of this series by the ``n``-th power 
        of the uniformizer.

        INPUT:

        - ``n`` -- an integer

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: f << 2  # indirect doctest
            (...000000000100)*x^3 + (...000000000100)*x + (...0000000001000)*x^2

        """
        cdef dict coeffs = { }
        cdef ETuple e
        cdef Element c
        cdef TateAlgebraElement ans = self._new_c()
        parent = self._parent
        base = parent.base_ring()
        if base.is_field():
            for (e,c) in self._poly.__repn.iteritems():
                coeffs[e] = c << n
            ans._prec = self._prec + n
        else:
            field = base.fraction_field()
            ngens = parent.ngens()
            for (e,c) in self._poly.__repn.iteritems():
                minval = ZZ(e.dotprod(<ETuple>parent._log_radii)).ceil()
                coeffs[e] = field(base(c) >> (minval-n)) << minval
            ans._prec = max(ZZ(0), self._prec - n)
        ans._poly = PolyDict(coeffs, self._poly.__zero)
        return ans

    def __lshift__(self, n):
        r"""
        Return the product of this series by the ``n``th power 
        of the uniformizer.

        INPUT:

        - ``n`` - an integer

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: f << 2  # indirect doctest
            (...000000000100)*x^3 + (...000000000100)*x + (...0000000001000)*x^2
            sage: f << -1  # indirect doctest
            (...000000000.1)*x^3 + (...000000000.1)*x + (...0000000001)*x^2

        If we're shifting by a negative number of digits over the ring of 
        integers of a Tate algebra, the result is truncated -- that is, the 
        output is the result of the integer division of the Tate series by 
        `\pi^{-n}` where `\pi` is a uniformizer.

            sage: AA = A.integer_ring()
            sage: AA(f) << -1
            (...0000000001)*x^2 + (...000000000)*x^3 + (...000000000)*x

        """
        return (<TateAlgebraElement>self)._lshift_c(n)

    def __rshift__(self, n):
        r"""
        Return the quotient in the division of this series by
        the ``n``-th power of the uniformizer.

        INPUT:

        - ``n`` - an integer

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: f << 2
            (...000000000100)*x^3 + (...000000000100)*x + (...0000000001000)*x^2
            sage: f << -1  # indirect doctest
            (...000000000.1)*x^3 + (...000000000.1)*x + (...0000000001)*x^2

        If we're working over the ring of integers of a Tate algebra, the 
        result is truncated -- that is, the output is the result of the integer 
        division of the Tate series by `\pi^n` where `\pi` is a uniformizer.

            sage: AA = A.integer_ring()
            sage: AA(f) << -1
            (...0000000001)*x^2 + (...000000000)*x^3 + (...000000000)*x

        """
        return (<TateAlgebraElement>self)._lshift_c(-n)

    def is_zero(self, prec=None):
        r"""
        Return ``True`` if this series is indistiguishable from zero.

        INPUT:

        - ``prec`` - an integer or ``None`` (default: ``None``), 
          the precision at which the series should be compared to zero

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2 + x^3; f
            (...0000000001)*x^3 + (...0000000001)*x + (...00000000010)*x^2
            sage: f.is_zero()
            False

            sage: g = f << 4; g
            (...00000000010000)*x^3 + (...00000000010000)*x + (...000000000100000)*x^2
            sage: g.is_zero()
            False
            sage: g.is_zero(5)
            False
            sage: g.is_zero(4)
            True
        
        """
        if prec is None:
            prec = self._prec
        return self.valuation() >= prec

    def inverse_of_unit(self):
        r"""
        Return the inverse of this series if it is invertible.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = A(1); f
            (...0000000001)
            sage: f.inverse_of_unit()
            (...0000000001) + O(2^10)

            sage: f = 2*x + 1; f
            (...0000000001) + (...00000000010)*x
            sage: f.inverse_of_unit()
            (...0000000001) + (...1111111110)*x + (...0000000100)*x^2 + (...1111111000)*x^3 + (...0000010000)*x^4 + (...1111100000)*x^5 + (...0001000000)*x^6 + (...1110000000)*x^7 + (...0100000000)*x^8 + (...1000000000)*x^9 + O(2^10)

        If the series is not invertible, an error is raised::

            sage: f = 1 + x; f
            (...0000000001)*x + (...0000000001)
            sage: f.inverse_of_unit()
            Traceback (most recent call last):
            ...        
            ValueError: this series in not invertible
        
        """
        if not self.is_unit():
            raise ValueError("this series in not invertible")
        t = self.leading_term()
        c = t.coefficient()
        parent = self._parent
        v, c = c.val_unit()
        cap = parent.precision_cap()
        inv = parent(c.inverse_of_unit()).add_bigoh(cap)
        x = (self >> v).add_bigoh(cap)
        prec = 1
        while prec < cap:
            prec *= 2
            inv = 2*inv - x*inv*inv
        return inv << v

    def is_unit(self):
        r"""
        Return ``True`` if this series in invertible.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x + 1; f
            (...0000000001) + (...00000000010)*x
            sage: f.is_unit()
            True

            sage: f = 1 + x; f
            (...0000000001)*x + (...0000000001)
            sage: f.is_unit()
            False

        Note that invertibility is tested in the parent of this series::

            sage: f = 4*x + 2
            sage: f.is_unit()
            True

            sage: AA = A.integer_ring()
            sage: AA(f).is_unit()
            False

        """
        if self.is_zero():
            return False
        t = self.leading_term()
        if max(t.exponent()) != 0:
            return False
        base = self.base_ring()
        if not base.is_field() and t.valuation() > 0:
            return False
        return True

    def restriction(self, log_radii):
        r"""
        Return the restriction of this series to a smaller domain.

        INPUT:

        - ``log_radii`` -- an integer or a tuple; the log-radii of
          convergence of the smaller domain (see :class:`TateAlgebra`
          for more details)

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x + y^2; f
            (...0000000001)*y^2 + (...00000000010)*x

            sage: g = f.restriction(-1); g
            (...0000000001)*y^2 + (...00000000010)*x
            sage: g.parent()
            Tate Algebra in x (val >= 1), y (val >= 1) over 2-adic Field with capped relative precision 10

        Note that restricting may change the order of the terms::

            sage: f.restriction([-1,-2])
            (...00000000010)*x + (...0000000001)*y^2

        """
        parent = self._parent
        from sage.rings.tate_algebra import TateAlgebra
        ring = TateAlgebra(self.base_ring(), names=parent.variable_names(), log_radii=log_radii, 
                           prec=parent.precision_cap(), order=parent.term_order())
        return ring(self)

    def terms(self):
        r"""
        Return a list of the terms of this series sorted in descending order.

        .. NOTE::

            The order on the terms is defined as follows: first we
            compare the valuation and, in case of equality, we compare
            the monomials with respect to the order given at the
            creation of the parent.

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x^2 + x
            sage: f.terms()
            [(...0000000001)*x, (...00000000010)*x^2]

        """
        if not self._is_normalized:
            self._normalize()
            self._terms = None
        return self._terms_c()

    cdef list _terms_c(self):
        r"""
        Return a list of the terms of this series sorted in descending order.

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x^2 + x
            sage: f.terms()  # indirect doctest
            [(...0000000001)*x, (...00000000010)*x^2]

        """
        if self._terms is not None:
            return self._terms
        cdef TateAlgebraTerm oneterm = self._parent._oneterm
        cdef TateAlgebraTerm term
        self._terms = []
        for (e,c) in self._poly.__repn.iteritems():
            term = oneterm._new_c()
            term._coeff = c
            term._exponent = e
            if term.valuation() < self._prec:
                self._terms.append(term)
        self._terms.sort(reverse=True)
        return self._terms

    def monomials(self):
        r"""
        Return a list of the monomials of this series.

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 2*x^2 + x
            sage: f.monomials()  # indirect doctest
            [(...0000000001)*x, (...0000000001)*x^2]

        """
        return [ t.monomial() for t in self.terms() ]

    def coefficients(self):
        r"""
        Return the list of coefficients of this series.

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2
            sage: f.coefficients()
            [...0000000001, ...00000000010]

        """
        return [ t.coefficient() for t in self.terms() ]

    def add_bigoh(self, n):
        r"""
        Return this series truncated at precision ``n``.

        INPUT:

        - ``n`` -- an integer

        EXAMPLES::

            sage: R = Zp(2,prec=10,print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 32*x + 64*x^2; f
            (...000000000100000)*x + (...0000000001000000)*x^2
            sage: f.add_bigoh(5)
            O(2^5)

            sage: g = f.add_bigoh(6); g 
            (...100000)*x + O(2^6)
            sage: g.precision_absolute()
            6

        """
        return self._parent(self, prec=n)

    def precision_absolute(self):
        r"""
        Return the maximal precision at which a term of this series is known.

        EXAMPLES::
        
            sage: R = Zp(2,prec=10,print_mode='digits')
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x + 2*x^2; f
            (...0000000001)*x + (...00000000010)*x^2
            sage: f.precision_absolute()
            +Infinity

            sage: g = f.add_bigoh(5); g
            (...00001)*x + (...00010)*x^2 + O(2^5)
            sage: g.precision_absolute()
            5

        The absolute precision may be higher than the precision of some
        individual coefficients::

            sage: g = f.add_bigoh(20); g
            (...0000000001)*x + (...00000000010)*x^2 + O(2^20)
            sage: g.precision_absolute()
            20

        """
        return self._prec

    cpdef valuation(self):
        r"""
        Return the valuation of this series.

        .. NOTE::

            The valuation of a series `f` is defined as the minimal
            valuation of `f(x)` for `x` varying in the domain of convergence
            (specified in the parent).

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^4 + 4*x*y + 1; f
            (...0000000001)*x^4 + (...0000000001) + (...000000000100)*x*y
            sage: f.valuation()
            0

            sage: g = 2*f; g
            (...00000000010)*x^4 + (...00000000010) + (...0000000001000)*x*y
            sage: g.valuation()
            1

        When the radius of convergence is not 1, the variables themselves
        have a nontrivial valuation::

            sage: A.<x,y> = TateAlgebra(R, log_radii=(1,2))
            sage: x.valuation()
            -1
            sage: y.valuation()
            -2

            sage: f = x^4 + 4*x*y + 1
            sage: f.valuation()
            -4
          
        """
        cdef TateAlgebraTerm t
        cdef list terms = self._terms_c()
        if terms:
            return min(terms[0].valuation(), self._prec)
        else:
            return self._prec

    def precision_relative(self):
        """
        Return the relative precision of this series.

        The relative precision is defined as the difference 
        between the absolute precision and the valuation.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10)
            sage: A.<x,y> = TateAlgebra(R.fraction_field())
            sage: f = x^4 + 4*x*y + 1; f
            (...0000000001)*x^4 + (...0000000001) + (...000000000100)*x*y
            sage: f.precision_relative()
            +Infinity

            sage: g = f.add_bigoh(5)
            sage: g.precision_relative()
            5
            sage: g.precision_absolute()
            5
            sage: g.valuation()
            0

            sage: h = g + 1/2 ; h
            (...00001.1) + (...00001)*x^4 + (...00100)*x*y + O(2^5)
            sage: h.precision_relative()
            6
            sage: h.precision_absolute()
            5
            sage: h.valuation()
            -1
        """
        return self._prec - self.valuation()

    def leading_term(self):
        r"""
        Return the leading term of this series.

        .. NOTE::

            The order on the terms is defined as follows: first we
            compare the valuation and, in case of equality, we compare
            the monomials with respect to the order given at the
            creation of the parent.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^4 + x*y + 1; f
            (...0000000001)*x^4 + (...0000000001)*x*y + (...0000000001)
            sage: f.leading_term()
            (...0000000001)*x^4

            sage: g = f + x^4; g
            (...0000000001)*x*y + (...0000000001) + (...0000000010)*x^4
            sage: g.leading_monomial()
            (...0000000001)*x*y

        Observe that the leading term may change after restriction::

            sage: f.restriction(-1).leading_term()
            (...0000000001)
        
        """
        terms = self.terms()
        if terms:
            return terms[0]
        else:
            raise ValueError("zero has no leading term")

    def leading_coefficient(self):
        """
        Return the leading coefficient of this series.

        .. NOTE::

            The leading coefficient is the coefficient of the leading term.

        EXAMPLES::

            sage: R = Zp(2, prec=10, print_mode="terse")
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^4 + 3*x*y + 1; f
            (1 + O(2^10))*x^4 + (3 + O(2^10))*x*y + (1 + O(2^10))
            sage: f.leading_coefficient()
            1 + O(2^10)

            sage: g = f + x^4; g
            (3 + O(2^10))*x*y + (1 + O(2^10)) + (2 + O(2^10))*x^4
            sage: g.leading_coefficient()
            3 + O(2^10)
        
        """
        terms = self.terms()
        if terms:
            return terms[0].coefficient()
        else:
            return self.base_ring()(0)

    def leading_monomial(self):
        """
        Return the leading coefficient of this series.

        .. NOTE::

            The leading monomial is the monomial of the leading term.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^4 + x*y + 1; f
            (...0000000001)*x^4 + (...0000000001)*x*y + (...0000000001)
            sage: f.leading_monomial()
            (...0000000001)*x^4

            sage: g = f + x^4; g
            (...0000000001)*x*y + (...0000000001) + (...0000000010)*x^4
            sage: g.leading_monomial()
            (...0000000001)*x*y

        """
        terms = self.terms()
        if terms:
            return terms[0].monomial()
        else:
            raise ValueError("zero has no leading monomial")

    cpdef TateAlgebraElement monic(self):
        r"""
        Return this series normalized so that it has valuation 0 
        and its leading coefficient is a power of the uniformizer.

        EXAMPLES:

        When the log radii of convergence are all zero, the
        leading coefficient of the returned series is `1`::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R, order='lex')
            sage: f = 3*x^2*y^2 + 4*x^2*y^3 + y^5; f
            (...0000000011)*x^2*y^2 + (...0000000001)*y^5 + (...000000000100)*x^2*y^3
            sage: f.monic()
            (...0000000001)*x^2*y^2 + (...1010101011)*y^5 + (...101010101100)*x^2*y^3

        However, when log radii do not vanish, behaviors might
        be different::

            sage: g = f.restriction(-1); g
            (...0000000011)*x^2*y^2 + (...0000000001)*y^5 + (...000000000100)*x^2*y^3
            sage: g.monic()
            (...000000.0001)*x^2*y^2 + (...101010.1011)*y^5 + (...10101010.11)*x^2*y^3
            sage: g.monic().valuation()
            0

        TESTS::

            sage: A(0).monic()
            0

        """
        cdef TateAlgebraElement ans = self._new_c()
        cdef TateAlgebraTerm t
        cdef long shi
        cdef pAdicGenericElement u
        cdef list terms = self._terms_c()
        if not terms:
            return self
        t = terms[0]
        shi, u = t._coeff.val_unit()
        shi -= t._exponent.dotprod(self._parent._log_radii)
        ans._poly = self._poly.scalar_lmult((~u) >> shi)
        ans._prec = self._prec - shi
        return ans

    def is_monic(self):
        """
        Return ``True`` if this series is monic, in the sense
        that it has valuation 0 and its leading coefficient is
        a power of the uniformizer.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x*y + 2; f
            (...0000000001)*x*y + (...00000000010)
            sage: f.is_monic()
            True

            sage: g = f.restriction(-1); g
            (...00000000010) + (...0000000001)*x*y
            sage: g.is_monic()
            False

        """
        if self.valuation() != 0:
            return False
        c = self.leading_coefficient()
        if c == 0 or c.unit_part() == 1:
            return True
        return False


    def weierstrass_degree(self):
        r"""
        Return the Weierstrass degree of this Tate series.

        .. NOTE::

            The Weierstrass degree is the total degree of the polynomial
            defined by the terms with least valuation in the series.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x*y + 2*x^4 + y; f
            (...0000000001)*x*y + (...0000000001)*y + (...00000000010)*x^4
            sage: f.weierstrass_degree()
            2
        
        """
        v = self.valuation()
        return self.residue(v+1).degree()

    def degree(self):
        r"""
        Return the Weierstrass degree of this Tate series.

        .. NOTE::

            The Weierstrass degree is the total degree of the polynomial
            defined by the terms with least valuation in the series.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x*y + 2*x^4 + y; f
            (...0000000001)*x*y + (...0000000001)*y + (...00000000010)*x^4
            sage: f.weierstrass_degree()
            2
        
        """
        return self.weierstrass_degree()

    def weierstrass_degrees(self):
        r"""
        Return the Weierstrass degrees of this Tate series.

        .. NOTE::

            The Weierstrass degrees are the partial degrees of the polynomial
            defined by the terms with least valuation in the series.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^2 + y^2 + 2*x^3 + y; f
            (...0000000001)*x^2 + (...0000000001)*y^2 + (...0000000001)*y + (...00000000010)*x^3
            sage: f.weierstrass_degrees()
            (2, 2)
        
        """
        v = self.valuation()
        return self.residue(v+1).degrees()

    def degrees(self):
        r"""
        Return the Weierstrass degrees of this series.

        .. NOTE::

            The Weierstrass degrees are the partial degrees of the polynomial
            defined by the terms with least valuation in the series.

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits',prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^2 + y^2 + 2*x^3 + y; f
            (...0000000001)*x^2 + (...0000000001)*y^2 + (...0000000001)*y + (...00000000010)*x^3
            sage: f.weierstrass_degrees()
            (2, 2)
        
        """
        return self.weierstrass_degrees()

    def residue(self, n=None):
        r"""
        Return this series modulo the ``n``-th power of the uniformizer.

        Note that by definition of Tate series, the output is a polynomial.

        INPUT:

        - ``n`` -- an integer (default: ``1``)

        EXAMPLES::

            sage: R = Zp(2, print_mode='digits', prec=10)
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^2 + y^2 + 6*x^3 + 3*y 
            sage: f.residue()
            x^2 + y^2 + y
            sage: f.residue().parent()
            Multivariate Polynomial Ring in x, y over Finite Field of size 2

            sage: f.residue(2)
            2*x^3 + x^2 + y^2 - y
            sage: f.residue(2).parent()
            Multivariate Polynomial Ring in x, y over Ring of integers modulo 4

        The residue can only be computed for series with non-negative valuation.
        
            sage: g = f >> 2; g
            (...00000000.01)*x^2 + (...00000000.01)*y^2 + (...00000000.11)*y + (...000000001.1)*x^3
            sage: g.residue()
            Traceback (most recent call last):
            ...
            ValueError: element must have non-negative valuation in order to compute residue.

        The residue is not implemented for series with convergence radius different from 1.

            sage: A.<x,y> = TateAlgebra(R, log_radii=(2,-1))
            sage: f = x^2 + y^2 + 6*x^3 + 3*y 
            sage: f.residue()
            Traceback (most recent call last):
            ...
            NotImplementedError: residues are only implemented for radius 1

        """
        for r in self._parent.log_radii():
            if r != 0:
                raise NotImplementedError("residues are only implemented for radius 1")
        if n is None:
            Rn = self.base_ring().residue_field()
        else:
            try:
                Rn = self.base_ring().residue_ring(n)
            except (AttributeError, NotImplementedError):
                Rn = self.base_ring().change(field=False, type="fixed-mod", prec=n)
        poly = self._parent._polynomial_ring(self._poly)
        return poly.change_ring(Rn)

    cdef _quo_rem_c(self, list divisors, bint quo, bint rem, bint integral):
        r"""
        Perform the division of this series by ``divisors``.

        INPUT:

        - ``divisors`` -- the list of divisors

        - ``quo`` -- a boolean, whether we should compute the quotients

        - ``rem`` -- a boolean, whether we should compute the remainder

        TESTS::

            sage: R = Zp(2, 5, print_mode="digits")
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 1 + 2*x*y + 3*x^2*y + 4*x*y^2
            sage: g = x^2
            sage: q, r = f.quo_rem(g)  # indirect doctest
            sage: q
            (...00011)*y
            sage: r
            (...00001) + (...00010)*x*y + (...00100)*x*y^2 + O(2^5)

        """
        cdef dict coeffs = { }
        cdef TateAlgebraElement f
        cdef TateAlgebraTerm lt
        cdef list ltds = [ (<TateAlgebraElement>d)._terms_c()[0] for d in divisors ]
        cdef list quos = [ ]
        cdef list terms = self._terms_c()
        cdef int index = 0

        if quo:
            for d in divisors:
                f = self._new_c()
                f._poly = PolyDict({})
                f._prec = d.precision_relative()
                quos.append(f)

        f = self._new_c()
        f._poly = PolyDict(self._poly.__repn)
        f._prec = self._prec
        while len(terms) > index:
            lt = terms[index]
            if lt._valuation_c() >= f._prec:
                break
            for i in range(len(divisors)):
                if (<TateAlgebraTerm>ltds[i])._divides_c(lt, integral=integral):
                    factor = lt._floordiv_c(<TateAlgebraTerm>ltds[i])
                    f = f - (<TateAlgebraElement>divisors[i])._term_mul_c(factor)
                    terms = f._terms_c()
                    index = 0
                    if quo:
                        quos[i] = (<TateAlgebraElement>quos[i]) + factor
                    break
            else:
                if rem:
                    if coeffs.has_key(lt._exponent):
                        coeffs[lt._exponent] += lt._coeff
                    else:
                        coeffs[lt._exponent] = lt._coeff
                del f._poly.__repn[lt._exponent]
                index += 1
        if rem:
            f._poly += PolyDict(coeffs, self._poly.__zero)
            f._terms = None
        return quos, f

    cdef _quo_rem_check(self, divisors, bint quo, bint rem):
        """
        Perform the division of this series by ``divisors``.

        INPUT:

        - ``divisors`` -- the list of divisors

        - ``quo`` -- a boolean, whether we should compute the quotients

        - ``rem`` -- a boolean, whether we should compute the remainder

        TESTS::

            sage: R = Zp(2, 10)
            sage: A.<u,v> = TateAlgebra(R)
            sage: AA = A.integer_ring()

            sage: f = AA(u^2 + 2*v^2)
            sage: f.quo_rem(u)  # indirect doctest
            ((1 + O(2^10))*u, (2 + O(2^10))*v^2 + O(2^10))

        We check that coercion works::

            sage: f % 2  # indirect doctest
            (1 + O(2^10))*u^2 + O(2^10)
            sage: f % [2,u]  # indirect doctest
            O(2^10)

            sage: S.<pi> = R.extension(x^2 - 2)
            sage: f % (pi*u)  # indirect doctest
            (pi^2 + O(pi^20))*v^2 + O(pi^20)
            sage: (pi*f) // (pi*u) == f // u
            True

       ::

            sage: s = Zp(3)(1)
            sage: f % s
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce all the elements to the same parent

            sage: f % ZZ
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce all the elements to the same parent

        """
        parent = self.parent()
        if self.precision_absolute() is Infinity:
            self = self.add_bigoh(self.valuation() + parent.precision_cap())
        if isinstance(divisors, (list, tuple)):
            onedivisor = False
        else:
            divisors = [divisors]
            onedivisor = True
        A = _pushout_family(divisors, self._parent)
        f = A(self)
        divisors = [A(d) for d in divisors]
        q, r = (<TateAlgebraElement>self)._quo_rem_c(divisors, quo, rem, False)
        if quo and onedivisor:
            q = q[0]
        if quo and rem:
            return q, r
        if quo:
            return q
        if rem:
            return r

    def quo_rem(self, divisors):
        """
        Return the quotient(s) and the remainder of the division of
        this series by ``divisors``.

        INPUT:

        - ``divisors`` -- a series, or a list of series

        NOTE::

        The condition on the remainder is that it has 

        - no term which is greater than the leading term of the 
          numerator and

        - no term which is divisible by the leading term of one
          divisor.

        EXAMPLES::

            sage: R = Zp(2, 5, print_mode="digits")
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 1 + 2*x*y + 3*x^2*y + 4*x*y^2
            sage: g = x^2
            sage: q, r = f.quo_rem(g)
            sage: q
            (...00011)*y
            sage: r
            (...00001) + (...00010)*x*y + (...00100)*x*y^2 + O(2^5)
            sage: f == g*q + r
            True

        We can also divide by a family of divisors::

            sage: g0 = x^2
            sage: g1 = x*y + 2*x
            sage: q, r = f.quo_rem([g0, g1])
            sage: q
            [(...00011)*y, (...11010) + (...00100)*y]
            sage: r
            (...00001) + (...01100)*x + O(2^5)
            sage: f == g0*q[0] + g1*q[1] + r
            True

        """
        return (<TateAlgebraElement>self)._quo_rem_check(divisors, True, True)

    def __mod__(self, divisors):
        """
        Return the remainder of the division of this series by 
        ``divisors``.

        INPUT:

        - ``divisors`` -- a series, or a list of series

        EXAMPLES::

            sage: R = Zp(2, 5, print_mode="digits")
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 1 + 2*x*y + 3*x^2*y + 4*x*y^2
            sage: g = x^2
            sage: f % g
            (...00001) + (...00010)*x*y + (...00100)*x*y^2 + O(2^5)

        We can also divide by a family of divisors::

            sage: g0 = x^2
            sage: g1 = x*y + 2*x
            sage: f % [g0, g1]
            (...00001) + (...01100)*x + O(2^5)

        """
        return (<TateAlgebraElement>self)._quo_rem_check(divisors, False, True)

    def __floordiv__(self, divisors):
        """
        Return the quotient(s) of the division of this series by 
        ``divisors``.

        INPUT:

        - ``divisors`` -- a series, or a list of series

        EXAMPLES::

            sage: R = Zp(2, 5, print_mode="digits")
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 1 + 2*x*y + 3*x^2*y + 4*x*y^2
            sage: g = x^2
            sage: f // g
            (...00011)*y

        We can also divide by a family of divisors::

            sage: g0 = x^2
            sage: g1 = x*y + 2*x
            sage: f // [g0, g1]
            [(...00011)*y, (...11010) + (...00100)*y]

        """
        return (<TateAlgebraElement>self)._quo_rem_check(divisors, True, False)

    def reduce(self, I):
        """
        Return a canonical representative of this series in the
        quotient of the Tate algebra (in which this series lives)
        by the ideal ``I``.

        EXAMPLES::

            sage: R = Zp(3, prec=10, print_mode="digits")
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = 3*x^2 + 5*x*y^2
            sage: g = 5*x^2*y + 3
            sage: I = A.ideal([f, g])

            sage: f.reduce(I)
            O(3^9)
            sage: h = (x^2 + 2*y)*f + (x^2*y^3 + 3*x*y^2 + 7)*g + 1
            sage: h.reduce(I)
            (...000000001) + O(3^9)

        TESTS::

            sage: s = I.random_element(integral=True)
            sage: s.reduce(I)
            O(3^9)

            sage: h = A.random_element()
            sage: (h + s).reduce(I) == h.reduce(I)
            True

        """
        return self % I.groebner_basis()

    @coerce_binop
    def Spoly(self, other):
        """
        Return the S-polynomial of this series and ``other``.

        INPUT:

        - ``other`` -- a Tate series

        EXAMPLES::

            sage: R = Zp(2, 5, print_mode="digits")
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^3*y + 2*x*y + 4*x^2
            sage: g = 2*x*y^2 + 2*x
            sage: h = f.Spoly(g); h
            (...111110)*x^3 + (...0000100)*x*y^2 + (...00001000)*x^2*y

            sage: h == 2*y*f - x^2*g
            True

        TESTS::

            sage: f.Spoly(0)
            Traceback (most recent call last):
            ...
            ValueError: the S-polynomial of zero is not defined

        """
        try:
            return self._Spoly_c(other)
        except IndexError:
            raise ValueError("the S-polynomial of zero is not defined")

    cdef TateAlgebraElement _Spoly_c(self, TateAlgebraElement other):
        """
        Return the S-polynomial of this series and ``other``.

        INPUT:

        - ``other`` -- a Tate series

        TESTS::

        We check that the S-polynomial of two monomials vanishes::

            sage: R = Zp(3, 5, print_mode="digits")
            sage: A.<x,y> = TateAlgebra(R)
            sage: f = x^3*y^2
            sage: g = x^2*y^3
            sage: f.Spoly(g)
            0

        """
        cdef TateAlgebraTerm st = self._terms_c()[0]
        cdef TateAlgebraTerm ot = other._terms_c()[0]
        cdef TateAlgebraTerm t = st._lcm_c(ot)
        return self._term_mul_c(t._floordiv_c(st)) - other._term_mul_c(t._floordiv_c(ot))
