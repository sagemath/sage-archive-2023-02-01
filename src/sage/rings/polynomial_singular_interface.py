"""
Polynomial Interfaces to Singular

AUTHORS:
     -- Martin Albrecht <malb@informatik.uni-bremen.de> (2006-04-21)

TESTS:
    sage: R = MPolynomialRing(GF(2**8,'a'),10,'x', order='revlex')
    sage: R == loads(dumps(R))
    True

"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import finite_field

from sage.interfaces.all import singular as singular_default, is_SingularElement
from complex_field import is_ComplexField
from real_mpfr import is_RealField
from complex_double import is_ComplexDoubleField
from real_double import is_RealDoubleField
from integer_ring import ZZ
import sage.rings.arith

class PolynomialRing_singular_repr:
    """
    Implements methods to convert polynomial rings to Singular.

    This class is a base class for all univariate and multivariate
    polynomial rings which support conversion from and to Singular
    rings.
    """
    def _singular_(self, singular=singular_default, force=False):
        """
        Returns a singular ring for this polynomial ring over a field.
        Currently QQ, GF(p), and GF(p^n), CC, and RR are supported.

        INPUT:
            singular -- Singular instance
            force -- polynomials over ZZ may be coerced to Singular by
                     treating them as polynomials over RR. This is
                     inexact but works for some cases where the
                     coeffients are not considered (default: False).

        OUTPUT:
            singular ring matching this ring

        EXAMPLES:
            sage: r = MPolynomialRing(GF(2**8,'a'),10,'x', order='revlex')
            sage: r._singular_()
            //   characteristic : 2
            //   1 parameter    : a
            //   minpoly        : (a^8+a^4+a^3+a^2+1)
            //   number of vars : 10
            //        block   1 : ordering rp
            //                  : names    x0 x1 x2 x3 x4 x5 x6 x7 x8 x9
            //        block   2 : ordering C
            sage: r = MPolynomialRing(GF(127),2,'x', order='revlex')
            sage: r._singular_()
            //   characteristic : 127
            //   number of vars : 2
            //        block   1 : ordering rp
            //                  : names    x0 x1
            //        block   2 : ordering C
            sage: r = MPolynomialRing(QQ,2,'x', order='revlex')
            sage: r._singular_()
            //   characteristic : 0
            //   number of vars : 2
            //        block   1 : ordering rp
            //                  : names    x0 x1
            //        block   2 : ordering C
            sage: r = PolynomialRing(QQ,'x')
            sage: r._singular_()
            //   characteristic : 0
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C
            sage: r = PolynomialRing(GF(127),'x')
            sage: r._singular_()
            //   characteristic : 127
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C
            sage: r = PolynomialRing(GF(2**8,'a'),'y')
            sage: r._singular_()
            //   characteristic : 2
            //   1 parameter    : a
            //   minpoly        : (a^8+a^4+a^3+a^2+1)
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    y
            //        block   2 : ordering C
            sage: R.<x,y> = PolynomialRing(CC,'x',2)
            sage: R._singular_()
            //   characteristic : 0 (complex:15 digits, additional 0 digits)
            //   1 parameter    : I
            //   minpoly        : (I^2+1)
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C
            sage: R.<x,y> = PolynomialRing(RealField(100),'x',2)
            sage: R._singular_()
            //   characteristic : 0 (real:29 digits, additional 0 digits)
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C

        WARNING:
           If the base ring is a finite extension field the ring will
           not only be returned but also be set as the current ring in
           Singular.

        NOTE:
            Singular represents precision of floating point numbers base 10
            while SAGE represents floating point precision base 2.
        """
        try:
            R = self.__singular
            if not (R.parent() is singular):
                raise ValueError
            R._check_valid()
            if self.base_ring() is ZZ or self.base_ring().is_prime_field():
                return R
            if self.base_ring().is_finite():
                R.set_ring() #sorry for that, but needed for minpoly
                if  singular.eval('minpoly') != self.__minpoly:
                    singular.eval("minpoly=%s"%(self.__minpoly))
            return R
        except (AttributeError, ValueError):
            return self._singular_init_(singular, force)

    def _singular_init_(self, singular=singular_default, force=False):
        """
        Return a newly created Singular ring matching this ring.
        """
        if not self._can_convert_to_singular() and not force:
            raise TypeError, "no conversion of this ring to a Singular ring defined"

        if self.ngens()==1:
            _vars = str(self.gen())
            if "*" in _vars: # 1.000...000*x
                _vars = _vars.split("*")[1]
            order = 'lp'
        else:
            _vars = str(self.gens())
            order = self.term_order().singular_str()

        if is_RealField(self.base_ring()):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            precision = self.base_ring().precision()
            digits = sage.rings.arith.integer_ceil((2*precision - 2)/7.0)
            self.__singular = singular.ring("(real,%d,0)"%digits, _vars, order=order, check=False)

        elif is_ComplexField(self.base_ring()):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            precision = self.base_ring().precision()
            digits = sage.rings.arith.integer_ceil((2*precision - 2)/7.0)
            self.__singular = singular.ring("(complex,%d,0,I)"%digits, _vars,  order=order, check=False)

        elif is_RealDoubleField(self.base_ring()):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            self.__singular = singular.ring("(real,15,0)", _vars, order=order, check=False)

        elif is_ComplexDoubleField(self.base_ring()):
            # singular converts to bits from base_10 in mpr_complex.cc by:
            #  size_t bits = 1 + (size_t) ((float)digits * 3.5);
            self.__singular = singular.ring("(complex,15,0,I)", _vars,  order=order, check=False)

        elif self.base_ring().is_prime_field() or (self.base_ring() is ZZ and force):
            self.__singular = singular.ring(self.characteristic(), _vars, order=order, check=False)

        elif self.base_ring().is_finite(): #must be extension field
            gen = str(self.base_ring().gen())
            r = singular.ring( "(%s,%s)"%(self.characteristic(),gen), _vars, order=order, check=False)
            self.__minpoly = "("+(str(self.base_ring().modulus()).replace("x",gen)).replace(" ","")+")"
            singular.eval("minpoly=%s"%(self.__minpoly) )

            self.__singular = r
        else:
            raise TypeError, "no conversion to a Singular ring defined"

        return self.__singular

    def _can_convert_to_singular(self):
        """
        Returns True if this rings base field/ring can be represented in
        Singular. If this is true then this polynomial ring can be
        represented in Singular.

        GF(p), GF(p^n), Rationals, Reals, and Complexes are supported.
        """
        base_ring = self.base_ring()
        return ( finite_field.is_FiniteField(base_ring)
                 or base_ring.is_prime_field()
                 or is_RealField(base_ring)
                 or is_ComplexField(base_ring)
                 or is_RealDoubleField(base_ring)
                 or is_ComplexDoubleField(base_ring)
                 or base_ring is ZZ )


class Polynomial_singular_repr:
    """
    Implements coercion of polynomials to Singular polynomials.

    This class is a base class for all (univariate and multivariate)
    polynomial classes which support conversion from and to
    Singular polynomials.
    """
    def _singular_(self, singular=singular_default, have_ring=False, force=False):
        """
        Return Singular polynomial matching this polynomial.

        INPUT:
            singular -- Singular instance to use

            have_ring -- if True we will not attempt to set this
                         element's ring as the current Singular
                         ring. This is useful to speed up a batch of
                         f._singular_() calls. However, it's dangerous
                         as it might lead to wrong results if another
                         ring is singluar.current_ring().  (default:
                         False)

            force -- polynomials over ZZ may be coerced to Singular by
                     treating them as polynomials over RR. This is
                     inexact but works for some cases where the
                     coeffients are not considered (default: False).


        EXAMPLES:
            sage: R.<x> = PolynomialRing(GF(7))
            sage: f = (x^3 + 2*x^2*x)^7; f
            3*x^21
            sage: h = f._singular_(); h
            3*x^21
            sage: R(h)
            3*x^21
            sage: R(h^20) == f^20
            True
            sage: R.<x,y> = PolynomialRing(GF(7), 2)
            sage: f = (x^3 + 2*y^2*x)^7; f
            x^21 + 2*x^7*y^14
            sage: h = f._singular_(); h
            x^21+2*x^7*y^14
            sage: R(h)
            x^21 + 2*x^7*y^14
            sage: R(h^20) == f^20
            True
        """
        if not have_ring:
            self.parent()._singular_(singular,force=force).set_ring() #this is expensive

        try:
            self.__singular._check_valid()
            if self.__singular.parent() is singular:
                return self.__singular
        except (AttributeError,ValueError):
            pass
        return self._singular_init_(singular,have_ring=have_ring)

    def _singular_init_(self, singular=singular_default, have_ring=False, force=False):
        """
        Return corresponding Singular polynomial but enforce that a new
        instance is created in the Singular interpreter.

        Use self._singular_() instead.
        """
        if not have_ring:
            self.parent()._singular_(singular,force=force).set_ring() #this is expensive

        self.__singular = singular(str(self))

        return self.__singular

    def lcm(self, right, have_ring=False):
        """
        Returns the least common multiple of this element and the right element.

        INPUT:
            right -- multivariate polynomial
            have_ring -- see self._singular_() (default:False)

        OUTPUT:
            multivariate polynomial representing the least common
            multiple of self and right

        ALGORITHM: Singular

        EXAMPLES:
            sage: r.<x,y> = MPolynomialRing(GF(2**8,'a'),2)
            sage: a = r.base_ring().0
            sage: f = (a^2+a)*x^2*y + (a^4+a^3+a)*y + a^5
            sage: f.lcm(x^4)
            (a^2 + a)*x^6*y + (a^4 + a^3 + a)*x^4*y + a^5*x^4
        """
        lcm = self._singular_(have_ring=have_ring).lcm(right._singular_(have_ring=have_ring))
        return lcm.sage_poly(self.parent())

    def diff(self, variable, have_ring=False):
        """
        Differentiates self with respect to the provided variable. This
        is completely symbolic so it is also defined over e.g. finite
        fields.

        INPUT:
            variable -- the derivative is taken with respect to variable
            have_ring -- see self._singular_() (default:False)

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(RR,2)
            sage: f = 3*x^3*y^2 + 5*y^2 + 3*x + 2
            sage: f.diff(x)
            9.00000000000000*x^2*y^2 + 3.00000000000000
            sage: f.diff(y)
            6.00000000000000*x^3*y + 10.0000000000000*y

            The derivate is also defined over finite fields:

            sage: R.<x,y> = PolynomialRing(GF(2**8, 'a'),2)
            sage: f = x^3*y^2 + y^2 + x + 2
            sage: f.diff(x)
            x^2*y^2 + 1

            The new coefficients are coerced to the base ring:

            sage: f.diff(y)
            0

        ALGORITHM: Singular

        """
        df = self._singular_(have_ring=have_ring).diff(variable._singular_(have_ring=have_ring))
        return df.sage_poly(self.parent())

##     def lt(self, have_ring=False):
##         """
##         Returns the leading (or initial) term of a polynomial
##         with respect to the monomial ordering.

##         INPUT:
##             have_ring -- see self.singular_() (default: False)

##         OUTPUT:
##             multivariate polynomial representing the lead term of self i.e.,
##             self.lc()*self.lm()

##         ALGORITHM: Singular

##         EXAMPLES:
##             sage: r=MPolynomialRing(GF(2**8),2,'x')
##             sage: x,y=r.gens()
##             sage: k=r.base_ring()
##             sage: f=k('a^2+a')*x^2 + k('a^4+a^3+a')*y^2 + k('a^5')
##             sage: f.lt()
##             (a^2 + a)*x0^2

##             sage: r=MPolynomialRing(GF(2**8),2,'x','deglex')
##             sage: x,y=r.gens()
##             sage: k=r.base_ring()
##             sage: f=k('a^2+a')*x^2 + k('a^4+a^3+a')*y^3 + k('a^5')
##             sage: f.lt()
##             (a^4 + a^3 + a)*x1^3

##         """
##         try:
##             return self.__lt
##         except AttributeError:
##             self.__lt = self._singular_(have_ring=have_ring).lead().sage_poly(self.parent())
##             return self.__lt

##     def lm(self, have_ring=False):
##         """
##         Returns the leading monomial of a multivariate polynomial as a
##         multivariate polynomial whose coefficient is one.

##         INPUT:
##             have_ring -- see self.singular_() (default: False)

##         OUTPUT:
##             multivariate polynomial representing the lead monomial of self

##         ALGORITHM: Singular

##         EXAMPLES:
##             sage: r=MPolynomialRing(GF(2**8),2,'x')
##             sage: x,y=r.gens()
##             sage: k=r.base_ring()
##             sage: f=k('a^2+a')*x^2 + k('a^4+a^3+a')*y^2 + k('a^5')
##             sage: f.lm()
##             x0^2

##             sage: r=MPolynomialRing(GF(2**8),2,'x','deglex')
##             sage: x,y=r.gens()
##             sage: k=r.base_ring()
##             sage: f=k('a^2+a')*x^2 + k('a^4+a^3+a')*y^3 + k('a^5')
##             sage: f.lm()
##             x1^3
##         """
##         try:
##             return self.__lm
##         except AttributeError:
##             self.__lm = self._singular_(have_ring=have_ring).leadmonom().sage_poly(self.parent())
##             return self.__lm

##     def lc(self, have_ring=False):
##         """
##         Returns the leading (or initial) coefficient of a polynomial
##         with respect to the monomial ordering.

##         INPUT:
##             have_ring -- see self.singular_() (default: False)

##         OUTPUT:
##             multivariate polynomial representing the lead coefficent of self

##         ALGORITHM: Singular

##         EXAMPLES:
##             sage: r=MPolynomialRing(GF(2**8),2,'x')
##             sage: x,y=r.gens()
##             sage: k=r.base_ring()
##             sage: f=k('a^2+a')*x^2 + k('a^4+a^3+a')*y^2 + k('a^5')
##             sage: f.lc()
##             a^2 + a

##             sage: r=MPolynomialRing(GF(2**8),2,'x','deglex')
##             sage: x,y=r.gens()
##             sage: k=r.base_ring()
##             sage: f=k('a^2+a')*x^2 + k('a^4+a^3+a')*y^3 + k('a^5')
##             sage: f.lc()
##             a^4 + a^3 + a

##             sage: R.<x,y,z> = PolynomialRing(QQ,3)
##             sage: f = (-1/3)*(1+x+y+z)
##             sage: (f^3).lc()
##             -1/27
##         """
##         try:
##             return self.__lc
##         except AttributeError:
##             c = self._singular_(have_ring=have_ring).leadcoef().sage_poly(self.parent())
##             self.__lc = self.base_ring()(c.constant_coefficient())
##             return self.__lc
