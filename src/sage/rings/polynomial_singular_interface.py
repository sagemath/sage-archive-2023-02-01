"""
Polynomial Interfaces to Singular

AUTHORS:
     -- Martin Albrecht <malb@informatik.uni-bremen.de> (2006-04-21)

"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
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
from real_field import is_RealField

class PolynomialRing_singular_repr:
    """
    Implements methods to convert polynomial rings to Singular.

    This class is a base class for all univariate and multivariate
    polynomial rings which support conversion from and to Singular
    rings.
    """
    def _singular_(self, singular=singular_default):
        """
        Returns a singular ring for this polynomial ring
        over a field. Currently QQ, GF(p), and GF(p^n) are
        supported. Singular also supports C and RR but these
        are unsupported right now.

        INPUT:
            singular -- Singular instance

        OUTPUT:
            singular ring matching this ring

        EXAMPLES:
            sage: r=MPolynomialRing(GF(2**8),10,'x')
            sage: r._singular_()
            //   characteristic : 2
            //   1 parameter    : a
            //   minpoly        : (a^8+a^4+a^3+a^2+1)
            //   number of vars : 10
            //        block   1 : ordering rp
            //                  : names    x0 x1 x2 x3 x4 x5 x6 x7 x8 x9
            //        block   2 : ordering C
            sage: r=MPolynomialRing(GF(127),2,'x')
            sage: r._singular_()
            //   characteristic : 127
            //   number of vars : 2
            //        block   1 : ordering rp
            //                  : names    x0 x1
            //        block   2 : ordering C
            sage: r=MPolynomialRing(QQ,2,'x')
            sage: r._singular_()
            //   characteristic : 0
            //   number of vars : 2
            //        block   1 : ordering rp
            //                  : names    x0 x1
            //        block   2 : ordering C
            sage: r=PolynomialRing(QQ)
            sage: r._singular_()
            //   characteristic : 0
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C
            sage: r=PolynomialRing(GF(127))
            sage: r._singular_()
            //   characteristic : 127
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C
            sage: r=PolynomialRing(GF(2**8),'y')
            sage: r._singular_()
            //   characteristic : 2
            //   1 parameter    : a
            //   minpoly        : (a^8+a^4+a^3+a^2+1)
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    y
            //        block   2 : ordering C

        WARNING:
           If the base ring is a finite extension field the ring will
           not only be returned but also be set as the current ring in
           Singular.
        """
        try:
            R = self.__singular
            if not (R.parent() is singular):
                raise ValueError
            R._check_valid()
            if self.base_ring().is_prime_field():
                return R
            if self.base_ring().is_finite():
                R.set_ring() #sorry for that, but needed for minpoly
                if  singular.eval('minpoly') != self.__minpoly:
                    singular.eval("minpoly=%s"%(self.__minpoly))
            return R
        except (AttributeError, ValueError):
            return self._singular_init_(singular)

    def _singular_init_(self, singular=singular_default):
        """
        Return a newly created Singular ring matching this ring.
        """
        if not self._can_convert_to_singular():
            raise TypeError, "no conversion of %s to a Singular ring defined"%self

        if self.ngens()==1:
            _vars = str(self.gen())
            order = 'lp'
        else:
            _vars = str(self.gens())
            order = self.term_order().singular_str()

        if is_RealField(self.base_ring()):
            # TODO: here we would convert the SAGE bit precision to a format
            # Singular understands, and would call
            # self.__singular = singular.ring("(real,<PREC>", _vars )
            self.__singular = singular.ring("real", _vars, order=order)

        elif is_ComplexField(self.base_ring()):
            # TODO: here we would convert the SAGE precision to a format
            # Singular understands, and would call
            # self.__singular = singular.ring("(complex,<PREC>,I", _vars )
            self.__singular = singular.ring("(complex,I)", _vars,  order=order)

        elif self.base_ring().is_prime_field():
            self.__singular = singular.ring(self.characteristic(), _vars, order=order)
            return self.__singular

        elif self.base_ring().is_finite(): #must be extension field
            gen = str(self.base_ring().gen())
            r = singular.ring( "(%s,%s)"%(self.characteristic(),gen), _vars, order=order)
            self.__minpoly = "("+(str(self.base_ring().modulus()).replace("x",gen)).replace(" ","")+")"
            singular.eval("minpoly=%s"%(self.__minpoly) )

            self.__singular = r
        else:
            raise TypeError, "no conversion of %s to a Singular ring defined"%self
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
                 or is_ComplexField(base_ring) )


class Polynomial_singular_repr:
    """
    Implements coercion of polynomials to Singular polynomials.

    This class is a base class for all (univariate and multivariate)
    polynomial classes which support conversion from and to
    Singular polynomials.
    """
    def _singular_(self, singular=singular_default):
        """
        Return Singular polynomial matching this polynomial.

        INPUT:
            singular -- Singular instance to use

        EXAMPLES:
            sage: R = PolynomialRing(GF(7))
            sage: x = R.gen()
            sage: f = (x^3 + 2*x^2*x)^7; f
            3*x^21
            sage: h = f._singular_(); h
            3*x^21
            sage: R(h)
            3*x^21
            sage: R(h^20) == f^20
            True
            sage: R = PolynomialRing(GF(7), 2, ['x','y'])
            sage: x, y = R.gens()
            sage: f = (x^3 + 2*y^2*x)^7; f
            2*x^7*y^14 + x^21
            sage: h = f._singular_(); h
            x^21+2*x^7*y^14
            sage: R(h)
            2*x^7*y^14 + x^21
            sage: R(h^20) == f^20
            True
        """
        self.parent()._singular_(singular).set_ring() #this is expensive
        try:
            if self.__singular.parent() is singular:
                return self.__singular
        except AttributeError:
            pass
        return self._singular_init_(singular)

    def _singular_init_(self,singular=singular_default):
        """
        Return corresponding Singular polynomial but enforce that a new
        instance is created in the Singular interpreter.
        """
        self.parent()._singular_(singular).set_ring() #this is expensive
        self.__singular = singular(str(self))
        return self.__singular
