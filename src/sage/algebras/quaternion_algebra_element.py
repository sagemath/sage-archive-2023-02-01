"""
Quaternion algebra elements
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import operator
from sage.misc.misc import repr_lincomb
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.rational_field import RationalField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.algebras.free_algebra_quotient_element import FreeAlgebraQuotientElement
from sage.rings.infinity import Infinity

class QuaternionAlgebraElement(FreeAlgebraQuotientElement):

    def __init__(self, H, x):
        """
        Create the element x of the quaternion algebra H.
        """
        FreeAlgebraQuotientElement.__init__(self, H, x)

    def __cmp__(self, x):
        return cmp(self.vector(), x.vector())

    def conjugate(self):
        """
        Return the conjugate of this element.

        EXAMPLES:
            sage: A.<i,j,k>=QuaternionAlgebra(QQ,-5,-2)
            sage: x=3*i-j+2
            sage: x.conjugate()
            2 - 3*i + j
        """
        return self.parent()(self.reduced_trace()) - self

    def reduced_trace(self):
        r"""
        Return the reduced trace of this element.

        \note{In a quaternion algebra $A$, every element $x$ is
        quadratic over the center, thus $x^2 = \Tr(x)*x - \Nr(x)$, so
        we solve for a linear relation $(1,-\Tr(x),\Nr(x))$ among
        $[x^2, x, 1]$ for the reduced trace of $x$.}

        EXAMPLES:
            sage: A.<i,j,k>=QuaternionAlgebra(QQ,-5,-2)
            sage: x=3*i-j+2
            sage: x.reduced_trace()
            4
        """
        v = self.vector()
        if v[1] == 0 and v[2] == 0 and v[3] == 0: return 2*v[0]
        u = (self**2).vector()
        K = self.parent().base_ring()
        A = MatrixSpace(K,3,4)
        M = A(list(u) + list(v) + [1,0,0,0]).kernel()
        w = M.gen(0)
        if w[0] == 1: return -w[1]
        return -w[1]/w[0]

    def reduced_norm(self):
        """
        Return the reduced norm of this element.

        EXAMPLES:
            sage: A.<i,j,k>=QuaternionAlgebra(QQ,-5,-2)
            sage: x=3*i-j+2
            sage: x.reduced_norm()
            51
        """
        x = self * self.conjugate()
        return x.vector()[0]

    def charpoly(self, var):
        """
        Return the characteristic polynomial of this element in terms
        of the given variable.

        EXAMPLES:
            sage: A.<i,j,k>=QuaternionAlgebra(QQ,-5,-2)
            sage: x=3*i-j+2
            sage: x.charpoly('t')
            t^2 - 4*t + 51
        """
        v = self.vector()
        if v[1] == 0 and v[2] == 0 and v[3] == 0:
            return 2*v[0]
        u = (self**2).vector()
        A = MatrixSpace(RationalField(),3,4)
        M = A(list(u) + list(v) + [1,0,0,0]).kernel()
        w = M.gen(0)
        P = PolynomialRing(self.parent().base_ring(), var)
        x = P.gen()
        if w[0] == 1:
            x**2 + w[1]*x + w[2]
            return x**2 + w[1]/w[0]*x + w[2]/w[0]

        return x**2 - self.reduced_trace()*x + self.reduced_norm()

    characteristic_polynomial = charpoly

    def minpoly(self, var):
        """
        Return the minimal polynomial of this element in terms
        of the given variable.

        EXAMPLES:
            sage: A.<i,j,k>=QuaternionAlgebra(QQ,-5,-2)
            sage: x=3*i-j+2
            sage: x.minpoly('t')
            t^2 - 4*t + 51
        """
        v = self.vector()
        if v[1] == 0 and v[2] == 0 and v[3] == 0:
            K = self.parent().base_ring()
            P = PolynomialRing(K, var)
            x = P.gen()
            return x - v[0]
        return self.charpoly(var)

    minimal_polynomial = minpoly

    def is_unit(self):
        """
        Return True if the element is an invertible element of the
        quaternion algebra.

        EXAMPLES:
            sage: A.<i,j,k> = QuaternionAlgebra(QQ,-1,-1)
            sage: i.is_unit()
            True
            sage: (i-5+j*k).is_unit()
            True
            sage: A(0).is_unit()
            False
        """
        if self.reduced_norm() == 0:
            return False
        else:
            return True

    def additive_order(self):
        if self.base_ring().is_finite():
            return self.base_ring().characteristic()
        else:
            return Infinity

    def is_scalar(self):
        """
        Return True is this element of a quaternion algebra is
        a scalar (i.e. lies in the base field).

        EXAMPLES:
            sage: A.<i,j,k> = QuaternionAlgebra(QQ,-1,-1)
            sage: i.is_scalar()
            False
            sage: (i-5+j*k).is_scalar()
            False
            sage: A(12).is_scalar()
            True
        """
        if (self.reduced_trace()-2*self).is_zero():
            return True
        else:
            return False

    def is_pure(self):
        """
        Return True is this element of a quaternion algebra is
        "pure" (i.e. has no scalar component, or has reduced-trace zero).

        EXAMPLES:
            sage: A.<i,j,k> = QuaternionAlgebra(QQ,-1,-1)
            sage: i.is_pure()
            True
            sage: (i-5+j*k).is_pure()
            False
            sage: A(12).is_pure()
            False
        """
        if self.is_zero():
            return False
        if self.reduced_trace() == 0:
            return True
        return False

    def scalar_part(self):
        """
        Return the part of the quaternion 'self' that lies in the base field/ring.
        This is given by the reduced trace (note: we assume characteristic not 2).
        We could cheat, using self.vector(), but we really don't know what basis
        is in place.  Do we?
        EXAMPLES:
            sage: A.<i,j,k> = QuaternionAlgebra(QQ,-1,-1)
            sage: i.scalar_part()
            0
            sage: x = A([1,-3/2,0,2])
            sage: x.scalar_part()
             1
        """
        return self.reduced_trace()/2

    def pure_part(self):
        """
        Return the part of the quaternion 'self' that lies in the vector subspace
        "<i,j,k>" (figuratively speaking).  We just strip off the scalar part...
        EXAMPLES:
            sage: A.<i,j,k> = QuaternionAlgebra(QQ,-1,-1)
            sage: x = A([1,-3/2,0,2])
            sage: x.pure_part()
             -3/2*i + 2*k
        """
        return self - self.scalar_part()

    def _div_(self, other):
        """
        Right division in the quaternion algebra

        EXAMPLES:
            sage: A.<i,j,k>=QuaternionAlgebra(QQ,-1,-1)
            sage: x=3*i-j+2
            sage: y=i-1
            sage: x/y
            1/2 - 5/2*i + 1/2*j - 1/2*k

        Note that 1/x will raise an AttributeError.  The way to get
        the inverse of x is
            sage: A(1)/x
            1/7 - 3/14*i + 1/14*j
        """
        if other.is_unit() == False:
            raise AttributeError, "The second operand must be a unit"
        H = self.parent()
        if H(other).is_scalar():
            return QuaternionAlgebraElement(H, [self.vector()[i]/(H(other).reduced_trace()/2) for i in range(4)])
        elif self.is_scalar():
            return self*other.conjugate()/other.reduced_norm()
        else:
            return self*(H(1)/other)

    def _backslash_(self, other):
        """
        Left division in the quaternion algebra

        EXAMPLES:
            sage: A.<i,j,k>=QuaternionAlgebra(QQ,-1,-1)
            sage: x=3*i-j+2
            sage: y=i-1
            sage: x\y
            1/14 + 5/14*i - 1/14*j - 1/14*k
        """
        if self.is_unit() == False:
            raise AttributeError, "The first operand must be a unit"
        else:
            H = self.parent()
            return (H(1)/self)*other

class QuaternionAlgebraElement_fast(QuaternionAlgebraElement):

    def __init__(self, H, x):
        FreeAlgebraQuotientElement.__init__(self, H, x)
        if isinstance(x, list):
            self._list = list(self.vector())
        else:
            self._list = x

    def _add_(self, other):
        return QuaternionAlgebraElement_fast(self.parent(), [self._list[i] + other._list[i] for i in range(4)])

    def _sub_(self, other):
        return QuaternionAlgebraElement_fast(self.parent(), [self._list[i] - other._list[i] for i in range(4)])

    def _neg_(self, other):
        return QuaternionAlgebraElement_fast(self.parent(), [-r for r in self._list])

    def _mul_(self, other):
        d1, d2 = self.parent().discriminants
        a = self._list
        b = other._list
        return QuaternionAlgebraElement_fast(self.parent(),
            [a[0]*b[0] + d1*a[1]*b[1] + d2*a[2]*b[2] - d1*d2*a[3]*b[3],
             a[0]*b[1] +    a[1]*b[0] + d2*a[3]*b[1] -    d2*a[2]*b[3],
             a[0]*b[2] +    a[2]*b[0] + d1*a[1]*b[3] -    d1*a[3]*b[1],
             a[0]*b[3] +    a[3]*b[0] +    a[1]*b[2] -       a[2]*b[1]])


    def conjugate(self):
        return QuaternionAlgebraElement_fast(self.parent(), [self._list[0], -self._list[1], -self._list[2], -self._list[3]])


