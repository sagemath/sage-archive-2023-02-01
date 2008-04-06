"""
Algebra order elements

AUTHOR: David Kohel, 2005-09
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
from sage.structure.element import bin_op
from sage.rings.ring import Ring
from sage.rings.integer import Integer
from sage.rings.ring_element import RingElement
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.modules.free_module_element import FreeModuleElement
from sage.algebras.algebra_element import AlgebraElement
from sage.algebras.free_algebra_quotient_element import FreeAlgebraQuotientElement
#from sage.algebras.free_algebra_quotient_element import is_FreeAlgebraQuotientElement

def is_AlgebraOrderElement(x):
    return isinstance(x, AlgebraOrderElement)

class AlgebraOrderElement(AlgebraElement):
    """
    """
    def __init__(self, A, x, check=True):
        """
        Create the element x of the algebra order A.
        """
        AlgebraElement.__init__(self, A)
        H = A.ambient_algebra()
        if isinstance(x, (int, long, Integer)):
            self.__ambient_algebra_element = H(x)
        elif isinstance(x, AlgebraOrderElement) and x.parent() is A:
            self.__ambient_algebra_element = x.__ambient_algebra_element
        elif isinstance(x, AlgebraElement) and x.parent() is H:
            if check and not x in A:
                raise ValueError, "Argument x (= %s) must be in the algebra order."%x
            self.__ambient_algebra_element = x
        elif isinstance(x, FreeModuleElement) and x.parent() is A.module():
            B = A.algebra_basis_elements()
            s = sum([ x[i] * B[i] for i in range(len(B)) ])
            self.__ambient_algebra_element = s
        elif isinstance(x, (list,tuple)):
            B = A.algebra_basis_elements()
            s = sum([ x[i] * B[i] for i in range(len(B)) ])
            self.__ambient_algebra_element = s
        elif isinstance(x, Ring):
            R = H.base_ring()
            try:
                x = R(x)
            except TypeError:
                raise TypeError, \
                      "Argument x (= %s) must be an ambient algebra element, list, or tuple"%x
            self.__ambient_algebra_element = H(x)
        else:
            raise TypeError, \
                  "Argument x (= %s) must be an ambient algebra element, list, or tuple"%x

    def __repr__(self):
        """
        """
        return str(self.__ambient_algebra_element)

    def ambient_algebra_element(self):
        return self.__ambient_algebra_element

    def trace(self):
        """
        The trace of the element with respect to its action
        by left or right multiplication on the order.
	"""
        R = self.parent().base_ring()
        return R(self.__ambient_algebra_element.trace())

    def norm(self):
        """
        The determinant (= norm) of the element with respect to
        its action by left or right multiplication on the order.
	"""
        R = self.parent().base_ring()
        return R(self.__ambient_algebra_element.norm())

    def charpoly(self, var):
        """
        The characteristic polynomial of the element with respect
        to its action by left or right multiplication on the order.
	"""
        R = self.parent().base_ring()
        P = PolynomialRing(R, var)
        return P(self.__ambient_algebra_element.charpoly('x'))

    characteristic_polynomial = charpoly

    def minpoly(self, var):
        """
        The minimal polynomial of the element with respect to its
        action by left or right multiplication on the order.
	"""
        R = self.parent().base_ring()
        P = PolynomialRing(R, var)
        return P(self.__element.minpoly())

    minimal_polynomial = minpoly

    def __add__(self,y):
        if not isinstance(y, (AlgebraOrderElement, FreeAlgebraQuotientElement)):
            return bin_op(self, y, operator.add)
        A = self.parent()
        x = self.__ambient_algebra_element
        if A is y.parent():
            return A(x + y.__ambient_algebra_element)
        Q = A.ambient_algebra()
        if Q is y.parent():
            return A(x + y)
        raise TypeError, "Argument y (= %s) is of the wrong type."%y

    def __radd__(self, y):
        """
        Can be deleted once this moves to the base ring class.
        """
        return bin_op(y, self, operator.add)

    def __neg__(self):
        return self.parent(-self.__ambient_algebra_element)

    def __sub__(self,y):
        if not isinstance(y, (AlgebraOrderElement, FreeAlgebraQuotientElement)):
            return bin_op(self, y, operator.sub)
        A = self.parent()
        x = self.__ambient_algebra_element
        if A is y.parent():
            return A(x - y.__ambient_algebra_element)
        Q = A.ambient_algebra()
        if Q is y.parent():
            return A(x - y)
        raise TypeError, "Argument y (= %s) is of the wrong type."%y

    def __mul__(self, y):
        if not isinstance(y, (AlgebraOrderElement, FreeAlgebraQuotientElement)):
            return bin_op(self, y, operator.mul)
        A = self.parent()
        x = self.__ambient_algebra_element
        if A is y.parent():
            return A(x * y.__ambient_algebra_element)
        Q = A.ambient_algebra()
        if Q is y.parent():
            return A(x * y)
        raise TypeError, "Argument y (= %s) is of the wrong type."%y




