"""
Quaternion order elements

AUTHORS:

- David Kohel (2005-09)
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

from sage.modules.free_module import FreeModule
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.algebras.algebra_order_element import AlgebraOrderElement

class QuaternionOrderElement(AlgebraOrderElement):
    """

    """
    def __init__(self, A, x, check=True):
        """
        Create the element x of the quaternion order A.
        """
        AlgebraOrderElement.__init__(self, A, x, check=check)

    def __eq__(self, x):
        H = self.parent()
        if not isinstance(x, QuaternionOrderElement) or not H is x.parent():
            try:
                x = H(x)
            except:
                return False
        return self.vector() == x.vector()

    def conjugate(self):
        return self.reduced_trace() - self

    def reduced_trace(self):
        r"""
	Return the reduced trace of this element.

	.. note::

           In a quaternion algebra `A`, every element
	   `x` is quadratic over the center, thus
           `x^2 = \mathrm{Tr}(x)*x - \mathrm{Nr}(x)`, so we solve for a linear
           relation `(1,-\mathrm{Tr}(x),\mathrm{Nr}(x))` among `[x^2, x, 1]`
           for the reduced trace of `x`.
	"""
        R = self.parent().base_ring()
        return R(self.ambient_algebra_element().reduced_trace())

    def reduced_norm(self):
        """

	"""
        R = self.parent().base_ring()
        return R(self.ambient_algebra_element().reduced_norm())

    def charpoly(self, var):
        """

	"""
        R = self.parent().base_ring()
        P = PolynomialRing(R, var)
        return P(self.ambient_algebra_element().charpoly('x'))

    characteristic_polynomial = charpoly

    def minpoly(self, var):
        """

	"""
        R = self.parent().base_ring()
        P = PolynomialRing(R, var)
        return P(self.ambient_algebra_element().minpoly())

    minimal_polynomial = minpoly

    def vector(self):
        """

	"""
        A = self.parent()
        R = A.base_ring()
        v = self.ambient_algebra_element().vector()
        M = A.inverse_embedding_matrix()
        V = FreeModule(R,4)
        return V([ R(x) for x in (v*M).list() ])

