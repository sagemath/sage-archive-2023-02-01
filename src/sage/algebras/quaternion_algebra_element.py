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
from sage.rings.polynomial.polynomial_ring import PolynomialRing
from sage.algebras.free_algebra_quotient_element import FreeAlgebraQuotientElement

class QuaternionAlgebraElement(FreeAlgebraQuotientElement):

    def __init__(self, H, x):
        """
        Create the element x of the quaternion algebra H.
        """
        FreeAlgebraQuotientElement.__init__(self, H, x)

    def __cmp__(self, x):
        return cmp(self.vector(), x.vector())

    def conjugate(self):
        return self.parent()(self.reduced_trace()) - self

    def reduced_trace(self):
        r"""
        Return the reduced trace of this element.

	\note{In a quaternion algebra $A$, every element $x$ is
	quadratic over the center, thus $x^2 = \Tr(x)*x - \Nr(x)$, so
	we solve for a linear relation $(1,-\Tr(x),\Nr(x))$ among
	$[x^2, x, 1]$ for the reduced trace of $x$.}
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
	"""
        x = self * self.conjugate()
	return x.vector()[0]

    def charpoly(self, var):
        """
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
	"""
	v = self.vector()
	if v[1] == 0 and v[2] == 0 and v[3] == 0:
	    K = self.parent().base_ring()
	    P = PolynomialRing(K, var)
	    x = P.gen()
	    return x - v[0]
	return self.charpoly(var)

    minimal_polynomial = minpoly
