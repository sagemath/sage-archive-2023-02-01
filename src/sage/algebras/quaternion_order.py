"""
Quaternion orders

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

from sage.rings.ring import Ring
from sage.matrix.matrix_space import MatrixSpace
from sage.algebras.free_algebra_quotient import FreeAlgebraQuotient
from sage.algebras.algebra_order import AlgebraOrder_generic
from sage.algebras.quaternion_algebra import QuaternionAlgebra_generic
from sage.algebras.quaternion_order_element import QuaternionOrderElement

def QuaternionOrderWithBasis(R, B):
    if not isinstance(R, Ring):
        raise TypeError, "Argument R (= %s) must be a ring."%(R)
    if not isinstance(B, (list,tuple)):
        raise TypeError, "Argument B (= %s) must be a list."%B
    if len(B) != 4:
        raise TypeError, "Argument B (= %s) must have length 4.."%B
    if not B[0] == 1:
        raise TypeError, "Argument B (= %s) must begin with 1."%B
    H = B[0].parent()
    if not isinstance(H, QuaternionAlgebra_generic):
        raise TypeError, "Argument B (= %s) must be a sequence of quaternions."%B
    return QuaternionOrder_generic(H, R, gens=B[1:], basis=B)

def QuaternionDefiningOrder(H, R):
    """
    Returns a hypothetical underlying order of H spanned by the basis
    over R.

    No checking is done to ensure that the ring R is an integral domain
    whose field of fractions is the base field of H. Nor is it checked
    that the basis of H in fact generates an R-module which is closed
    under multiplication.
    """
    if not isinstance(H, QuaternionAlgebra_generic):
        raise TypeError, "Argument H (= %s) must be a quaternion algebra."%H
    if not isinstance(R, Ring):
        raise TypeError, "Argument R (= %s) must be a ring."%R
    B = H.basis()
    return QuaternionOrder_generic(H, R, gens=B[1:], basis=B)

class QuaternionOrder_generic(AlgebraOrder_generic):
    """
    An order in a quaternion algebra.
    """
    def __init__(self, H, R, gens, basis=None):
        """
	An order in a quaternion algebra.
	"""
        AlgebraOrder_generic.__init__(self, H, R, gens, basis=basis, rank=4)

    def __call__(self, x, check=True):
        if isinstance(x, QuaternionOrderElement) and x.parent() is self:
	    return x
        return QuaternionOrderElement(self, x, check=check)

    def __repr__(self):
        return "Quaternion order over %s with basis %s"%(
            self.base_ring(), self.basis())

    def discriminant(self):
        return self.gram_matrix().determinant()

    def gram_matrix(self):
        L = self.module()
        if not L._FreeModule_generic__inner_product_is_dot_product:
            return L.gram_matrix()
        R = self.base_ring()
        M = MatrixSpace(R,4,4)(0)
        B = self.basis()
        for i in range(4):
            x = B[i]
            M[i,i] = 2*(x.reduced_norm())
            for j in range(i+1,4):
                y = B[j]
                c = (x * y.conjugate()).reduced_trace()
                M[i,j] = c
                M[j,i] = c
        return M

    def inner_product_matrix(self):
        return self.gram_matrix()

    def random_element(self):
        R = self.base_ring()
        return self([ R.random_element() for _ in range(4) ])
