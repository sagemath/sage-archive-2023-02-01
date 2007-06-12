"""nodoctest
Jacobian of a General Hyperelliptic Curve
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import is_Ring
from sage.schemes.jacobians.abstract_jacobian import Jacobian_generic
import sage.schemes.generic.homset as homset
import sage.schemes.generic.morphism as morphism
from hyperelliptic_generic import is_HyperellipticCurve
import jacobian_homset
import jacobian_morphism

class HyperellipticJacobian_generic(Jacobian_generic):
    def __call__(self, *args, **kwds):
        """
        EXAMPLES:
	    sage: FF = FiniteField(2003)
	    sage: R.<x> = PolynomialRing(FF)
	    sage: f = x**5 + 1184*x**3 + 1846*x**2 + 956*x + 560
	    sage: C = HyperellipticCurve(f)
            sage: J = Jacobian(C)
	    sage: a = x**2 + 376*x + 245; b = 1015*x + 1368
	    sage: X = J(FF)
	    sage: D = X([a,b])
	    sage: D
	    (x^2 + 376*x + 245, y + 988*x + 635)
	    sage: J(0)
	    (1)
	    sage: D == J([a,b])
	    True
	    sage: D == D + J(0)
	    True
  	"""
        if len(args) == 1:
            if is_Ring(args[0]):
                return jacobian_homset.JacobianHomset_divisor_classes(self, args[0])
	    return jacobian_homset.JacobianHomset_divisor_classes(self, self.base_ring())(args[0])
        raise TypeError, "Arguments must be a coefficient ring or Mumford divisor."

    def dimension(self):
        return self.__curve.genus()

    def _homset_class(self, *args, **kwds):
        return jacobian_homset.JacobianHomset_divisor_classes(*args, **kwds)

    def _point_class(self, *args, **kwds):
        return jacobian_morphism.JacobianMorphism_divisor_class_field(*args, **kwds)

    def _cmp_(self,other):
        if self.curve() == other.curve():
            return 0
        else:
            return -1

