"""
Jacobian of a hyperelliptic curve
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.schemes.jacobians.abstract_jacobian import Jacobian
import sage.schemes.generic.homset as homset
import sage.schemes.generic.morphism as morphism
from hyperelliptic_generic import is_HyperellipticCurve
import jacobian_homset

class HyperellipticJacobian_generic(Jacobian):
    def dimension(self):
        return self.__curve.genus()

    def _homset_class(self, *args, **kwds):
        return jacobian_homset.JacobianHomset_divisor_classes(*args, **kwds)

    def _point_morphism_class(self, *args, **kwds):
        return jacobian_morphism.JacobianMorphism_divisor_class(*args, **kwds)

    def _cmp_(self,other):
        if self.curve() == other.curve():
            return 0
        else:
            return -1

