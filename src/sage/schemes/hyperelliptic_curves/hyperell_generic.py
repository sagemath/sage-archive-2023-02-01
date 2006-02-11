"""
Hyperelliptic curves over a general ring

AUTHOR: 2005-11-13, David Joyner <wdj@usna.edu>
AUTHOR: 2005-11-13, William Stein <wstein@ucsd.edu>

EXAMPLE:
    sage: x, y, z = MPolynomialRing(GF(5), 3, 'xyz').gens()
    sage: HyperellipticCurve(y^2*z^7 - x^9 - x*z^8)
    Hyperelliptic Curve over Finite Field of size 5 defined by y^2*z^7 + 4*x*z^8 + 4*x^9

EXAMPLE:
    sage: F = GF(5)
    sage: x,y,z = (F['x,y,z']).gens()
    sage: f = y^2*z^7 - x^9 - x*z^8
    sage: C = HyperellipticCurve(f); C
    Hyperelliptic Curve over Finite Field of size 5 defined by y^2*z^7 + 4*x*z^8 + 4*x^9
    sage: C.genus()
    4

"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#                     2005 David Joyner <wdj@usna.edu>
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

from sage.rings.all import (factor, FractionField, MPolynomial, MPolynomialRing)
from sage.misc.all import latex

import sage.schemes.plane_curves.projective_curve as plane_curve
from sage.schemes.generic.all import ProjectiveSpace

def is_affine_hyperelliptic_curve_equation(f):
    if not isinstance(f, MPolynomial):
        return False
    if f.parent().ngens() != 2:
        return False
    if f.degree(f.parent().gen(1)) != 2:
        return False
    return True

class HyperellipticCurve_generic(plane_curve.ProjectiveCurve_generic):
    def _repr_type(self):
        """
        String representation of this sort of curve is returned by this function.

        EXAMPLES:
            sage: F = GF(5)
            sage: x,y,z = (F['x,y,z']).gens()
            sage: f = y^2*z^7 - x^9 - x*z^8
            sage: C = HyperellipticCurve(f); C
            Hyperelliptic Curve over Finite Field of size 5 defined by y^2*z^7 + 4*x*z^8 + 4*x^9
        """
        return "Hyperelliptic"

