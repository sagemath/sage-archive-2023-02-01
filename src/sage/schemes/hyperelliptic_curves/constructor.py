"""
Hyperelliptic curves over a general ring

Author: 2005-11-13, David Joyner, (wdj@usna.edu)
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


from sage.schemes.generic.all import ProjectiveSpace
import hyperell_generic, hyperell_finite_field
from sage.rings.all import MPolynomial

def HyperellipticCurve(f):
    if not isinstance(f, MPolynomial):
        raise TypeError, "f (=%s) must be a multi-variate polynomial"%f
    R = f.parent()
    if R.ngens() != 3:
        raise TypeError, "f (=%s) must be a 3-variable polynomial"%f
    A = ProjectiveSpace(2, R.base_ring())
    A._coordinate_ring = R
    if R.base_ring().is_finite() and R.base_ring().is_field():
        return hyperell_finite_field.HyperellipticCurve_finite_field(A, f)
    return hyperell_generic.HyperellipticCurve_generic(A, f)

