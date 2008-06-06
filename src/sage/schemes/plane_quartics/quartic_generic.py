"""
Plane quartic curves over a general ring.  These are generic genus 3 curves,
as distinct from hyperelliptic curves of genus 3.

EXAMPLE:
    sage: PP, (X,Y,Z) = ProjectiveSpace(2,QQ,'X,Y,Z').objgens()
    sage: f = X^4 + Y^4 + Z^4 - 3*X*Y*Z*(X+Y+Z)
    sage: C = QuarticCurve(f); C
    Quartic Curve over Rational Field defined by X^4 + Y^4 - 3*X^2*Y*Z - 3*X*Y^2*Z - 3*X*Y*Z^2 + Z^4
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import PolynomialRing
from sage.misc.all import latex

import sage.schemes.plane_curves.projective_curve as projective_curve
from sage.schemes.generic.all import ProjectiveSpace

def is_QuarticCurve(C):
    return isinstance(C,QuarticCurve_generic)

class QuarticCurve_generic(projective_curve.ProjectiveCurve_generic):
    # DRK: Note that we should check whether the curve is

    def _repr_type(self):
        return "Quartic"

    def genus(self):
        return 3
