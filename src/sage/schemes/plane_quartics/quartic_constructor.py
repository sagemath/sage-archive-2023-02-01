"""
Quartic curve constructor
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.schemes.generic.all import is_ProjectiveSpace, ProjectiveSpace
from sage.rings.all import is_MPolynomial

from quartic_generic import QuarticCurve_generic

#from sage.rings.all import is_FiniteField, is_RationalField

def QuarticCurve(F,PP=None,check=False):
    """
    Returns the quartic curve defined by F.
    """
    if not PP is None:
        if not is_ProjectiveSpace(PP) and PP.dimension == 2:
            raise TypeError, "Argument PP (=%s) must be a projective plane"%PP
    elif not is_MPolynomial(F):
        raise TypeError, \
              "Argument F (=%s) must be a homogeneous multivariate polynomial"%F
    else:
        if not F.is_homogeneous():
            "Argument F (=%s) must be a homogeneous polynomial"%F
        P = F.parent()
        if P.ngens() != 3:
            "Argument F (=%s) must be a homogeneous multivariate polynomial in 3 variables"%F
        PP = ProjectiveSpace(P)
    if check:
        raise TypeError, "Argument checking (for nonsingularity) is not implemented."
    return QuarticCurve_generic(PP, F)
