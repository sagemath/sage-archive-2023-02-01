"""
Hypersurfaces in affine and projective space

AUTHOR: 2005-12-08, William Stein <wstein@ucsd.edu>
                    David Kohel <kohel@maths.usyd.edu.au>
"""

#*****************************************************************************
#  Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import is_MPolynomialRingElement
from algebraic_scheme import AlgebraicScheme
from affine_scheme import AffineScheme_generic
from affine_space import AffineSpace
from projective_scheme import ProjectiveScheme
from projective_space import ProjectiveSpace

def is_Hypersurface(self):
    return isinstance(self,Hypersurface)

class Hypersurface(AlgebraicScheme):
    """
    A hypersurface defined by a polynomial in affine or projective space
    """
    def __init__(self):
        pass

    def defining_polynomial(self):
        return self.defining_polynomials()[0]

class ProjectiveHypersurface(ProjectiveScheme, HypersurfaceScheme):
    """
    The projective hypersurface defined by the given polynomial.
    """
    def __init__(self, ambient, poly):
        if not is_MPolynomialRingElement(poly):
            raise TypeError, \
                  "Defining polynomial (=%s) must be a multivariate polynomial."%poly
        if not poly.is_homogeneous():
            raise TypeError, "Defining polynomial (=%s) must be homogeneous."%poly
        if ambient == None:
            P = poly.parent()
            ambient = ProjectiveSpace(P.base_ring(), P.ngens()-1, coordinate_ring = P)
        ProjectiveScheme.__init__(self, ambient, [poly])
        HypersurfaceScheme.__init__(self)

    def _repr_(self):
        return "Projective hypersurface defined by %s in %s"%(
            self.defining_polynomial(), self.ambient_space())

class AffineHypersurface(AffineScheme_generic, HypersurfaceScheme):
    """
    The affine hypersurface defined by the given polynomial.
    """
    def __init__(self, ambient, poly):
        if not is_MPolynomialRingElement(poly):
            raise TypeError, "Defining polynomial (= %s) must be a multivariate polynomial"%poly
        if ambient == None:
            P = poly.parent()
            ambient = AffineSpace(P.base_ring(), P.ngens(), coordinate_ring = P)
        AffineScheme.__init__(self, ambient, [poly])
        HypersurfaceScheme.__init__(self)

    def _repr_(self):
        return "Affine hypersurface defined by %s in %s"%(
            self.defining_polynomial(), self.ambient_space())




