"""
Curve points.

EXAMPLES:

We can create points on projective curves::

    sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
    sage: C = Curve([x^3 - 2*x*z^2 - y^3, z^3 - w^3 - x*y*z], P)
    sage: Q = C([1,1,0,0])
    sage: type(Q)
    <class 'sage.schemes.curves.point.ProjectiveCurvePoint_field'>
    sage: Q.parent()
    Set of rational points of Projective Curve over Rational Field defined
    by x^3 - y^3 - 2*x*z^2, -x*y*z + z^3 - w^3

or on affine curves::

    sage: A.<x,y> = AffineSpace(GF(23), 2)
    sage: C = Curve([y - y^4 + 17*x^2 - 2*x + 22], A)
    sage: Q = C([22,21])
    sage: type(Q)
    <class 'sage.schemes.curves.point.AffinePlaneCurvePoint_finite_field'>
    sage: Q.parent()
    Set of rational points of Affine Plane Curve over Finite Field of size
    23 defined by -y^4 - 6*x^2 - 2*x + y - 1

regardless of whether the curve is a plane curve or not.

AUTHORS:

- Grayson Jorgenson (2016-6)
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.schemes.affine.affine_point import (SchemeMorphism_point_affine_field,
                                              SchemeMorphism_point_affine_finite_field)
from sage.schemes.projective.projective_point import (SchemeMorphism_point_projective_field,
                                                      SchemeMorphism_point_projective_finite_field)

class ProjectiveCurvePoint_field(SchemeMorphism_point_projective_field):

    def curve(self):
        r"""
        Return the curve that this point is on.

        OUTPUT:

        - the projective curve that is the codomain of this point.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(GF(17), 3)
            sage: C = Curve([x^3 - y^2*z, z^2 - x*y])
            sage: Q = C([0,16,0,1])
            sage: Q.curve()
            Projective Curve over Finite Field of size 17 defined by x^3 - y^2*z, -x*y + z^2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y^3 - x^3 - x^2*z])
            sage: Q = C([-1,0,1])
            sage: Q.curve()
            Projective Plane Curve over Rational Field defined by -x^3 + y^3 - x^2*z
        """
        return self.codomain()

    def multiplicity(self):
        r"""
        Return the multiplicity of this point with respect to the projective curve it is on.

        OUTPUT: Integer.

        EXAMPLES::

            sage: P.<x,y,z,w,t> = ProjectiveSpace(GF(7), 4)
            sage: C = Curve([y^3 - 2*x^3 - z^3, x^3 - w^3, t - z])
            sage: Q = C([0,2,1,0,1])
            sage: Q.multiplicity()
            3
        """
        return self.curve().multiplicity(self)

    def is_singular(self):
        r"""
        Return whether this point is or is not a singular point of the projective curve it is on.

        OUTPUT: Boolean.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: C = Curve([x^2 - y^2, z - w], P)
            sage: Q1 = C([0,0,1,1])
            sage: Q1.is_singular()
            True
            sage: Q2 = C([1,1,1,1])
            sage: Q2.is_singular()
            False
        """
        return self.curve().is_singular(self)

class ProjectivePlaneCurvePoint_field(ProjectiveCurvePoint_field):

    def tangents(self):
        r"""
        Return the tangents at this point of the projective plane curve this point is on.

        OUTPUT:

        - a list of polynomials in the coordinate ring of the ambient space of the curve this point is on.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y^2*z^3 - x^5 + 18*y*x*z^3])
            sage: Q = C([0,0,1])
            sage: Q.tangents()
            [y, 18*x + y]
        """
        return self.curve().tangents(self)

    def is_ordinary_singularity(self):
        r"""
        Return whether this point is an ordinary singularity of the projective plane curve it is on.

        OUTPUT: Boolean.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([z^6 - x^6 - x^3*z^3 - x^3*y^3])
            sage: Q = C([0,1,0])
            sage: Q.is_ordinary_singularity()
            False

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^2 - 3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: C = P.curve([x^2*y^3*z^4 - y^6*z^3 - 4*x^2*y^4*z^3 - 4*x^4*y^2*z^3 + 3*y^7*z^2 + 10*x^2*y^5*z^2\
            + 9*x^4*y^3*z^2 + 5*x^6*y*z^2 - 3*y^8*z - 9*x^2*y^6*z - 11*x^4*y^4*z - 7*x^6*y^2*z - 2*x^8*z + y^9\
            + 2*x^2*y^7 + 3*x^4*y^5 + 4*x^6*y^3 + 2*x^8*y])
            sage: Q = C([-1/2, 1/2, 1])
            sage: Q.is_ordinary_singularity()
            True
        """
        return self.curve().is_ordinary_singularity(self)

class ProjectivePlaneCurvePoint_finite_field(ProjectivePlaneCurvePoint_field, SchemeMorphism_point_projective_finite_field):
    pass

class AffineCurvePoint_field(SchemeMorphism_point_affine_field):

    def curve(self):
        r"""
        Return the curve that this point is on.

        OUTPUT:

        - the affine curve that is the codomain of this point.

        EXAMPLES::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: C = A.curve([y - x^5 - z^2, x - w^3 - z*y, y - x - w])
            sage: Q = C([0,1,-1,1])
            sage: Q.curve()
            Affine Curve over Rational Field defined by -x^5 - z^2 + y, -w^3 - y*z + x, -x + y - w

        ::

            sage: A.<x,y> = AffineSpace(GF(11), 2)
            sage: C = Curve([y - 6*x^4 - y^2], A)
            sage: Q = C([8,7])
            sage: Q.curve()
            Affine Plane Curve over Finite Field of size 11 defined by 5*x^4 - y^2 + y
        """
        return self.codomain()

    def multiplicity(self):
        r"""
        Return the multiplicity of this point with respect to the affine curve it is on.

        OUTPUT: Integer.

        EXAMPLES::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: C = Curve([y^4 - 17*x^2 - x^3 + z^3, z^3 - x^2])
            sage: Q = C([0,0,0])
            sage: Q.multiplicity()
            6
        """
        return self.curve().multiplicity(self)

    def is_singular(self):
        r"""
        Return whether this point is or is not a singular point of the affine curve it is on.

        OUTPUT: Boolean.

        EXAMPLES::

            sage: K = QuadraticField(-1)
            sage: A.<x,y,z> = AffineSpace(K, 3)
            sage: C = Curve([(x^4 + 2*z + 2)*y, z - y + 1])
            sage: Q1 = C([0,0,-1])
            sage: Q1.is_singular()
            True
            sage: Q2 = C([-K.gen(),0,-1])
            sage: Q2.is_singular()
            False
        """
        return self.curve().is_singular(self)

class AffinePlaneCurvePoint_field(AffineCurvePoint_field):

    def tangents(self):
        r"""
        Return the tangents at this point of the affine plane curve this point is on.

        OUTPUT:

        - a list of polynomials in the coordinate ring of the ambient space of the curve this point is on.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve([x^5 - x^3*y^2 + 5*x^4 - x^3*y - 3*x^2*y^2 + x*y^3 + 10*x^3 - 3*x^2*y -\
            3*x*y^2 + y^3 + 10*x^2 - 3*x*y - y^2 + 5*x - y + 1])
            sage: Q = C([-1,0])
            sage: Q.tangents()
            [y, -x + y - 1, x + 1, x + y + 1]
        """
        return self.curve().tangents(self)

    def is_ordinary_singularity(self):
        r"""
        Return whether this point is an ordinary singularity of the affine plane curve it is on.

        OUTPUT: Boolean.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = A.curve([x^5 - x^3*y^2 + 5*x^4 - x^3*y - 3*x^2*y^2 + x*y^3 + 10*x^3 - 3*x^2*y -\
            3*x*y^2 + y^3 + 10*x^2 - 3*x*y - y^2 + 5*x - y + 1])
            sage: Q = C([-1,0])
            sage: Q.is_ordinary_singularity()
            True

        ::

            sage: A.<x,y> = AffineSpace(GF(7), 2)
            sage: C = A.curve([y^2 - x^7 - 6*x^3])
            sage: Q = C([0,0])
            sage: Q.is_ordinary_singularity()
            False
        """
        return self.curve().is_ordinary_singularity(self)

class AffinePlaneCurvePoint_finite_field(AffinePlaneCurvePoint_field, SchemeMorphism_point_affine_finite_field):
    pass
