r"""
Set of homomorphisms between two proejctive schemes

For schemes `X` and `Y`, this module implements the set of morphisms
`Hom(X,Y)`. This is done by :class:`SchemeHomset_generic`.

As a special case, the Hom-sets can also represent the points of a
scheme. Recall that the `K`-rational points of a scheme `X` over `k`
can be identified with the set of morphisms `Spec(K) \to X`. In Sage
the rational points are implemented by such scheme morphisms. This is
done by :class:`SchemeHomset_points` and its subclasses.

.. note::

    You should not create the Hom-sets manually. Instead, use the
    :meth:`~sage.structure.parent.Hom` method that is inherited by all
    schemes.

AUTHORS:

- William Stein (2006): initial version.

- Volker Braun (2011-08-11): significant improvement and refactoring.

- Ben Hutz (June 2012): added support for projective ring
"""


#*****************************************************************************
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import ZZ
from sage.schemes.generic.homset import SchemeHomset_points

from sage.rings.rational_field import is_RationalField
from sage.rings.finite_rings.constructor import is_FiniteField

#*******************************************************************
# Projective varieties
#*******************************************************************
class SchemeHomset_points_projective_field(SchemeHomset_points):
    """
    Set of rational points of a projective variety over a field.

    INPUT:

    See :class:`SchemeHomset_generic`.

    EXAMPLES::

        sage: from sage.schemes.projective.projective_homset import SchemeHomset_points_projective_field
        sage: SchemeHomset_points_projective_field(Spec(QQ), ProjectiveSpace(QQ,2))
        Set of rational points of Projective Space of dimension 2 over Rational Field
    """
    def points(self, B=0):
        """
        Return some or all rational points of a projective scheme.

        INPUT:

        - `B` -- integer (optional, default=0). The bound for the
          coordinates.

        OUTPUT:

        A list of points. Over a finite field, all points are
        returned. Over an infinite field, all points satisfying the
        bound are returned.

        EXAMPLES::

            sage: P1 = ProjectiveSpace(GF(2),1)
            sage: F.<a> = GF(4,'a')
            sage: P1(F).points()
            [(0 : 1), (1 : 0), (1 : 1), (a : 1), (a + 1 : 1)]
        """
        from sage.schemes.projective.projective_rational_point import enum_projective_rational_field
        from sage.schemes.projective.projective_rational_point import enum_projective_finite_field
        R = self.value_ring()
        if is_RationalField(R):
            if not B > 0:
                raise TypeError, "A positive bound B (= %s) must be specified."%B
            return enum_projective_rational_field(self,B)
        elif is_FiniteField(R):
            return enum_projective_finite_field(self.extended_codomain())
        else:
            raise TypeError, "Unable to enumerate points over %s."%R

class SchemeHomset_points_projective_ring(SchemeHomset_points):
    """
    Set of rational points of a projective variety over a commutative ring.

    INPUT:

    See :class:`SchemeHomset_generic`.

    EXAMPLES::

        sage: from sage.schemes.projective.projective_homset import SchemeHomset_points_projective_ring
        sage: SchemeHomset_points_projective_ring(Spec(ZZ), ProjectiveSpace(ZZ,2))
        Set of rational points of Projective Space of dimension 2 over Integer Ring
    """

    def points(self, B=0):
        """
        Return some or all rational points of a projective scheme.

        INPUT:

        - `B` -- integer (optional, default=0). The bound for the
          coordinates.

        EXAMPLES::

            sage: from sage.schemes.projective.projective_homset import SchemeHomset_points_projective_ring
            sage: H = SchemeHomset_points_projective_ring(Spec(ZZ), ProjectiveSpace(ZZ,2))
            sage: H.points(3)
            [(0 : 0 : 1), (0 : 1 : -3), (0 : 1 : -2), (0 : 1 : -1), (0 : 1 : 0), (0
            : 1 : 1), (0 : 1 : 2), (0 : 1 : 3), (0 : 2 : -3), (0 : 2 : -1), (0 : 2 :
            1), (0 : 2 : 3), (0 : 3 : -2), (0 : 3 : -1), (0 : 3 : 1), (0 : 3 : 2),
            (1 : -3 : -3), (1 : -3 : -2), (1 : -3 : -1), (1 : -3 : 0), (1 : -3 : 1),
            (1 : -3 : 2), (1 : -3 : 3), (1 : -2 : -3), (1 : -2 : -2), (1 : -2 : -1),
            (1 : -2 : 0), (1 : -2 : 1), (1 : -2 : 2), (1 : -2 : 3), (1 : -1 : -3),
            (1 : -1 : -2), (1 : -1 : -1), (1 : -1 : 0), (1 : -1 : 1), (1 : -1 : 2),
            (1 : -1 : 3), (1 : 0 : -3), (1 : 0 : -2), (1 : 0 : -1), (1 : 0 : 0), (1
            : 0 : 1), (1 : 0 : 2), (1 : 0 : 3), (1 : 1 : -3), (1 : 1 : -2), (1 : 1 :
            -1), (1 : 1 : 0), (1 : 1 : 1), (1 : 1 : 2), (1 : 1 : 3), (1 : 2 : -3),
            (1 : 2 : -2), (1 : 2 : -1), (1 : 2 : 0), (1 : 2 : 1), (1 : 2 : 2), (1 :
            2 : 3), (1 : 3 : -3), (1 : 3 : -2), (1 : 3 : -1), (1 : 3 : 0), (1 : 3 :
            1), (1 : 3 : 2), (1 : 3 : 3), (2 : -3 : -3), (2 : -3 : -2), (2 : -3 :
            -1), (2 : -3 : 0), (2 : -3 : 1), (2 : -3 : 2), (2 : -3 : 3), (2 : -2 :
            -3), (2 : -2 : -1), (2 : -2 : 1), (2 : -2 : 3), (2 : -1 : -3), (2 : -1 :
            -2), (2 : -1 : -1), (2 : -1 : 0), (2 : -1 : 1), (2 : -1 : 2), (2 : -1 :
            3), (2 : 0 : -3), (2 : 0 : -1), (2 : 0 : 1), (2 : 0 : 3), (2 : 1 : -3),
            (2 : 1 : -2), (2 : 1 : -1), (2 : 1 : 0), (2 : 1 : 1), (2 : 1 : 2), (2 :
            1 : 3), (2 : 2 : -3), (2 : 2 : -1), (2 : 2 : 1), (2 : 2 : 3), (2 : 3 :
            -3), (2 : 3 : -2), (2 : 3 : -1), (2 : 3 : 0), (2 : 3 : 1), (2 : 3 : 2),
            (2 : 3 : 3), (3 : -3 : -2), (3 : -3 : -1), (3 : -3 : 1), (3 : -3 : 2),
            (3 : -2 : -3), (3 : -2 : -2), (3 : -2 : -1), (3 : -2 : 0), (3 : -2 : 1),
            (3 : -2 : 2), (3 : -2 : 3), (3 : -1 : -3), (3 : -1 : -2), (3 : -1 : -1),
            (3 : -1 : 0), (3 : -1 : 1), (3 : -1 : 2), (3 : -1 : 3), (3 : 0 : -2), (3
            : 0 : -1), (3 : 0 : 1), (3 : 0 : 2), (3 : 1 : -3), (3 : 1 : -2), (3 : 1
            : -1), (3 : 1 : 0), (3 : 1 : 1), (3 : 1 : 2), (3 : 1 : 3), (3 : 2 : -3),
            (3 : 2 : -2), (3 : 2 : -1), (3 : 2 : 0), (3 : 2 : 1), (3 : 2 : 2), (3 :
            2 : 3), (3 : 3 : -2), (3 : 3 : -1), (3 : 3 : 1), (3 : 3 : 2)]
        """
        R = self.value_ring()
        if R == ZZ:
            if not B > 0:
                raise TypeError, "A positive bound B (= %s) must be specified."%B
            from sage.schemes.projective.projective_rational_point import enum_projective_rational_field
            return enum_projective_rational_field(self,B)
        else:
            raise TypeError, "Unable to enumerate points over %s."%R


#*******************************************************************
# Abelian varieties
#*******************************************************************
class SchemeHomset_points_abelian_variety_field(SchemeHomset_points_projective_field):
    r"""
    Set of rational points of an abelian variety.

    INPUT:

    See :class:`SchemeHomset_generic`.

    TESTS:

    The bug reported at trac #1785 is fixed::

        sage: K.<a> = NumberField(x^2 + x - (3^3-3))
        sage: E = EllipticCurve('37a')
        sage: X = E(K)
        sage: X
        Abelian group of points on Elliptic Curve defined by
        y^2 + y = x^3 + (-1)*x over Number Field in a with
        defining polynomial x^2 + x - 24
        sage: P = X([3,a])
        sage: P
        (3 : a : 1)
        sage: P in E
        False
        sage: P in E.base_extend(K)
        True
        sage: P in X.codomain()
        False
        sage: P in X.extended_codomain()
        True
    """

    def _element_constructor_(self, *v, **kwds):
        """
        The element contstructor.

        INPUT:

        - ``v`` -- anything that determines a scheme morphism in the
          Hom-set.

        OUTPUT:

        The scheme morphism determined by ``v``.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: X = E(QQ)
            sage: P = X([0,1,0]);  P
            (0 : 1 : 0)
            sage: type(P)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_number_field'>

        TESTS::

            sage: X._element_constructor_([0,1,0])
            (0 : 1 : 0)
        """
        if len(v) == 1:
            v = v[0]
        return self.codomain()._point(self.extended_codomain(), v, **kwds)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: X = E(QQ)
            sage: X._repr_()
            'Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field'
        """
        s = 'Abelian group of points on ' + str(self.extended_codomain())
        return s

    def base_extend(self, R):
        """
        Extend the base ring.

        This is currently not implemented except for the trivial case
        ``R==ZZ``.

        INPUT:

        - ``R`` -- a ring.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: Hom = E.point_homset();  Hom
            Abelian group of points on Elliptic Curve defined
            by y^2 + y = x^3 - x over Rational Field
            sage: Hom.base_ring()
            Integer Ring
            sage: Hom.base_extend(QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError: Abelian variety point sets are not
            implemented as modules over rings other than ZZ.
        """
        if R is not ZZ:
            raise NotImplementedError('Abelian variety point sets are not '
                            'implemented as modules over rings other than ZZ.')
        return self


from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.schemes.generic.homset',
                           'SchemeHomsetModule_abelian_variety_coordinates_field',
                           SchemeHomset_points_abelian_variety_field)

