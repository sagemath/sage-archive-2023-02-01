"""
Set of homomorphisms between two schemes
"""

#*****************************************************************************
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.structure.parent_old as parent_old
from sage.categories.homset import HomsetWithBase
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.constructor import is_FiniteField
from sage.rings.rational_field import is_RationalField
from sage.rings.morphism import is_RingHomomorphism

import spec
import morphism
from sage.schemes.generic.toric_morphism import SchemeMorphism_toric_coordinates_field




def is_SchemeHomset(H):
    return isinstance(H, SchemeHomset_generic)

def SchemeHomset(R, S, category=None, check=True):
    if spec.is_Spec(R) and spec.is_Spec(S):
        return SchemeHomset_spec(R, S, category=category, check=check)
    else:
        return SchemeHomset_generic(R, S, category=category, check=check)

class SchemeHomset_generic(parent_old.Parent, HomsetWithBase):
    def __init__(self, X, Y, category=None, check=True, base=ZZ):
        HomsetWithBase.__init__(self, X, Y, category=category, check=check, base=base)

    def has_coerce_map_from_impl(self, S):
        if self == S:   # an obvious case
            return True
        # Todo -- implement more cases.
        return False

    def _repr_(self):
        try:
            return "Set of points of %s defined over %s"%(self.codomain(), self.domain().coordinate_ring())
        except ValueError:
            return "Set of morphisms from %s to %s"%(self.domain(), self.codomain())

    def natural_map(self):
        X = self.domain()
        Y = self.codomain()
        if spec.is_Spec(Y) and Y.coordinate_ring() == X.base_ring():
            return morphism.SchemeMorphism_structure_map(self)
        raise NotImplementedError

    def __call__(self, x, check=True):
        """
        INPUT:

            - `x` -- a ring morphism, or a list or tuple of that
              define a ring morphism.

            - ``check`` -- (default: True) passed onto functions
              called by this to be more careful about input argument
              type checking

        EXAMPLES::

            sage: f = ZZ.hom(QQ); f
            Ring Coercion morphism:
              From: Integer Ring
              To:   Rational Field

        ::

            sage: H = Hom(Spec(QQ,ZZ), Spec(ZZ)); H
            Set of points of Spectrum of Integer Ring defined over Rational Field

        ::

            sage: phi = H(f); phi
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Ring Coercion morphism:
                      From: Integer Ring
                      To:   Rational Field

        TESTS::

        We illustrate input type checking::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: C = A.subscheme(x*y-1)
            sage: H = C.Hom(C); H
            Set of points of Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x*y - 1 defined over Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x*y - 1)
            sage: H(1)
            Traceback (most recent call last):
            ...
            TypeError: x must be a ring homomorphism, list or tuple
        """
        if isinstance(x, (list, tuple)):
            return self.domain()._point_morphism_class(self, x, check=check)

        if is_RingHomomorphism(x):
            return morphism.SchemeMorphism_spec(self, x, check=check)

        raise TypeError, "x must be a ring homomorphism, list or tuple"

class SchemeHomset_spec(SchemeHomset_generic):
    pass

class SchemeHomset_coordinates(SchemeHomset_generic):
    """
    Set of points on X defined over the base ring of X, and given by
    explicit tuples.
    """
    def __init__(self, X, S):
        R = X.base_ring()
        if R != S:
            X = X.base_extend(S)
        SchemeHomset_generic.__init__(self, spec.Spec(S, R), X)

    def _repr_(self):
        S = self.domain()
        if S == self.codomain().base_scheme():
            return "Set of Rational Points of %s"%self.codomain()
        if spec.is_Spec(S):
            S = S.coordinate_ring()
        return "Set of Rational Points over %s of %s"%(S, self.codomain())

    def value_ring(self):
        """
        Returns S for a homset X(T) where T = Spec(S).
        """
        T = self.domain()
        if spec.is_Spec(T):
            return T.coordinate_ring()
        else:
            raise TypeError, "Domain of argument must be of the form Spec(S)."


#*******************************************************************
# Affine varieties
#*******************************************************************
class SchemeHomset_affine_coordinates(SchemeHomset_coordinates):
    """
    Set of points on X defined over the base ring of X, and given by
    explicit tuples.
    """
    def __call__(self, *v):
        if len(v) == 1:
            v = v[0]
        return morphism.SchemeMorphism_affine_coordinates(self, v)

    def points(self, B=0):
        r"""
        Return the set of points given by coordinate tuples with coordinates
        in the base ring.

        INPUT:

        - ``B`` -- an integer.

        OUTPUT:

        - If the base ring is a finite field: the set of points given by
          coordinate tuples.

        - If the base ring is `\QQ` or `\ZZ`: the subset of points whose
          coordinates have height ``B`` or less.

        EXAMPLES: The bug reported at #11526 is fixed::

            sage: R = ZZ
            sage: A.<x,y> = R[]
            sage: I = A.ideal(x^2-y^2-1)
            sage: V = AffineSpace(R,2)
            sage: X = V.subscheme(I)
            sage: M = X(R)
            sage: M.points(1)
            [(-1, 0), (1, 0)]
        """
        try:
            R = self.value_ring()
        except TypeError:
            raise TypeError, "Domain of argument must be of the form Spec(S)."
        if is_RationalField(R) or R == ZZ:
            if not B > 0:
                raise TypeError, "A positive bound B (= %s) must be specified."%B
            from sage.schemes.generic.rational_point import enum_affine_rational_field
            return enum_affine_rational_field(self,B)
        elif is_FiniteField(R):
            from sage.schemes.generic.rational_point import enum_affine_finite_field
            return enum_affine_finite_field(self)
        else:
            raise TypeError, "Unable to enumerate points over %s."%R


#*******************************************************************
# Projective varieties
#*******************************************************************
class SchemeHomset_projective_coordinates_field(SchemeHomset_coordinates):
    """
    Set of points on X defined over the base ring of X, and given by
    explicit tuples.
    """
    def __call__(self, *v):
        if len(v) == 1:
            v = v[0]
        X = self.codomain()
        try:
            return X._point_class(X, v)
        except AttributeError:  # should be very rare
            return morphism.SchemeMorphism_projective_coordinates_field(self, v)

    def points(self, B=0):
        from sage.schemes.generic.rational_point import enum_projective_rational_field
        from sage.schemes.generic.rational_point import enum_projective_finite_field
        try:
            R = self.value_ring()
        except TypeError:
            raise TypeError, "Domain of argument must be of the form Spec(S)."
        if is_RationalField(R):
            if not B > 0:
                raise TypeError, "A positive bound B (= %s) must be specified."%B
            return enum_projective_rational_field(self,B)
        elif is_FiniteField(R):
            return enum_projective_finite_field(self)
        else:
            raise TypeError, "Unable to enumerate points over %s."%R

class SchemeHomset_projective_coordinates_ring(SchemeHomset_coordinates):
    """
    Set of points on X defined over the base ring of X, and given by
    explicit tuples.
    """
    def __call__(self, *v):
        raise NotImplementedError

    def points(self, B=0):
        raise NotImplementedError # fixed when __call__ is defined.
        try:
            R = self.value_ring()
        except TypeError:
            raise TypeError, "Domain of argument must be of the form Spec(S)."
        if R == ZZ:
            if not B > 0:
                raise TypeError, "A positive bound B (= %s) must be specified."%B
            from sage.schemes.generic.rational_points import enum_projective_rational_field
            return enum_projective_rational_field(self,B)
        else:
            raise TypeError, "Unable to enumerate points over %s."%R


#*******************************************************************
# Abelian varieties
#*******************************************************************
class SchemeHomsetModule_abelian_variety_coordinates_field(SchemeHomset_projective_coordinates_field):
    def __init__(self, X, S, category=None, check=True):
        r"""
        EXAMPLES: The bug reported at trac #1785 is fixed::

            sage: K.<a> = NumberField(x^2 + x - (3^3-3))
            sage: E = EllipticCurve('37a')
            sage: X = E(K)
            sage: X
            Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 + (-1)*x over Number Field in a with defining polynomial x^2 + x - 24
            sage: P = X([3,a])
            sage: P
            (3 : a : 1)
            sage: P in E
            False
            sage: P in E.base_extend(K)
            True
        """
        R = X.base_ring()
        if R != S:
            X = X.base_extend(S)
        Y = spec.Spec(S, R)
        HomsetWithBase.__init__(self, Y, X, category=category,
                                check = check, base = ZZ)

    def _repr_(self):
        return "Abelian group of points on %s"%self.codomain()

    def base_extend(self, R):
        if R != sage.rings.integer_ring.ZZ:
            raise NotImplementedError, "Abelian variety point sets not implemented as modules over rings other than ZZ."
        return self


#*******************************************************************
# Toric varieties
#*******************************************************************
class SchemeHomset_toric_coordinates_field(SchemeHomset_coordinates):
    """
    Construct the `Hom`-space of morphisms given on coordinates.

    .. WARNING::

        You should not create objects of this class directly.


    INPUT:

    - same as for :class:`SchemeHomset_coordinates`.

    OUPUT:

    - :class:`SchemeHomset_toric_coordinates_field`.

    TESTS::

        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: P1xP1 = ToricVariety(fan)
        sage: import sage.schemes.generic.homset as HOM
        sage: HOM.SchemeHomset_toric_coordinates_field(P1xP1, QQ)
        Set of Rational Points of 2-d toric variety
        covered by 4 affine patches

    A better way to construct the same `Hom`-space as above::

        sage: P1xP1(QQ)
        Set of Rational Points of 2-d toric variety
        covered by 4 affine patches
    """
    # Mimicking SchemeHomset_projective_coordinates_field,
    # affine spaces implement only "except" case
    def __call__(self, *arg):
        r"""
        Construct a morphism from given parameters.

        INPUT:

        - data determining a morphism.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: import sage.schemes.generic.homset as HOM
            sage: H = HOM.SchemeHomset_toric_coordinates_field(P1xP1, QQ)
            sage: H(1,2,3,4)
            [1 : 2 : 3 : 4]
        """
        # This may break for one-dimensional varieties.
        if len(arg) == 1:
            arg = arg[0]
        X = self.codomain()
        try:
            return X._point_class(X, arg)
        except AttributeError:  # should be very rare
            return SchemeMorphism_toric_coordinates_field(self, arg)

