"""
Set of homomorphisms between two schemes
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.all import Homset, Schemes
from sage.rings.all      import is_RingHomomorphism

import spec

import morphism

SCH = Schemes()

def is_SchemeHomset(H):
    return isinstance(H, SchemeHomset_generic)

def SchemeHomset(R, S, cat=None, check=True):
    if spec.is_Spec(R) and spec.is_Spec(S):
        return SchemeHomset_spec(R, S, cat=cat, check=check)
    else:
        return SchemeHomset_generic(R, S, cat=cat, check=check)

class SchemeHomset_generic(Homset):
    def __init__(self, X, Y, cat=None, check=True):
        Homset.__init__(self, X, Y, cat=cat, check=check)

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
        EXAMPLES:
            sage: f = Z.hom(Q); f
            Coercion morphism:
              From: Integer Ring
              To:   Rational Field

            sage: H = Hom(Spec(Q), Spec(Z)); H
            Set of points of Spectrum of Integer Ring defined over Rational Field

            sage: phi = H(f); phi
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Coercion morphism:
                      From: Integer Ring
                      To:   Rational Field
        """
        if isinstance(x, (list, tuple)):
            return self.codomain()._point_morphism_class(self, x, check=check)

        if is_RingHomomorphism(x):
            return morphism.SchemeMorphism_spec(self, x, check=check)

class SchemeHomset_spec(SchemeHomset_generic):
    pass

class SchemeHomset_coordinates(SchemeHomset_generic):
    """
    Set of points on X defined over the base ring of X, and given
    by explicit tuples.
    """
    def __init__(self, X):
        R = X.base_ring()
        SchemeHomset_generic.__init__(self, spec.Spec(R, R), X)

    def _repr_(self):
        return "Set of Rational Points of %s"%self.codomain()

class SchemeHomset_affine_coordinates(SchemeHomset_coordinates):
    """
    Set of points on X defined over the base ring of X, and given
    by explicit tuples.
    """
    def __call__(self, *v):
        if len(v) == 1:
            v = v[0]
        return morphism.SchemeMorphism_affine_coordinates(self.codomain(), v)

class SchemeHomset_projective_coordinates_field(SchemeHomset_coordinates):
    """
    Set of points on X defined over the base ring of X, and given
    by explicit tuples.
    """
    def __call__(self, *v):
        if len(v) == 1:
            v = v[0]
        return morphism.SchemeMorphism_projective_coordinates_field(self.codomain(), v)

class SchemeHomset_projective_coordinates_ring(SchemeHomset_coordinates):
    """
    Set of points on X defined over the base ring of X, and given
    by explicit tuples.
    """
    def __call__(self, *v):
        raise NotImplementedError
