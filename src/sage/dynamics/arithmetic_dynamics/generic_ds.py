r"""
Dynamical Systems of Schemes

AUTHORS:

- Ben Hutz (July 2017): initial version
"""

#*****************************************************************************
#       Copyright (C) 2017 Ben Hutz <bn4941@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function
from sage.categories.homset import End
from sage.schemes.generic.morphism import SchemeMorphism_polynomial
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField
from sage.rings.fraction_field import is_FractionField
from sage.schemes.affine.affine_space import is_AffineSpace
from sage.schemes.affine.affine_space import AffineSpace
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_affine
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_projective
from sage.schemes.projective.projective_space import is_ProjectiveSpace
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.schemes.product_projective.space import is_ProductProjectiveSpaces
from sage.symbolic.ring import SR

from sage.categories.fields import Fields
_Fields = Fields()

def is_DynamicalSystem(f):
    """
    Test whether ``f`` is a dynamical system.

    INPUT:

    - ``f`` -- anything.

    OUTPUT:

    Boolean. Return ``True`` if ``f`` is a dynamical system

    EXAMPLES::

        sage: A.<x,y> = AffineSpace(QQ,2)
        sage: f = DynamicalSystem_affine([y,x^2+y]); f
        Dynamical System of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (y, x^2 + y)
        sage: from sage.dynamics.arithmetic_dynamics.generic_ds import is_DynamicalSystem
        sage: is_DynamicalSystem(f)
        True
    """
    return isinstance(f, DynamicalSystem_generic)

def DynamicalSystem_affine(polys, domain=None, check=True):
    if domain is None:
        if isinstance(polys, SchemeMorphism_polynomial):
            domain = polys.domain()
            R = polys.base_ring()
        else:
            R = polys[0].base_ring()
            CR = polys[0].parent()
            if CR is SR:
                raise TypeError("cannot have Symbolic Ring as base ring")
            if is_FractionField(CR):
                domain = AffineSpace(CR.ring())
            else:
                domain = AffineSpace(CR)
    else:
        R = domain.base_ring()
        CR = domain.coordinate_ring()
    if is_AffineSpace(domain) or isinstance(domain, AlgebraicScheme_subscheme_affine):
        from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine_ring
        from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine_field
        from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine_finite_field
        if R in _Fields:
            if is_FiniteField(R):
                return DynamicalSystem_affine_finite_field(polys, domain=domain, check=check)
            else:
                return DynamicalSystem_affine_field(polys, domain=domain, check=check)
        return DynamicalSystem_affine_ring(polys, domain=domain, check=check)
    else:
        raise TypeError("must be affine or affine subscheme")

def DynamicalSystem_projective(polys, domain=None, check=True):
    if domain is None:
        if isinstance(polys, SchemeMorphism_polynomial):
            domain = polys.domain()
            R = polys.base_ring()
        else:
            CR = polys[0].parent()
            if CR is SR:
                raise TypeError("cannot have Symbolic Ring as base ring")
            R = CR.base_ring()
            domain = ProjectiveSpace(CR)
    elif is_ProductProjectiveSpaces(domain):
        from sage.dynamics.arithmetic_dynamics.product_projective_ds import DynamicalSystem_product_projective_ring
        return DynamicalSystem_product_projective_ring(polys, domain, check=check)
    else:
        R = domain.base_ring()
        CR = domain.coordinate_ring()
    if is_ProjectiveSpace(domain) or isinstance(domain, AlgebraicScheme_subscheme_projective):
        from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective_ring
        from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective_field
        from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective_finite_field
        if R in _Fields:
            if is_FiniteField(R):
                return DynamicalSystem_projective_finite_field(polys, domain=domain, check=check)
            else:
                return DynamicalSystem_projective_field(polys, domain=domain, check=check)
        return DynamicalSystem_projective_ring(polys, domain=domain, check=check)
    else:
        raise TypeError("must be projective or projective subscheme")

DynamicalSystem = DynamicalSystem_projective
##todo this should have logic for when the domain is specified to choose the right one

class DynamicalSystem_generic(SchemeMorphism_polynomial):
    """
    Base class for dynamical systems of schemes

    INPUT:

    - ``parent`` -- the parent of the morphism.

    EXAMPLES::

        sage: A.<x> = AffineSpace(QQ,1)
        sage: f = DynamicalSystem_affine([x^2+1])
        sage: type(f)
        <class 'sage.dynamics.arithmetic_dynamics.affine_ds.DynamicalSystem_affine_field'>
    """

    def __init__(self, polys, domain=None, check=True):
        """
        The Python constructor.

        INPUT:

        - ``polys`` -- SchemeMorphisms or list or tuple of polynomial or rational functions.

        - ``domain`` -- scheme or subscheme

        - ``check`` -- Boolean.

        EXAMPLES::

        sage: P.<x,y> = ProjectiveSpace(QQ,1)
        sage: f = DynamicalSystem_projective([x^2+y^2, y^2])
        sage: type(f)
        <class 'sage.dynamics.arithmetic_dynamics.projective_ds.DynamicalSystem_projective_field'>
        """
        from sage.schemes.generic.morphism import SchemeMorphism_polynomial
        from sage.schemes.affine.affine_space import is_AffineSpace, AffineSpace
        if isinstance(polys, SchemeMorphism_polynomial):
            source = polys.domain().ambient_space()
            target = polys.codomain().ambient_space()
            if not source is target:
                raise ValueError("domain and codomain must be the same")
            parent = polys.parent()
            polys = list(polys)
            if domain is None:
                domain = polys.domain()
        else:
            if domain is None:
                from sage.rings.fraction_field import is_FractionField
                R = polys[0].parent()
                if is_FractionField(R):
                    R = R.ring()
                domain = AffineSpace(R, names=R.variable_names())
        parent = End(domain)
        #from sage.schemes.affine.affine_morphism import SchemeMorphism_polynomial_affine_space
        #SchemeMorphism_polynomial_affine_space.__init__(self, parent, polys, check=check)
        SchemeMorphism_polynomial.__init__(self, parent, polys, False)

    # We copy methods of sage.categories.map.Map, to make
    # a future transition of SchemeMorphism to a sub-class of Morphism
    # easier.
    def __call__(self, x, *args, **kwds):
        """
        Do not override this method!

        For implementing application of maps, implement a method
        ``_call_(self, x)`` and/or a method ``_call_with_args(x, args, kwds)`.
        In these methods, you can assume that ``x`` belongs to the domain of
        this morphism, ``args`` is a tuple and ``kwds`` is a dict.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: f = DynamicalSystem_affine([y,x^2+y])
            sage: f([2,3])    # indirect doctest
            (3, 7)

        An example with optional arguments::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3, x*y^2])
            sage: P = PS(0,1)
            sage: f(P, check=False)     # indirect doctest
            (0 : 0)
        """
        P = parent(x)
        D = self.domain()
        if P is D: # we certainly want to call _call_/with_args
            if not args and not kwds:
                return self._call_(x)
            return self._call_with_args(x, args, kwds)
        # Is there coercion?
        converter = D._internal_coerce_map_from(P)
        if converter is None:
            try:
                return self.pushforward(x,*args,**kwds)
            except (AttributeError, TypeError, NotImplementedError):
                pass # raise TypeError, "%s must be coercible into %s"%(x, self.domain())
            # Here, we would like to do
            ##try:
            ##    x = D(x).
            ##except (TypeError, NotImplementedError):
            ##    raise TypeError, "%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(x, self.domain())
            # However, this would involve a test whether x.codomain() ==
            # self. This would trigger a Groebner basis computation, that
            # (1) could be slow and (2) could involve an even slower toy
            # implementation, resulting in a warning.
            #
            # Contract: If x is a scheme morphism point, then _call_ knows
            # what to do with it (e.g., use the _coords attribute). Otherwise,
            # we can try a conversion into the domain (e.g., if x is a list),
            # WITHOUT to trigger a Groebner basis computation.
            if kwds.get('check', True):
                if not isinstance(x, SchemeMorphism_point):
                    try:
                        x = D(x)
                    except (TypeError, NotImplementedError):
                        raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(x, self.domain()))
                elif self.domain()!=x.codomain():
                    raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(x, self.domain()))
        else:
            x = converter(x)
        if not args and not kwds:
            return self._call_(x)
        return self._call_with_args(x, args, kwds)


    def _repr_type(self):
        r"""
        Return a string representation of the type of ``self``.

        OUTPUT: string.

        EXAMPLES::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3, x*y^2])
            sage: f._repr_type()
            'Dynamical System'
        """
        return "Dynamical System"

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: String.

        EXAMPLES::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3, x*y^2])
            sage: f._repr_()
            'Dynamical System of Projective Space of dimension 1 over Rational Field\n
              Defn: Defined on coordinates by sending (x : y) to\n        (x^3 : x*y^2)'
        """
        s = "%s of %s"%(self._repr_type(), self.domain())
        d = self._repr_defn()
        if d != '':
            s += "\n  Defn: %s"%('\n        '.join(self._repr_defn().split('\n')))
        return s


    def as_scheme_morphism(self):
        """
        Return this dynamical system as :class:`SchemeMorphism_polynomial`.

        OUTPUT:

        - :class:`SchemeMorphism_polynomial`.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space'>

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2-y^2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space_field'>

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(5), 1)
            sage: f = DynamicalSystem_projective([x^2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space_finite_field'>

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: f = DynamicalSystem_affine([x^2-2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.affine.affine_morphism.SchemeMorphism_polynomial_affine_space'>

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: f = DynamicalSystem_affine([x^2-2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.affine.affine_morphism.SchemeMorphism_polynomial_affine_space_field'>

        ::

            sage: A.<x,y> = AffineSpace(GF(3), 2)
            sage: f = DynamicalSystem_affine([x^2-2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.affine.affine_morphism.SchemeMorphism_polynomial_affine_space_finite_field'>
        """
        H = End(self.domain())
        return H(list(self))

    def change_ring(self, R, check=True):
        r"""
        Returns a new :class:DynamicalSystem_projective` which is this map coerced to ``R``.

        If ``check`` is ``True``, then the initialization checks are performed.

        INPUT:

        - ``R`` -- ring or morphism.

        - ``check`` -- Boolean

        OUTPUT:

        - A new :class:`DynamicalSystem_projective` which is this map coerced to ``R``.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: f = DynamicalSystem_projective([3*x^2, y^2])
            sage: f.change_ring(GF(5))
            Dynamical System of Projective Space of dimension 1 over Finite Field of size 5
              Defn: Defined on coordinates by sending (x : y) to
                    (-2*x^2 : y^2)
        """
        f = self.as_scheme_morphism()
        F = f.change_ring(R)
        return F.as_dynamical_system()

    def specialization(self, D=None, phi=None, homset=None):
        r"""
        Specialization of this dynamical system.

        Given a family of maps defined over a polynomial ring. A specialization
        is a particular member of that family. The specialization can be specified either
        by a dictionary or a :class:`SpecializationMorphism`.

        INPUT:

        - ``D`` -- dictionary (optional)

        - ``phi`` -- SpecializationMorphism (optional)

        - ``homset`` -- homset of specialized map (optional)

        OUTPUT: :class:`DynamicalSystem`

        EXAMPLES::

            sage: R.<c> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: f = DynamicalSystem_projective([x^2 + c*y^2,y^2], domain=P)
            sage: f.specialization({c:1})
            Dynamical System of Projective Space of dimension 1 over Rational Field
                  Defn: Defined on coordinates by sending (x : y) to
                        (x^2 + y^2 : y^2)
        """
        f = self.as_scheme_morphism()
        F = f.specialization(D, phi, homset)
        return F.as_dynamical_system()
