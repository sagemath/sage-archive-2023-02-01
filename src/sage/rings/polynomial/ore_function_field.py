r"""
AUTHOR:

- Xavier Caruso (2020-04)

"""
# ***************************************************************************
#    Copyright (C) 2020 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

import sage

from sage.structure.category_object import normalize_names
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.algebras import Algebras
from sage.categories.fields import Fields
from sage.rings.ring import Algebra

from sage.rings.morphism import RingHomomorphism
from sage.categories.homset import Hom
from sage.categories.map import Section

from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing
from sage.rings.polynomial.ore_function_element import OreFunctionBaseringInjection

WORKING_CENTER_MAX_TRIES = 1000


# Generic implementation of Ore polynomial rings
#################################################

class OreFunctionField(Algebra, UniqueRepresentation):
    Element = None

    def __init__(self, ring, category=None):
        if self.Element is None:
            import sage.rings.polynomial.ore_function_element
            self.Element = sage.rings.polynomial.ore_function_element.OreFunction
        if not isinstance(ring, OrePolynomialRing):
            raise TypeError("not a Ore Polynomial Ring")
        if ring.base_ring() not in Fields():
            raise TypeError("the base ring must be a field")
        try:
            _ = ring.twisting_morphism(-1)
            self._simplification = True
        except (TypeError, NotImplementedError):
            self._simplification = False
        self._ring = ring
        base = ring.base_ring()
        category = Algebras(base).or_subcategory(category)
        Algebra.__init__(self, base, names=ring.variable_name(), normalize=True, category=category)

    def _element_constructor_(self, *args, **kwds):
        return self.Element(self, *args, **kwds)

    def _coerce_map_from_base_ring(self):
        return OreFunctionBaseringInjection(self.base_ring(), self)

    def _coerce_map_from_(self, P):
        if isinstance(P, OreFunctionField):
            return P._ring.has_coerce_map_from(self._ring)
        if isinstance(P, Parent):
            return P.has_coerce_map_from(self._ring)

    def _repr_(self):
        s = "Ore Function Field in %s over %s twisted by " % (self.variable_name(), self.base_ring())
        morphism = self.twisting_morphism()
        derivation = self.twisting_derivation()
        if derivation is None:
            s += morphism._repr_short()
        else:
            if morphism is not None:
                s += "%s and " % morphism._repr_short()
            s += derivation._repr_()
        return s

    def change_var(self, var):
        return OreFunctionField(self._ring.change_var(var)) 

    def characteristic(self):
        r"""
        Return the characteristic of the base ring of ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: R['x',sigma].characteristic()
            0

            sage: k.<u> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: k['y',Frob].characteristic()
            5
        """
        return self.base_ring().characteristic()


    def twisting_morphism(self, n=1):
        return self._ring.twisting_morphism(n)

    def twisting_derivation(self):
        return self._ring.twisting_derivation()

    def gen(self, n=0):
        return self(self._ring.gen(n))

    parameter = gen

    def gens_dict(self):
        r"""
        Return a {name: variable} dictionary of the generators of
        this Ore polynomial ring.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = SkewPolynomialRing(R,sigma)
            sage: S.gens_dict()
            {'x': x}
        """
        return dict(zip(self.variable_names(), self.gens()))

    def is_finite(self):
        return False

    def is_exact(self):
        return self.base_ring().is_exact()

    def is_sparse(self):
        return self._ring.is_sparse()

    def ngens(self):
        return 1

    def random_element(self, degree=2, monic=True, *args, **kwds):
        numerator = self._ring.random_element(degree, monic, *args, **kwds)
        denominator = self._ring.random_element(degree, True, *args, **kwds)
        return self(numerator, denominator)
    
    def is_commutative(self):
        return self._ring.is_commutative()

    def is_field(self):
        return True

    def fraction_field(self):
        return self


# Special classes for twisting morphisms with finite order
##########################################################

class SectionOreFunctionCenterInjection(Section):
    def __init__(self, embed):
        Section.__init__(self, embed)
        self._ringsection = embed._ringembed.section()
        self._simplify = not embed._codomain._simplification

    def _call_(self, x):
        numerator = x._numerator
        denominator = x._denominator
        if self._simplify:
            D = numerator.right_gcd(denominator)
            numerator, _ = numerator.right_quo_rem(D)
            denominator, _ = denominator.right_quo_rem(D)
        return self._ringsection(numerator) / self._ringsection(denominator)

    def _richcmp_(self, right, op):
        if op == op_EQ:
            return (self.domain() is right.domain()) and (self.codomain() is right.codomain())
        return NotImplemented


class OreFunctionCenterInjection(RingHomomorphism):
    def __init__(self, domain, codomain, ringembed):
        RingHomomorphism.__init__(self, Hom(domain, codomain))
        self._codomain = codomain
        self._ringembed = ringembed
        self._section = SectionOreFunctionCenterInjection(self)

    def _repr_(self):
        return "Embedding of the center of %s into this field" % self._codomain

    def _call_(self, x):
        numerator = self._ringembed(x.numerator())
        denominator = self._ringembed(x.denominator())
        return self._codomain(numerator, denominator, simplify=False)

    def _richcmp_(self, right, op):
        if op == op_EQ:
            return (self.domain() is right.domain()) and (self.codomain() is right.codomain())
        return NotImplemented

    def section(self):
        return self._section


class OreFunctionField_finite_order(OreFunctionField):
    def __init__(self, ring, category=None):
        if self.Element is None:
            self.Element = sage.rings.polynomial.ore_function_element.OreFunction_finite_order
        OreFunctionField.__init__(self, ring, category)
        self._center = {}
        self._center_variable_name = 'z'
        for i in range(WORKING_CENTER_MAX_TRIES):
            try:
                self._working_center = self.center()
                self._center_variable_name = None
                break
            except ValueError:
                self._center_variable_name = "z%s_" % i
        if self._center_variable_name is not None:
            raise NotImplementedError("unable to create the center")

    def center(self, name=None, names=None, default=False):
        if name is not None and names is not None:
            raise ValueError("you must specify the name of the variable")
        if names is None:
            if name is None:
                name = self._center_variable_name
            if name is None:
                name = 'z'
            names = (name,)
        names = normalize_names(1, names)
        name = names[0]
        if name in self._center:
            center = self._center[name]
        else:
            ring = self._ring
            ringcenter = ring.center(name, default=False)
            ringembed = ring.coerce_map_from(ringcenter)
            center = ringcenter.fraction_field()
            embed = OreFunctionCenterInjection(center, self, ringembed)
            try:
                assert not self.has_coerce_map_from(center)
                self.register_coercion(embed)
                center.register_conversion(embed.section())
            except AssertionError:
                raise ValueError("creation of coercion map fails; consider using another variable name")
            self._center[name] = center
        if default or (self._center_variable_name is None):
            self._center_variable_name = name
        return center
