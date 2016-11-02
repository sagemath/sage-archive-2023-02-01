# -*- coding: utf-8 -*-

#*****************************************************************************
#       Copyright (C) 2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .valuation_space import DiscretePseudoValuationSpace
from .trivial_valuation import TrivialValuation, TrivialPseudoValuation
from .padic_valuation import pAdicValuation
from .gauss_valuation import GaussValuation
from .value_group import DiscreteValuationCodomain, DiscreteValueGroup
from .function_field_valuation import FunctionFieldValuation
from .augmented_valuation import AugmentedValuation

# fix unpickling and type checks of classes (otherwise, the instances of the
# local file and the instances that come from the mac_lane import define
# different types)
from .trivial_valuation import TrivialDiscreteValuation, TrivialDiscretePseudoValuation
from .function_field_valuation import FunctionFieldValuation_base, RationalFunctionFieldValuation_base, InducedFunctionFieldValuation_base, ClassicalFunctionFieldValuation_base, FunctionFieldFromLimitValuation
from .limit_valuation import LimitValuation, MacLaneLimitValuation, FiniteExtensionFromInfiniteValuation, FiniteExtensionFromLimitValuation, LimitValuation_generic
from .augmented_valuation import FiniteAugmentedValuation, InfiniteAugmentedValuation
from .gauss_valuation import GaussValuation_generic
from .valuation import DiscretePseudoValuation, DiscreteValuation, InfiniteDiscretePseudoValuation
from .padic_valuation import pAdicValuation_base, pAdicValuation_int, pAdicValuation_padic, pAdicFromLimitValuation

# =================
# MONKEY PATCH SAGE
# =================
import sage

# Implement Qp/Zp.valuation
sage.rings.padics.padic_generic.pAdicGeneric.valuation = lambda self: pAdicValuation(self)

# Fix contains check of rational fuction fields
r"""
sage: K.<x> = FunctionField(QQ)
sage: K(1) in QQ
True
sage: K(x) in K._ring
True
"""
def to_polynomial(self, x):
    R = x.parent()._ring
    K = x.parent().constant_base_field()
    if x.denominator() in K:
        return x.numerator()/K(x.denominator())
    raise ValueError("Only polynomials can be converted to the underlying polynomial ring")

sage.rings.function_field.function_field.RationalFunctionField._to_polynomial = to_polynomial
sage.rings.function_field.function_field.RationalFunctionField.__old_init__ = sage.rings.function_field.function_field.RationalFunctionField.__init__
def __init__(self, *args, **kwargs):
    self.__old_init__(*args, **kwargs)
    from sage.categories.morphism import SetMorphism
    self._ring.register_conversion(SetMorphism(self.Hom(self._ring), self._to_polynomial))

sage.rings.function_field.function_field.RationalFunctionField.__init__ = __init__
del(__init__)
del(to_polynomial)

# implement principal_part for newton polygons
import sage.geometry.newton_polygon
sage.geometry.newton_polygon.NewtonPolygon_element.principal_part = lambda self: sage.geometry.newton_polygon.NewtonPolygon(self.vertices(), last_slope=0)
sage.geometry.newton_polygon.NewtonPolygon_element.sides = lambda self: zip(self.vertices(), self.vertices()[1:])

# implement coercion of function fields that comes from coercion of their base fields
def _coerce_map_from_(target, source):
    from sage.categories.function_fields import FunctionFields
    if source in FunctionFields():
        if source.base_field() is source and target.base_field() is target:
            if source.variable_name() == target.variable_name():
                # source and target are rational function fields in the same variable
                base_coercion = target.constant_field().coerce_map_from(source.constant_field())
                if base_coercion is not None:
                    return source.hom([target.gen()], base_morphism=base_coercion)
        if source.base_field() is not source and target.base_field() is not target:
            # source and target are extensions of rational function fields
            base_coercion = target.base_field().coerce_map_from(source.base_field())
            if base_coercion is not None:
                # The base field of source coerces into the base field of target.
                if source.polynomial().map_coefficients(base_coercion)(target.gen()) == 0:
                    # The defining polynomial of source has a root in target,
                    # therefore there is a map. To be sure that it is
                    # canonical, we require a root of the defining polynomial
                    # of target to be a root of the defining polynomial of
                    # source (and that the variables are named equally):
                    if source.variable_name() == target.variable_name():
                        return source.hom([target.gen()], base_morphism=base_coercion)

sage.rings.function_field.function_field.FunctionField._coerce_map_from_ = _coerce_map_from_
del(_coerce_map_from_)

# patch is_injective() for many morphisms
def patch_is_injective(method, patch_map):
    def patched(*args, **kwargs):
        ret = method(*args, **kwargs)
        if type(ret) in patch_map:
            ret = patch_map[type(ret)](ret)
        return ret
    return patched

# a ring homomorphism from a field into a ring is injective (as it respects inverses)
class RingHomomorphism_coercion_patched(sage.rings.morphism.RingHomomorphism_coercion):
    def is_injective(self):
        from sage.categories.fields import Fields
        if self.domain() in Fields(): return True
        coercion = self.codomain().coerce_map_from(self.domain())
        if coercion is not None:
            return coercion.is_injective()
        raise NotImplementedError
sage.rings.homset.RingHomset_generic.natural_map = patch_is_injective(sage.rings.homset.RingHomset_generic.natural_map, {sage.rings.morphism.RingHomomorphism_coercion: (lambda coercion: RingHomomorphism_coercion_patched(coercion.parent()))})

# a morphism of polynomial rings which is induced by a ring morphism on the base is injective if the morphis on the base is
class PolynomialRingHomomorphism_from_base_patched(sage.rings.polynomial.polynomial_ring_homomorphism.PolynomialRingHomomorphism_from_base):
    def is_injective(self):
        return self.underlying_map().is_injective()         
sage.rings.polynomial.polynomial_ring.PolynomialRing_general._coerce_map_from_ = patch_is_injective(sage.rings.polynomial.polynomial_ring.PolynomialRing_general._coerce_map_from_, {sage.rings.polynomial.polynomial_ring_homomorphism.PolynomialRingHomomorphism_from_base: (lambda morphism: PolynomialRingHomomorphism_from_base_patched(morphism.parent(), morphism.underlying_map()))})

# morphisms of number fields are injective
class Q_to_quadratic_field_element_patched(sage.rings.number_field.number_field_element_quadratic.Q_to_quadratic_field_element):
    def is_injective(self): return True
class Z_to_quadratic_field_element_patched(sage.rings.number_field.number_field_element_quadratic.Z_to_quadratic_field_element):
    def is_injective(self): return True
sage.rings.number_field.number_field.NumberField_quadratic._coerce_map_from_ = patch_is_injective(sage.rings.number_field.number_field.NumberField_quadratic._coerce_map_from_, {sage.rings.number_field.number_field_element_quadratic.Q_to_quadratic_field_element: (lambda morphism: Q_to_quadratic_field_element_patched(morphism.codomain())), sage.rings.number_field.number_field_element_quadratic.Z_to_quadratic_field_element: (lambda morphism: Z_to_quadratic_field_element_patched(morphism.codomain()))})

# the integers embed into the rationals
class Z_to_Q_patched(sage.rings.rational.Z_to_Q):
    def is_injective(self): return True
from sage.rings.all import QQ
QQ.coerce_map_from = patch_is_injective(QQ.coerce_map_from, {sage.rings.rational.Z_to_Q: (lambda morphism: Z_to_Q_patched())})

# the integers embed into their extensions in number fields
class DefaultConvertMap_unique_patched(sage.structure.coerce_maps.DefaultConvertMap_unique):
    def is_injective(self): return True
def _coerce_map_from_patched(self, domain):
    from sage.rings.all import ZZ
    if domain is ZZ or domain is int or domain is long:
        return DefaultConvertMap_unique_patched(domain, self)
    return False
from sage.rings.number_field.order import Order
Order._coerce_map_from_ = _coerce_map_from_patched
del(_coerce_map_from_patched)

# register modules at some standard places so imports work as exepcted
r"""
sage: from sage.rings.valuation.gauss_valuation import GaussValuation
"""
import imp, sys
sage.rings.valuation = sys.modules['sage.rings.valuation'] = imp.new_module('sage.rings.valuation')
sage.rings.valuation.gauss_valuation = sys.modules['sage.rings.valuation.gauss_valuation'] = gauss_valuation
sage.rings.valuation.valuation = sys.modules['sage.rings.valuation.valuation'] = valuation
sage.rings.valuation.valuation_space = sys.modules['sage.rings.valuation.valuation_space'] = valuation_space
sage.rings.valuation.augmented_valuation = sys.modules['sage.rings.valuation.augmented_valuation'] = augmented_valuation
sage.rings.function_field.function_field_valuation = sys.modules['sage.rings.function_field.function_field_valuation'] = function_field_valuation

# fix unpickling of factories
from sage.structure.factory import register_factory_unpickle
register_factory_unpickle("pAdicValuation", pAdicValuation)
register_factory_unpickle("GaussValuation", GaussValuation)
register_factory_unpickle("TrivialValuation", TrivialValuation)
register_factory_unpickle("TrivialPseudoValuation", TrivialPseudoValuation)
register_factory_unpickle("FunctionFieldValuation", FunctionFieldValuation)
register_factory_unpickle("AugmentedValuation", AugmentedValuation)
register_factory_unpickle("LimitValuation", LimitValuation)
