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
from .developing_valuation import DevelopingValuation

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

# factorization in polynomial quotient fields
def _factor_univariate_polynomial(self, f):
    from sage.structure.factorization import Factorization

    if f.is_zero():
        raise ValueError("factorization of 0 is not defined")
    elif f.degree() <= 1:
        return Factorization([(f,1)])

    from_absolute_field, to_absolute_field, absolute_field = self.absolute_extension()

    F = f.map_coefficients(lambda c:to_absolute_field(c), absolute_field).factor()
    return Factorization([(g.map_coefficients(lambda c:from_absolute_field(c), self), e) for g,e in F], unit=from_absolute_field(F.unit()))
sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRing_generic._factor_univariate_polynomial = _factor_univariate_polynomial
del(_factor_univariate_polynomial)

# factorization needs to go to the absolute field and back
from sage.misc.cachefunc import cached_method
@cached_method
def absolute_extension(self):
    """
    Return a ring isomorphic to this ring which is not a
    :class:`PolynomialQuotientRing` but of a type which offers more
    functionality.

    INUPT:

    - ``name`` -- a list of strings or ``None`` (default: ``None``), the
      name of the generator of the absolute extension. If ``None``, this
      will be the same as the name of the generator of this ring.

    EXAMPLES::

        sage: k.<a> = GF(4)
        sage: R.<b> = k[]
        sage: l.<b> = k.extension(b^2+b+a); l
        Univariate Quotient Polynomial Ring in b over Finite Field in a of size 2^2 with modulus b^2 + b + a
        sage: from_ll,to_ll, ll = l.absolute_extension(); ll
        Finite Field in z4 of size 2^4
        sage: all([to_ll(from_ll(ll.gen()**i)) == ll.gen()**i for i in range(ll.degree())])
        True

        sage: R.<c> = l[]
        sage: m.<c> = l.extension(c^2+b*c+b); m
        Univariate Quotient Polynomial Ring in c over Univariate Quotient Polynomial Ring in b over Finite Field in a of size 2^2 with modulus b^2 + b + a with modulus c^2 + b*c + b
        sage: from_mm, to_mm, mm = m.absolute_extension(); mm
        Finite Field in z8 of size 2^8
        sage: all([to_mm(from_mm(mm.gen()**i)) == mm.gen()**i for i in range(mm.degree())])
        True

    """
    from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
    if not self.is_field():
        raise NotImplementedError("absolute_extension() only implemented for fields")

    if self.is_finite():
        if self.base_ring().is_prime_field():
            if self.modulus().degree() == 1:
                ret = self.base_ring()
                from sage.categories.homset import Hom
                from sage.categories.morphism import SetMorphism
                to_ret = SetMorphism(Hom(self, ret), lambda x: x.lift()[0])
                from_ret = self.coerce_map_from(ret)
                return from_ret, to_ret, ret
            else:
                raise NotImplementedError

        if isinstance(self.base_ring(), PolynomialQuotientRing_generic):
            abs_base_to_base, base_to_abs_base, abs_base = self.base_ring().absolute_extension()
            modulus_over_abs_base = self.modulus().map_coefficients(lambda c:base_to_abs_base(c), abs_base)
            new_self = modulus_over_abs_base.parent().quo(modulus_over_abs_base)
            ret_to_new_self, new_self_to_ret, ret = new_self.absolute_extension()
            from_ret = ret.hom([ret_to_new_self(ret.gen()).lift().map_coefficients(abs_base_to_base, self.base_ring())(self.gen())], check=False)
            to_ret = lambda x: x.lift().map_coefficients(lambda c: new_self_to_ret(base_to_abs_base(c)), ret)(new_self_to_ret(new_self.gen()))
            from sage.categories.homset import Hom
            from sage.categories.morphism import SetMorphism
            to_ret = SetMorphism(Hom(self, ret), to_ret)
            return from_ret, to_ret, ret
        else:
            N = self.cardinality()
            from sage.rings.all import GF
            ret = GF(N,prefix='v')
            base_to_ret = self.base_ring().hom([self.base_ring().modulus().change_ring(ret).roots()[0][0]])
            im_gen = self.modulus().map_coefficients(lambda c:base_to_ret(c), ret).roots()[0][0]
            to_ret = lambda x: x.lift().map_coefficients(base_to_ret, ret)(im_gen)
            from sage.categories.homset import Hom
            from sage.categories.morphism import SetMorphism
            to_ret = SetMorphism(Hom(self, ret), to_ret)

            basis = [self.gen()**i*self.base_ring().gen()**j for i in range(self.degree()) for j in range(self.base_ring().degree())]
            assert len(basis) == ret.degree()
            basis_in_ret = [to_ret(b)._vector_() for b in basis]
            from sage.matrix.constructor import matrix
            A = matrix(basis_in_ret)
            assert A.is_square()
            x = A.solve_left(A.column_space().basis()[1])
            from_ret = ret.hom([sum(c*b for c,b in zip(x.list(),basis))], check=False)
            return from_ret, to_ret, ret
    else:
        raise NotImplementedError
sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRing_generic.absolute_extension = absolute_extension
del(absolute_extension)

# factorization needs some linear algebra (it seems)
def vector_space(self):
    if not self.base().base_ring().is_field():
        raise ValueError

    return self.base().base_ring()**self.modulus().degree()
sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRing_generic.vector_space = vector_space
del(vector_space)

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
