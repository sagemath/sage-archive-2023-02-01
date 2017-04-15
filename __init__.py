# -*- coding: utf-8 -*-
r"""
Monkey patches to make the MacLane code work in standalone mode, i.e., without
modifying the sage source code at build time.
"""
#*****************************************************************************
#       Copyright (C) 2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Fix doctests so they work in standalone mode (when invoked with sage -t, they run within the mac_lane/ directory)
import sys, os
if hasattr(sys.modules['__main__'], 'DC') and 'standalone' in sys.modules['__main__'].DC.options.optional:
    sys.path.append(os.getcwd())
    sys.path.append(os.path.dirname(os.getcwd()))

import valuation_space
from valuation_space import DiscretePseudoValuationSpace
import trivial_valuation
from trivial_valuation import TrivialValuation, TrivialPseudoValuation
import padic_valuation
from padic_valuation import pAdicValuation
import gauss_valuation
from gauss_valuation import GaussValuation
import value_group
from value_group import DiscreteValuationCodomain, DiscreteValueGroup, DiscreteValueSemigroup
import function_field_valuation
from function_field_valuation import FunctionFieldValuation
import augmented_valuation
from augmented_valuation import AugmentedValuation
import scaled_valuation
from scaled_valuation import ScaledValuation

# fix unpickling and type checks of classes (otherwise, the instances of the
# local file and the instances that come from the mac_lane import define
# different types)
from trivial_valuation import TrivialDiscreteValuation, TrivialDiscretePseudoValuation
from function_field_valuation import FunctionFieldValuation_base, DiscreteFunctionFieldValuation_base, RationalFunctionFieldValuation_base, InducedFunctionFieldValuation_base, ClassicalFunctionFieldValuation_base, FunctionFieldFromLimitValuation, InfiniteRationalFunctionFieldValuation, FiniteRationalFunctionFieldValuation, NonClassicalRationalFunctionFieldValuation, InfiniteRationalFunctionFieldValuation, FunctionFieldMappedValuation_base, FunctionFieldExtensionMappedValuation, RationalFunctionFieldMappedValuation
from limit_valuation import LimitValuation, MacLaneLimitValuation, LimitValuation_generic
from mapped_valuation import MappedValuation_base, FiniteExtensionFromLimitValuation, FiniteExtensionFromInfiniteValuation, MappedValuation_base
from augmented_valuation import FiniteAugmentedValuation, InfiniteAugmentedValuation
from gauss_valuation import GaussValuation_generic
from valuation import DiscretePseudoValuation, DiscreteValuation, InfiniteDiscretePseudoValuation
from padic_valuation import pAdicValuation_base, pAdicValuation_int, pAdicValuation_padic, pAdicFromLimitValuation
from developing_valuation import DevelopingValuation
from augmented_valuation import AugmentedValuation_base, FinalAugmentedValuation, NonFinalAugmentedValuation, FinalFiniteAugmentedValuation, NonFinalFiniteAugmentedValuation
from inductive_valuation import FiniteInductiveValuation, FinalInductiveValuation, InfiniteInductiveValuation, NonFinalInductiveValuation
from scaled_valuation import ScaledValuation_generic

# =================
# MONKEY PATCH SAGE
# =================
import sage

# Implement Qp/Zp.valuation
sage.rings.padics.padic_generic.pAdicGeneric.valuation = lambda self: pAdicValuation(self)

# Fix contains check of rational fuction fields
def to_polynomial(self, x):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: K(x) in K._ring # indirect doctest
        True

    """
    R = x.parent()._ring
    K = x.parent().constant_base_field()
    if x.denominator() in K:
        return x.numerator()/K(x.denominator())
    raise ValueError("Only polynomials can be converted to the underlying polynomial ring")

def to_constant(self, x):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: K(1) in QQ # indirect doctest
        True

    """
    K = x.parent().constant_base_field()
    if x.denominator() in K and x.numerator() in K:
        return K(x.numerator()) / K(x.denominator())
    raise ValueError("only constants can be converted to the underlying constant field")

sage.rings.function_field.function_field.RationalFunctionField._to_polynomial = to_polynomial
sage.rings.function_field.function_field.RationalFunctionField._to_constant = to_constant
if not hasattr(sage.rings.function_field.function_field.RationalFunctionField, "__old_init__"):
    sage.rings.function_field.function_field.RationalFunctionField.__old_init__ = sage.rings.function_field.function_field.RationalFunctionField.__init__
def __init__(self, *args, **kwargs):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: K(1/2) in QQ
        True
        
    """
    self.__old_init__(*args, **kwargs)
    from sage.categories.morphism import SetMorphism
    self._ring.register_conversion(SetMorphism(self.Hom(self._ring), self._to_polynomial))
    try:
        self.constant_base_field().register_conversion(SetMorphism(self.Hom(self.constant_base_field()), self._to_constant))
    except AssertionError:
        # since #21872 there is already such a conversion
        pass

sage.rings.function_field.function_field.RationalFunctionField.__init__ = __init__
del(__init__)
del(to_polynomial)

# implement principal_part for newton polygons
r"""
TESTS::

    sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
    sage: NP = sage.geometry.newton_polygon.NewtonPolygon([(0,1),(1,0),(2,1)])
    sage: NP.principal_part()
    Infinite Newton polygon with 2 vertices: (0, 1), (1, 0) ending by an infinite line of slope 0

"""
import sage.geometry.newton_polygon
sage.geometry.newton_polygon.NewtonPolygon_element.principal_part = lambda self: sage.geometry.newton_polygon.NewtonPolygon(self.vertices(), last_slope=0)
sage.geometry.newton_polygon.NewtonPolygon_element.sides = lambda self: zip(self.vertices(), self.vertices()[1:])

# implement coercion of function fields that comes from coercion of their base fields

# Frac(K[x]) injects into K(x)
class DefaultConvertMap_unique_patched2(sage.structure.coerce_maps.DefaultConvertMap_unique):
    def is_injective(self):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: K.<x> = FunctionField(QQ)
            sage: R.fraction_field().is_subring(K) # indirect doctest
            True

        """
        from sage.categories.fields import Fields
        if self.domain() in Fields():
            return True
        raise NotImplementedError

def _coerce_map_from_(target, source):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: L.<x> = FunctionField(GaussianIntegers().fraction_field())
        sage: L.has_coerce_map_from(K)
        True

        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^3 + 1)
        sage: K.<x> = FunctionField(GaussianIntegers().fraction_field())
        sage: R.<y> = K[]
        sage: M.<y> = K.extension(y^3 + 1)
        sage: M.has_coerce_map_from(L) # not tested, base morphism is not implemented
        True

        sage: K.<x> = FunctionField(QQ)
        sage: R.<I> = K[]
        sage: L.<I> = K.extension(I^2 + 1)
        sage: M.<x> = FunctionField(GaussianIntegers().fraction_field())
        sage: M.has_coerce_map_from(L) # not tested, base_morphism is not implemented
        True

    """
    from sage.categories.function_fields import FunctionFields
    if source in FunctionFields():
        if source.base_field() is source:
            if target.base_field() is target:
                # source and target are rational function fields
                if source.variable_name() == target.variable_name():
                    # ... in the same variable
                    base_coercion = target.constant_field().coerce_map_from(source.constant_field())
                    if base_coercion is not None:
                        return source.hom([target.gen()], base_morphism=base_coercion)
        else:
            # source is an extensions of rational function fields
            base_coercion = target.coerce_map_from(source.base_field())
            if base_coercion is not None:
                # the base field of source coerces into the base field of target
                target_polynomial = source.polynomial().map_coefficients(base_coercion)
                # try to find a root of the defining polynomial in target
                if target_polynomial(target.gen()) == 0:
                    # The defining polynomial of source has a root in target,
                    # therefore there is a map. To be sure that it is
                    # canonical, we require a root of the defining polynomial
                    # of target to be a root of the defining polynomial of
                    # source (and that the variables are named equally):
                    if source.variable_name() == target.variable_name():
                        return source.hom([target.gen()], base_morphism=base_coercion)

                roots = target_polynomial.roots()
                for root, _ in roots:
                    if target_polynomial(root) == 0:
                        # The defining polynomial of source has a root in target,
                        # therefore there is a map. To be sure that it is
                        # canonical, we require the names of the roots to match
                        if source.variable_name() == repr(root):
                            return source.hom([root], base_morphism=base_coercion)
    if source is target._ring:
        return DefaultConvertMap_unique_patched2(source, target)
    if source is target._ring.fraction_field():
        return DefaultConvertMap_unique_patched2(source, target)

sage.rings.function_field.function_field.FunctionField._coerce_map_from_ = _coerce_map_from_
del(_coerce_map_from_)

# patch is_injective() for many morphisms
def patch_is_injective(method, patch_map):
    r"""
    Patch ``method`` to return ``patch_map[type]`` if it returned a result of
    ``type``.

    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: QQ.coerce_map_from(ZZ).is_injective() # indirect doctest
        True

    """
    def patched(*args, **kwargs):
        ret = method(*args, **kwargs)
        if type(ret) in patch_map:
            ret = patch_map[type(ret)](ret)
        return ret
    return patched

# a ring homomorphism from a field into a ring is injective (as it respects inverses)
class RingHomomorphism_coercion_patched(sage.rings.morphism.RingHomomorphism_coercion):
    def is_injective(self):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: QQ.coerce_map_from(ZZ).is_injective() # indirect doctest
            True

        """
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
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: QQ['x'].coerce_map_from(ZZ['x']).is_injective() # indirect doctest
            True
    
        This should be fixed in
        `sage.rings.padics.qadic_flint_CA.pAdicCoercion_CA_frac_field`
        instead::

            sage: R.<a> = ZqCA(9)
            sage: R['x'].is_subring(R.fraction_field()['x'])
            True

        """
        if self.underlying_map().codomain() is self.underlying_map().domain().fraction_field():
            # fix this in pAdicCoercion_CA_frac_field and similar
            return True
        return self.underlying_map().is_injective()         
sage.rings.polynomial.polynomial_ring.PolynomialRing_general._coerce_map_from_ = patch_is_injective(sage.rings.polynomial.polynomial_ring.PolynomialRing_general._coerce_map_from_, {sage.rings.polynomial.polynomial_ring_homomorphism.PolynomialRingHomomorphism_from_base: (lambda morphism: PolynomialRingHomomorphism_from_base_patched(morphism.parent(), morphism.underlying_map()))})

# morphisms of number fields are injective
class Q_to_quadratic_field_element_patched(sage.rings.number_field.number_field_element_quadratic.Q_to_quadratic_field_element):
    def is_injective(self):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: GaussianIntegers().fraction_field().coerce_map_from(QQ).is_injective()
            True

        """
        return True
class Z_to_quadratic_field_element_patched(sage.rings.number_field.number_field_element_quadratic.Z_to_quadratic_field_element):
    def is_injective(self):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: GaussianIntegers().fraction_field().coerce_map_from(ZZ).is_injective()
            True

        """
        return True
sage.rings.number_field.number_field.NumberField_quadratic._coerce_map_from_ = patch_is_injective(sage.rings.number_field.number_field.NumberField_quadratic._coerce_map_from_, {sage.rings.number_field.number_field_element_quadratic.Q_to_quadratic_field_element: (lambda morphism: Q_to_quadratic_field_element_patched(morphism.codomain())), sage.rings.number_field.number_field_element_quadratic.Z_to_quadratic_field_element: (lambda morphism: Z_to_quadratic_field_element_patched(morphism.codomain()))})

# the integers embed into the rationals
class Z_to_Q_patched(sage.rings.rational.Z_to_Q):
    def is_injective(self):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: QQ.coerce_map_from(ZZ).is_injective()
            True

        """
        return True
from sage.rings.all import QQ
QQ.coerce_map_from = patch_is_injective(QQ.coerce_map_from, {sage.rings.rational.Z_to_Q: (lambda morphism: Z_to_Q_patched())})

# the integers embed into their extensions in number fields
class DefaultConvertMap_unique_patched(sage.structure.coerce_maps.DefaultConvertMap_unique):
    def is_injective(self):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: CyclotomicField(5).maximal_order().coerce_map_from(ZZ).is_injective()
            True

        """
        from sage.rings.all import ZZ
        if self.domain() is ZZ or domain is int or domain is long:
            return True
        return super(DefaultConvertMap_unique, self).is_injective()

def _coerce_map_from_patched(self, domain):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: CyclotomicField(5).maximal_order().coerce_map_from(ZZ).is_injective() # indirect doctest
        True

    """
    from sage.rings.all import ZZ
    if domain is ZZ or domain is int or domain is long:
        return DefaultConvertMap_unique_patched(domain, self)
    return False

sage.rings.number_field.order.Order._coerce_map_from_ = _coerce_map_from_patched
del(_coerce_map_from_patched)

# quotient rings embed if their underlying rings do
class DefaultConvertMap_unique_patched3(sage.structure.coerce_maps.DefaultConvertMap_unique):
    def is_injective(self):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = ZZ[]
            sage: S.<x> = QQ[]
            sage: S.quo(x^2 + 1).coerce_map_from(R.quo(x^2 + 1)).is_injective()
            True

        """
        if self.codomain().base().coerce_map_from(self.domain().base()).is_injective():
            return True
        raise NotImplementedError

sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRing_generic._coerce_map_from_original = sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRing_generic._coerce_map_from_
def _coerce_map_from_patched(self, domain):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: R.<x> = ZZ[]
        sage: S.<x> = QQ[]
        sage: S.quo(x^2 + 1).coerce_map_from(R.quo(x^2 + 1)).is_injective() # indirect doctest
        True

    """
    from sage.rings.polynomial.polynomial_quotient_ring import is_PolynomialQuotientRing
    if is_PolynomialQuotientRing(domain) and domain.modulus() == self.modulus():
        if self.base().has_coerce_map_from(domain.base()):
            return DefaultConvertMap_unique_patched3(domain, self)
    from sage.rings.fraction_field import is_FractionField
    if is_FractionField(domain):
        # this should be implemented on a much higher level:
        # if there is a morphism R -> K then there is a morphism Frac(R) -> K
        if self.has_coerce_map_from(domain.base()):
            return True
    return self._coerce_map_from_original(domain)
sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRing_generic._coerce_map_from_ = _coerce_map_from_patched
del(_coerce_map_from_patched)

# a ring embeds into its field of fractions
class CallableConvertMap_patched(sage.rings.fraction_field.CallableConvertMap):
    def is_injective(self):
        r"""
        TESTS::
        
            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: R.is_subring(R.fraction_field()) # indirect doctest
            True

        """
        return True

    def is_surjective(self):
        r"""
        TESTS::
        
            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: R.fraction_field().coerce_map_from(R).is_surjective()
            False

        """
        return False

sage.rings.fraction_field.CallableConvertMap = CallableConvertMap_patched

# inverses of quotient ring elements
def inverse_of_unit(self):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: R.<x> = ZZ[]
        sage: S = R.quo(x^2+x+1)
        sage: S(1).inverse_of_unit()
        1

    """
    inverse = ~self
    if inverse.parent() is self.parent():
        return inverse
    raise NotImplementedError

sage.rings.polynomial.polynomial_quotient_ring_element.PolynomialQuotientRingElement.inverse_of_unit = inverse_of_unit
del(inverse_of_unit)

# factorization in polynomial quotient fields
def _factor_univariate_polynomial(self, f):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: K = GF(2)
        sage: R.<x> = K[]
        sage: L.<x> = K.extension(x^2 + x + 1)
        sage: R.<y> = L[]
        sage: L.<y> = L.extension(y^2 + y + x)
        sage: R.<T> = L[]
        sage: (T^2 + T + x).factor() # indirect doctest
        (T + y) * (T + y + 1)

    """
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

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: k.<a> = GF(4)
        sage: R.<b> = k[]
        sage: l.<b> = k.extension(b^2+b+a); l
        Univariate Quotient Polynomial Ring in b over Finite Field in a of size 2^2 with modulus b^2 + b + a
        sage: from_ll,to_ll, ll = l.absolute_extension(); ll
        Finite Field in v4 of size 2^4
        sage: all([to_ll(from_ll(ll.gen()**i)) == ll.gen()**i for i in range(ll.degree())])
        True

        sage: R.<c> = l[]
        sage: m.<c> = l.extension(c^2+b*c+b); m
        Univariate Quotient Polynomial Ring in c over Univariate Quotient Polynomial Ring in b over Finite Field in a of size 2^2 with modulus b^2 + b + a with modulus c^2 + b*c + b
        sage: from_mm, to_mm, mm = m.absolute_extension(); mm
        Finite Field in v8 of size 2^8
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
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: K = GF(2)
        sage: R.<x> = K[]
        sage: L.<x> = K.extension(x^2 + x + 1)
        sage: R.<y> = L[]
        sage: L.<y> = L.extension(y^2 + y + x)
        sage: L.vector_space()
        Vector space of dimension 2 over Finite Field in x of size 2^2

    """
    if not self.base().base_ring().is_field():
        raise ValueError

    return self.base().base_ring()**self.modulus().degree()
sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRing_generic.vector_space = vector_space
del(vector_space)

# make some_elements() non-trivial for number fields
def some_elements(self):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: K = GaussianIntegers().fraction_field()
        sage: list(K.some_elements())
        [I, 0, 1, 1/2, 2*I, -I, -2, 0, 0]

    """
    for element in self.polynomial_ring().some_elements():
        yield element(self.gen())
sage.rings.number_field.number_field.NumberField_generic.some_elements = some_elements
del(some_elements)

# make some_elements() deterministic for function fields
def some_elements(self):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: list(K.some_elements()) == list(K.some_elements())
        True

    """
    for num in self._ring.some_elements():
        for den in self._ring.some_elements():
            if den != 0:
                yield self(num) / self(den)
sage.rings.function_field.function_field.RationalFunctionField.some_elements = some_elements
del(some_elements)

def some_elements(self):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: list(L.some_elements()) == list(L.some_elements())
        True

    """
    for element in self._ring.some_elements():
        yield self(element)
sage.rings.function_field.function_field.FunctionField_polymod.some_elements = some_elements
del(some_elements)

# make some_elements() non-trivial for fraction fields
def some_elements(self):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: R.<x> = QQ[]
        sage: K = R.fraction_field()
        sage: len(list(K.some_elements()))
        72

    """
    for num in self.ring().some_elements():
        for den in self.ring().some_elements():
            if den != 0:
                yield self(num) / self(den)
sage.rings.fraction_field.FractionField_generic.some_elements = some_elements

# make some_elements() non-trivial for orders in number fields
def some_elements(self):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: R = GaussianIntegers()
        sage: list(R.some_elements())
        [I, 0, 1, 2*I, -I, -2, 0, 0]

    """
    for element in self.fraction_field().some_elements():
        if element in self:
            yield self(element)
sage.rings.number_field.order.Order.some_elements = some_elements
del(some_elements)

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
register_factory_unpickle("ScaledValuation", ScaledValuation)
