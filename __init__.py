from .valuation_space import DiscretePseudoValuationSpace, DiscreteValuationSpace
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
from .function_field_valuation import FunctionFieldValuation_base, RationalFunctionFieldValuation_base, InducedFunctionFieldValuation_base

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

import imp, sys
# register modules at some standard places so imports work as exepcted
r"""
sage: from sage.rings.valuation.gauss_valuation import GaussValuation
"""
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
