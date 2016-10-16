from .valuation_space import DiscretePseudoValuationSpace, DiscreteValuationSpace
from .trivial_valuation import TrivialDiscreteValuation, TrivialDiscretePseudoValuation
from .padic_valuation import pAdicValuation
from .gauss_valuation import GaussValuation
from .value_group import DiscreteValuationCodomain, DiscreteValueGroup
import valuation, padic_valuation, mac_lane

# =================
# MONKEY PATCH SAGE
# =================
import sage

# Implement Qp/Zp.valuation
sage.rings.padics.padic_generic.pAdicGeneric.valuation = lambda self: pAdicValuation(self)

# fix unpickling of factories
from sage.structure.factory import register_factory_unpickle
register_factory_unpickle("pAdicValuation", pAdicValuation)
register_factory_unpickle("GaussValuation", GaussValuation)
