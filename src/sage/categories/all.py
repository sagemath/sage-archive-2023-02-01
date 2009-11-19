from category import    is_Category, Category, HomCategory, AbstractCategory

from category_types import(
                        Elements,
                        Sequences,
                        SimplicialComplexes,
                        ChainComplexes,
)

from tensor     import CategoryWithTensorProduct, TensorCategory, tensor
from cartesian_product import CategoryWithCartesianProduct, CartesianProductCategory, cartesian_product
from dual       import DualityCategory

from functor  import (is_Functor,
                      ForgetfulFunctor,
                      IdentityFunctor)

from homset   import (Hom, hom, is_Homset,
                      End, end, is_Endset,
                      Homset, HomsetWithBase)

from morphism import Morphism, is_Morphism

from basic import *

from g_sets import GSets
from pointed_sets import PointedSets

from sets_with_partial_maps import SetsWithPartialMaps
from groupoid import Groupoid

# enumerated sets
from enumerated_sets import EnumeratedSets
from finite_enumerated_sets import FiniteEnumeratedSets
from infinite_enumerated_sets import InfiniteEnumeratedSets

# fields
from quotient_fields import QuotientFields
from finite_fields import FiniteFields
from number_fields import NumberFields

# modules
from left_modules import LeftModules
from right_modules import RightModules
from bimodules import Bimodules

from modules import Modules
RingModules = Modules
from vector_spaces import VectorSpaces

# (hopf) algebra structures
from algebras import Algebras
from commutative_algebras import CommutativeAlgebras
from coalgebras import Coalgebras
from bialgebras import Bialgebras
from hopf_algebras import HopfAlgebras

# specific algebras
from monoid_algebras import MonoidAlgebras
from group_algebras import GroupAlgebras
from matrix_algebras import MatrixAlgebras

# ideals
from ring_ideals import RingIdeals
Ideals = RingIdeals
from commutative_ring_ideals import CommutativeRingIdeals
from algebra_modules import AlgebraModules
from algebra_ideals import AlgebraIdeals
from commutative_algebra_ideals import CommutativeAlgebraIdeals

# schemes and varieties
from modular_abelian_varieties import ModularAbelianVarieties
from schemes import Schemes

# * with basis
from modules_with_basis import ModulesWithBasis
FreeModules = ModulesWithBasis
from hecke_modules            import HeckeModules
from algebras_with_basis      import AlgebrasWithBasis
from coalgebras_with_basis    import CoalgebrasWithBasis
from bialgebras_with_basis    import BialgebrasWithBasis
from hopf_algebras_with_basis import HopfAlgebrasWithBasis

# finite dimensional * with basis
from finite_dimensional_modules_with_basis       import FiniteDimensionalModulesWithBasis
from finite_dimensional_algebras_with_basis      import FiniteDimensionalAlgebrasWithBasis
from finite_dimensional_coalgebras_with_basis    import FiniteDimensionalCoalgebrasWithBasis
from finite_dimensional_bialgebras_with_basis    import FiniteDimensionalBialgebrasWithBasis
from finite_dimensional_hopf_algebras_with_basis import FiniteDimensionalHopfAlgebrasWithBasis

# graded *
from graded_modules       import GradedModules
from graded_algebras      import GradedAlgebras
from graded_coalgebras    import GradedCoalgebras
from graded_bialgebras    import GradedBialgebras
from graded_hopf_algebras import GradedHopfAlgebras

# graded * with basis
from graded_modules_with_basis       import GradedModulesWithBasis
from graded_algebras_with_basis      import GradedAlgebrasWithBasis
from graded_coalgebras_with_basis    import GradedCoalgebrasWithBasis
from graded_bialgebras_with_basis    import GradedBialgebrasWithBasis
from graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis




