import catalog as crystals

def CrystalOfLetters(cartan_type, element_print_style=None, dual=None):
    """
    TESTS::

        sage: C = CrystalOfLetters(['A',5])
        doctest:...: DeprecationWarning: this being removed from the global namespace. Use crystals.Letters instead
        See http://trac.sagemath.org/15882 for details.
    """
    from sage.misc.superseded import deprecation
    from letters import CrystalOfLetters
    deprecation(15882,'this being removed from the global namespace. Use crystals.Letters instead')
    return CrystalOfLetters(cartan_type, element_print_style, dual)

def FastCrystal(cartan_type, shape, format="string"):
    """
    TESTS::

        sage: C = FastCrystal(['A',2],shape=[4,1])
        doctest:...: DeprecationWarning: this being removed from the global namespace. Use crystals.FastRankTwo instead
        See http://trac.sagemath.org/15882 for details.
    """
    from sage.misc.superseded import deprecation
    from fast_crystals import FastCrystal
    deprecation(15882,'this being removed from the global namespace. Use crystals.FastRankTwo instead')
    return FastCrystal(cartan_type, shape, format)

def CrystalOfTableaux(cartan_type, shapes=None, shape=None):
    """
    TESTS::

        sage: C = CrystalOfTableaux(['A',2],shape=[4,1])
        doctest:...: DeprecationWarning: this being removed from the global namespace. Use crystals.Tableaux instead
        See http://trac.sagemath.org/15882 for details.
    """
    from sage.misc.superseded import deprecation
    from tensor_product import CrystalOfTableaux
    deprecation(15882,'this being removed from the global namespace. Use crystals.Tableaux instead')
    return CrystalOfTableaux(cartan_type, shapes, shape)

def HighestWeightCrystal(dominant_weight, model=None):
    """
    TESTS::

        sage: C = CartanType(['E',6])
        sage: La = C.root_system().weight_lattice().fundamental_weights()
        sage: T = HighestWeightCrystal(La[1])
        doctest:...: DeprecationWarning: this being removed from the global namespace. Use crystals.HighestWeight instead
        See http://trac.sagemath.org/15882 for details.
    """
    from sage.misc.superseded import deprecation
    from highest_weight_crystals import HighestWeightCrystal
    deprecation(15882,'this being removed from the global namespace. Use crystals.HighestWeight instead')
    return HighestWeightCrystal(dominant_weight, model)

def DirectSumOfCrystals(crystals, **options):
    """
    TESTS::

        sage: C = crystals.Letters(['A',2])
        sage: C2 = DirectSumOfCrystals([C,C])
        doctest:...: DeprecationWarning: this being removed from the global namespace. Use crystals.DirectSum instead
        See http://trac.sagemath.org/15882 for details.
    """
    from sage.misc.superseded import deprecation
    from direct_sum import DirectSumOfCrystals
    deprecation(15882,'this being removed from the global namespace. Use crystals.DirectSum instead')
    return DirectSumOfCrystals(crystals, **options)

def TensorProductOfCrystals(*crystals, **options):
    """
    TESTS::

        sage: C = crystals.Letters(['A',2])
        sage: C2 = TensorProductOfCrystals(C,C)
        doctest:...: DeprecationWarning: this being removed from the global namespace. Use crystals.TensorProduct instead
        See http://trac.sagemath.org/15882 for details.
    """
    from sage.misc.superseded import deprecation
    from tensor_product import TensorProductOfCrystals
    deprecation(15882,'this being removed from the global namespace. Use crystals.TensorProduct instead')
    return TensorProductOfCrystals(*crystals, **options)

from spins import CrystalOfSpins
from spins import CrystalOfSpinsPlus
from spins import CrystalOfSpinsMinus
from kyoto_path_model import KyotoPathModel
from affine import AffineCrystalFromClassical
from affine import AffineCrystalFromClassicalAndPromotion
from kirillov_reshetikhin import KirillovReshetikhinCrystal
from alcove_path import CrystalOfAlcovePaths
from alcove_path import ClassicalCrystalOfAlcovePaths
from littelmann_path import CrystalOfLSPaths, CrystalOfProjectedLevelZeroLSPaths
from generalized_young_walls import InfinityCrystalOfGeneralizedYoungWalls, CrystalOfGeneralizedYoungWalls
from infinity_crystals import InfinityCrystalOfTableaux
from elementary_crystals import TCrystal, RCrystal, ElementaryCrystal, ComponentCrystal
from monomial_crystals import InfinityCrystalOfNakajimaMonomials, CrystalOfNakajimaMonomials

