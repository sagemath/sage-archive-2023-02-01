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

def CrystalOfSpins(ct):
    """
    TESTS::

        sage: C = CrystalOfSpins(['B',3])
        doctest:...: DeprecationWarning: this being removed from the global namespace. Use crystals.Spins instead
        See http://trac.sagemath.org/15882 for details.
    """
    from sage.misc.superseded import deprecation
    from spins import CrystalOfSpins
    deprecation(15882,'this being removed from the global namespace. Use crystals.Spins instead')
    return CrystalOfSpins(ct)

def CrystalOfSpinsPlus(ct):
    """
    TESTS::

        sage: C = CrystalOfSpinsPlus(['D',4])
        doctest:...: DeprecationWarning: this being removed from the global namespace. Use crystals.SpinsPlus instead
        See http://trac.sagemath.org/15882 for details.
    """
    from sage.misc.superseded import deprecation
    from spins import CrystalOfSpinsPlus
    deprecation(15882,'this being removed from the global namespace. Use crystals.SpinsPlus instead')
    return CrystalOfSpinsPlus(ct)

def CrystalOfSpinsMinus(ct):
    """
    TESTS::

        sage: C = CrystalOfSpinsMinus(['D',4])
        doctest:...: DeprecationWarning: this being removed from the global namespace. Use crystals.SpinsMinus instead
        See http://trac.sagemath.org/15882 for details.
    """
    from sage.misc.superseded import deprecation
    from spins import CrystalOfSpinsMinus
    deprecation(15882,'this being removed from the global namespace. Use crystals.SpinsMinus instead')
    return CrystalOfSpinsMinus(ct)

def KyotoPathModel(crystals, weight):
    """
    TESTS::

        sage: B = crystals.kirillov_reshetikhin.KashiwaraNakashimaTableaux(['A',2,1], 1,1)
        sage: L = RootSystem(['A',2,1]).weight_space()
        sage: C = KyotoPathModel(B, L.fundamental_weight(0))
        doctest:...: DeprecationWarning: this being removed from the global namespace. Use crystals.KyotoPathModel instead
        See http://trac.sagemath.org/15882 for details.
    """
    from sage.misc.superseded import deprecation
    from kyoto_path_model import KyotoPathModel
    deprecation(15882,'this being removed from the global namespace. Use crystals.KyotoPathModel instead')
    return KyotoPathModel(crystals, weight)

def AffineCrystalFromClassical(cartan_type, *args, **options):
    """
    TESTS::

        sage: B = crystals.Tableaux(['A',2],shape=[1])
        sage: C = AffineCrystalFromClassical(['A',2,1], B)
        doctest:...: DeprecationWarning: this being removed from the global namespace. Use crystals.AffineFromClassical instead
        See http://trac.sagemath.org/15882 for details.
    """
    from sage.misc.superseded import deprecation
    from affine import AffineCrystalFromClassical
    deprecation(15882,'this being removed from the global namespace. Use crystals.AffineFromClassical instead')
    return AffineCrystalFromClassical(cartan_type, *args, **options)

def AffineCrystalFromClassicalAndPromotion(cartan_type, *args, **options):
    """
    TESTS::

        sage: n = 2
        sage: B = crystals.Tableaux(['A',n],shape=[1])
        sage: pr = attrcall("promotion")
        sage: pr_inverse = attrcall("promotion_inverse")
        sage: C = AffineCrystalFromClassicalAndPromotion(['A',n,1],B,pr,pr_inverse,1)
        doctest:...: DeprecationWarning: this being removed from the global namespace. Use crystals.AffineFromClassicalAndPromotion instead
        See http://trac.sagemath.org/15882 for details.
    """
    from sage.misc.superseded import deprecation
    from affine import AffineCrystalFromClassicalAndPromotion
    deprecation(15882,'this being removed from the global namespace. Use crystals.AffineFromClassicalAndPromotion instead')
    return AffineCrystalFromClassicalAndPromotion(cartan_type, *args, **options)

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


from kirillov_reshetikhin import KirillovReshetikhinCrystal
from alcove_path import CrystalOfAlcovePaths
from alcove_path import ClassicalCrystalOfAlcovePaths
from littelmann_path import CrystalOfLSPaths, CrystalOfProjectedLevelZeroLSPaths
from generalized_young_walls import InfinityCrystalOfGeneralizedYoungWalls, CrystalOfGeneralizedYoungWalls
from infinity_crystals import InfinityCrystalOfTableaux
from elementary_crystals import TCrystal, RCrystal, ElementaryCrystal, ComponentCrystal
from monomial_crystals import InfinityCrystalOfNakajimaMonomials, CrystalOfNakajimaMonomials

