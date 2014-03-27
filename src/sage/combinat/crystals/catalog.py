"""
Catalog of Crystals

This is a catalog of crystals that are currently in Sage:

* :class:`Letters <sage.combinat.crystals.letters.CrystalOfLetters>`
* :class:`Tableaux <sage.combinat.crystals.tensor_product.CrystalOfTableaux>`
* :func:`HighestWeight <sage.combinat.crystals.highest_weight_crystals.HighestWeightCrystal>`
* :class:`NakajimaMonomials <sage.combinat.crystals.monomial_crystals.CrystalOfNakajimaMonomials>`
* :class:`RiggedConfigurations <sage.combinat.rigged_configurations.rc_crystal.CrystalOfRiggedConfigurations>`
* :class:`GeneralizedYoungWalls <sage.combinat.crystals.generalized_young_walls.CrystalOfGeneralizedYoungWalls>`
* :class:`KyotoPathModel <sage.combinat.crystals.kyoto_path_model.KyotoPathModel>`
* :class:`Spins <sage.combinat.crystals.spins.CrystalOfSpins>`
* :class:`SpinsPlus <sage.combinat.crystals.spins.CrystalOfSpinsPlus>`
* :class:`SpinsMinus <sage.combinat.crystals.spins.CrystalOfSpinsMinus>`
* :class:`FastRankTwo <sage.combinat.crystals.fast_crystals.FastCrystal>`
* :class:`~sage.combinat.crystals.affine.AffineCrystalFromClassical`
* :class:`~sage.combinat.crystals.affine.AffineCrystalFromClassicalAndPromotion`

Functorial constructions:

* :class:`tensor products <sage.combinat.crystals.tensor_product.TensorProductOfCrystals>`
* :class:`direct sums <sage.combinat.crystals.direct_sum.DirectSumOfCrystals>`

Subcatalogs:

* :mod:`Kirillov-Reshetihkin crystals <sage.combinat.crystals.catalog_kirillov_reshetikhin>`
* `B(\infty)` :mod:`(infinity) crystals <sage.combinat.crystals.catalog_infinity_crystals>`
* :mod:`Elementary crystals <sage.combinat.crystals.catalog_elementary_crystals>`
"""
from letters import CrystalOfLetters as Letters
from spins import CrystalOfSpins as Spins
from spins import CrystalOfSpinsPlus as SpinsPlus
from spins import CrystalOfSpinsMinus as SpinsMinus
from tensor_product import CrystalOfTableaux as Tableaux
from fast_crystals import FastCrystal as FastRankTwo
from affine import AffineCrystalFromClassical as AffineFromClassical
from affine import AffineCrystalFromClassicalAndPromotion as AffineFromClassicalAndPromotion
from affine_factorization import AffineFactorizationCrystal as AffineFactorization
from highest_weight_crystals import HighestWeightCrystal as HighestWeight
from alcove_path import CrystalOfAlcovePaths as AlcovePaths
from littelmann_path import CrystalOfLSPaths as LSPaths
from littelmann_path import CrystalOfProjectedLevelZeroLSPaths as ProjectedLevelZeroLSPaths
from kyoto_path_model import KyotoPathModel
from generalized_young_walls import CrystalOfGeneralizedYoungWalls as GeneralizedYoungWalls
from monomial_crystals import CrystalOfNakajimaMonomials as NakajimaMonomials
from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import TensorProductOfKirillovReshetikhinTableaux
from sage.combinat.rigged_configurations.rc_crystal import CrystalOfRiggedConfigurations as RiggedConfigurations

from tensor_product import TensorProductOfCrystals as TensorProduct
from direct_sum import DirectSumOfCrystals as DirectSum

def KirillovReshetikhin(cartan_type, r, s, model='KN'):
    """
    Return the Kirillov-Reshetikhin crystal of the given model.

    INPUT:

    - ``cartan_type`` -- the Cartan type
    - ``r`` -- an index of the classical Dynkin diagram
    - ``s`` -- a positive integer
    - ``model`` -- (default: ``'KN'``) can be one of the following:

      * ``'KN'`` or ``'KashiwaraNakashimaTableaux'`` - use the
        Kashiwara-Nakashima tableaux model
      * ``'KR'`` or ``'KirillovReshetkihinTableaux'`` - use the
        Kirillov-Reshetkihin tableaux model
      * ``'RC'`` or ``'RiggedConfiguration'`` - use the rigged
        configuration model
      * ``'LSPaths'`` - use the LS path model

    EXAMPLES::

        sage: KN = crystals.KirillovReshetikhin(['D',4,1], 2, 1)
        sage: KR = crystals.KirillovReshetikhin(['D',4,1], 2, 1, model='KR')
        sage: RC = crystals.KirillovReshetikhin(['D',4,1], 2, 1, model='RC')
        sage: LS = crystals.KirillovReshetikhin(['D',4,1], 2, 1, model='LSPaths')
        sage: G = KN.digraph()
        sage: G.is_isomorphic(KR.digraph(), edge_labels=True)
        True
        sage: G.is_isomorphic(RC.digraph(), edge_labels=True)
        True
        sage: G.is_isomorphic(LS.digraph(), edge_labels=True)
        True

    TESTS::

        sage: KN = crystals.KirillovReshetikhin(['D',4,1], 2, 1)
        sage: KN2 = crystals.KirillovReshetikhin(['D',4,1], 2, 1, model='KN')
        sage: KN3 = crystals.KirillovReshetikhin(['D',4,1], 2, 1, model='KashiwaraNakashimaTableaux')
        sage: KN is KN2 and KN is KN3
        True
    """
    if model in ['KN', 'KashiwaraNakashimaTableaux']:
        from sage.combinat.crystals.kirillov_reshetikhin import KirillovReshetikhinCrystal
        return KirillovReshetikhinCrystal(cartan_type, r, s)
    if model in ['KR', 'KirillovReshetikhinTableaux']:
        from sage.combinat.rigged_configurations.kr_tableaux import KirillovReshetikhinTableaux
        return KirillovReshetikhinTableaux(cartan_type, r, s)
    if model in ['RC', 'RiggedConfigurations']:
        from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
        return RiggedConfigurations(cartan_type, [[r,s]])
    if model == 'LSPaths':
        from sage.combinat.crystals.kirillov_reshetikhin import KirillovReshetikhinCrystalFromLSPaths
        return KirillovReshetikhinCrystalFromLSPaths(cartan_type, r, s)
    raise ValueError("invalid model")

import catalog_kirillov_reshetikhin as kirillov_reshetikhin
import catalog_infinity_crystals as infinity
import catalog_elementary_crystals as elementary

