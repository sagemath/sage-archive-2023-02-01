"""
Catalog of Crystals

This is a catalog of crystals that are currently in Sage:

* :class:`Letters <sage.combinat.crystals.letters.CrystalOfLetters>`
* :class:`Tableaux <sage.combinat.crystals.tensor_product.CrystalOfTableaux>`
* :func:`HighestWeight <sage.combinat.crystals.highest_weight_crystals.HighestWeightCrystal>`
* :class:`NakajimaMonomials <sage.combinat.crystals.monomial_crystals.CrystalOfNakajimaMonomials>`
* :class:`GeneralizedYoungWalls <sage.combinat.crystals.generalized_young_walls.CrystalOfGeneralizedYoungWalls>`
* :class:`KyotoPathModel <sage.combinat.crystals.kyoto_path_model.KyotoPathModel>`
* :class:`Spins <sage.combinat.crystals.spins.CrystalOfSpins>`
* :class:`SpinsPlus <sage.combinat.crystals.spins.CrystalOfSpinsPlus>`
* :class:`SpinsMinus <sage.combinat.crystals.spins.CrystalOfSpinsMinus>`
* :class:`FastRankTwo <sage.combinat.crystals.fast_crystals.FastCrystal>`
* :class:`~sage.combinat.crystals.affine.AffineCrystalFromClassical`
* :class:`~sage.combinat.crystals.affine.AffineCrystalFromClassicalAndPromotion`
* :func:`KirillovReshetikhin <sage.combinat.crystals.kirillov_reshetikhin.KirillovReshetikhinCrystal>`

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
from sage.combinat.crystals.kirillov_reshetikhin import KirillovReshetikhinCrystal as KirillovReshetikhin

from tensor_product import TensorProductOfCrystals as TensorProduct
from direct_sum import DirectSumOfCrystals as DirectSum

import catalog_kirillov_reshetikhin as kirillov_reshetikhin
import catalog_infinity_crystals as infinity
import catalog_elementary_crystals as elementary

