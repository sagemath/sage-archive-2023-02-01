r"""
Catalog Of Crystals

Definition of a Crystal
-----------------------

Let `C` be a CartanType with index set `I`, and `P` be
the corresponding weight lattice of the type `C`. Let `\alpha_i`
and `\alpha^{\vee}_i` denote the corresponding simple roots
and coroots respectively. Let us give the axiomatic definition
of a crystal.

A type `C` crystal `\mathcal{B}` is a non-empty set with maps
`\operatorname{wt} : \mathcal{B} \to P`,
`e_i, f_i : \mathcal{B} \to \mathcal{B} \cup \{0\}`, and
`\varepsilon_i, \varphi_i : \mathcal{B} \to \ZZ \cup \{-\infty\}`
for `i \in I` satisfying the following properties for all `i \in I`:

- `\varphi_i(b) = \varepsilon_i(b) + \langle \alpha^{\vee}_i,
  \operatorname{wt}(b) \rangle`,

- if `e_i b \in \mathcal{B}`, then:

  * `\operatorname{wt}(e_i x) = \operatorname{wt}(b) + \alpha_i`,
  * `\varepsilon_i(e_i b) = \varepsilon_i(b) - 1`,
  * `\varphi_i(e_i b) = \varphi_i(b) + 1`,

- if `f_i b \in \mathcal{B}`, then:

  * `\operatorname{wt}(f_i b) = \operatorname{wt}(b) - \alpha_i`,
  * `\varepsilon_i(f_i b) = \varepsilon_i(b) + 1`,
  * `\varphi_i(f_i b) = \varphi_i(b) - 1`,

- `f_i b^{\prime} = b` if and only if `e_i b = b^{\prime}`
  for `b, b^{\prime} \in \mathcal{B}`,

- if `\varphi_i(b) = -\infty` for `b \in \mathcal{B}`,
  then `e_i b = f_i b = 0`.

.. SEEALSO::

    - :mod:`sage.categories.crystals`
    - :mod:`sage.combinat.crystals.crystals`

Catalog
-------

This is a catalog of crystals that are currently in Sage:

* :class:`~sage.combinat.crystals.affine.AffineCrystalFromClassical`
* :class:`~sage.combinat.crystals.affine.AffineCrystalFromClassicalAndPromotion`
* :class:`AffineFactorization <sage.combinat.crystals.affine_factorization.AffineFactorizationCrystal>`
* :class:`AlcovePaths <sage.combinat.crystals.alcove_path.CrystalOfAlcovePaths>`
* :class:`FastRankTwo <sage.combinat.crystals.fast_crystals.FastCrystal>`
* :class:`GeneralizedYoungWalls
  <sage.combinat.crystals.generalized_young_walls.CrystalOfGeneralizedYoungWalls>`
* :func:`HighestWeight <sage.combinat.crystals.highest_weight_crystals.HighestWeightCrystal>`
* :func:`KirillovReshetikhin <sage.combinat.crystals.kirillov_reshetikhin.KirillovReshetikhinCrystal>`
* :class:`KyotoPathModel <sage.combinat.crystals.kyoto_path_model.KyotoPathModel>`
* :class:`Letters <sage.combinat.crystals.letters.CrystalOfLetters>`
* :class:`LSPaths <sage.combinat.crystals.littelmann_path.CrystalOfLSPaths>`
* :class:`NakajimaMonomials <sage.combinat.crystals.monomial_crystals.CrystalOfNakajimaMonomials>`
* :class:`ProjectedLevelZeroLSPaths
  <sage.combinat.crystals.littelmann_path.CrystalOfProjectedLevelZeroLSPaths>`
* :class:`RiggedConfigurations
  <sage.combinat.rigged_configurations.rc_crystal.CrystalOfRiggedConfigurations>`
* :class:`Spins <sage.combinat.crystals.spins.CrystalOfSpins>`
* :class:`SpinsPlus <sage.combinat.crystals.spins.CrystalOfSpinsPlus>`
* :class:`SpinsMinus <sage.combinat.crystals.spins.CrystalOfSpinsMinus>`
* :class:`Tableaux <sage.combinat.crystals.tensor_product.CrystalOfTableaux>`

Functorial constructions:

* :class:`DirectSum <sage.combinat.crystals.direct_sum.DirectSumOfCrystals>`
* :class:`TensorProduct <sage.combinat.crystals.tensor_product.TensorProductOfCrystals>`

Subcatalogs:

* `B(\infty)` :mod:`(infinity) crystals <sage.combinat.crystals.catalog_infinity_crystals>`
* :mod:`Elementary crystals <sage.combinat.crystals.catalog_elementary_crystals>`
* :mod:`Kirillov-Reshetihkin crystals <sage.combinat.crystals.catalog_kirillov_reshetikhin>`
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
from sage.combinat.rigged_configurations.rc_crystal import CrystalOfRiggedConfigurations as RiggedConfigurations

from tensor_product import TensorProductOfCrystals as TensorProduct
from direct_sum import DirectSumOfCrystals as DirectSum

import catalog_kirillov_reshetikhin as kirillov_reshetikhin
import catalog_infinity_crystals as infinity
import catalog_elementary_crystals as elementary

