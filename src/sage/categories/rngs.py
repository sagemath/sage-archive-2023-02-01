r"""
Rngs
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2012 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.misc.lazy_import import LazyImport
from sage.categories.distributive_magmas_and_additive_magmas import DistributiveMagmasAndAdditiveMagmas

class Rngs(CategoryWithAxiom):
    """
    The category of rngs

    A rng `(S,+,*)` is similar to a ring but not necessarilly
    unital. In other words, it is a combination of a commutative
    additive group `(S,+)` and a multiplicative semigroup `(S,*)`,
    where `*` distributes over `+`.

    EXAMPLES::

      sage: C = Rngs(); C
      Category of rngs
      sage: sorted(C.super_categories(), key=str)
      [Category of additive commutative additive associative additive inverse additive unital distributive magmas and additive magmas,
       Category of associative additive commutative additive associative additive unital distributive magmas and additive magmas]
      sage: sorted(C.axioms())
      ['AdditiveAssociative', 'AdditiveCommutative', 'AdditiveInverse', 'AdditiveUnital', 'Associative']
      sage: C is DistributiveMagmasAndAdditiveMagmas().Associative().AdditiveAssociative().AdditiveCommutative().AdditiveUnital().AdditiveInverse()
      True
      sage: C.Unital()
      Category of rings

    TESTS::

        sage: TestSuite(C).run()
    """

    _base_category_class_and_axiom = (DistributiveMagmasAndAdditiveMagmas.AdditiveAssociative.AdditiveCommutative.AdditiveUnital.AdditiveInverse, "Associative")

    Unital = LazyImport('sage.categories.rings', 'Rings', at_startup=True)
