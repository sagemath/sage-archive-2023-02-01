r"""
Manifolds
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import LazyImport
from sage.categories.category import Category
from sage.categories.category_types import Category_over_base_ring
from sage.categories.sets_cat import Sets
from sage.categories.fields import Fields
from sage.rings.all import RR

class Manifolds(Category_over_base_ring):
    r"""
    The category of maniolds.

    Let `k` be a field. A `d`-dimensional `k`-*manifold* `M` is a
    second countable Hausdorff space such that the neighborhood
    of any point `x \in M` is homeomorphic to `k^d`.

    INPUT:

    - ``k`` -- (default: ``RR``) the field `k`

    EXAMPLES::

        sage: from sage.categories.manifolds import Manifolds
        sage: C = Manifolds(); C
        Category of manifolds over Real Field
        sage: C.super_categories()
        [Category of additive commutative additive associative additive unital distributive magmas and additive magmas,
         Category of modules over Integer Ring]
        sage: C = Manifolds(CC); C
        Category of manifolds over Complex Field

    TESTS::

        sage: TestSuite(C).run()
    """
    def __init__(self, k=RR):
        """
        Initialize ``self``.
        """
        if k not in Fields():
            raise ValueError("must be defined over a field")
        Category_over_base_ring.__init__(self, k)

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.magmatic_algebras import MagmaticAlgebras
            sage: MagmaticAlgebras(ZZ).super_categories()
            [Category of additive commutative additive associative additive unital distributive magmas and additive magmas, Category of modules over Integer Ring]

            sage: from sage.categories.additive_semigroups import AdditiveSemigroups
            sage: MagmaticAlgebras(ZZ).is_subcategory((AdditiveSemigroups() & Magmas()).Distributive())
            True
        """
        return [Sets().Topological()]

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of manifolds defines no new
        structure: a morphism of metric spaces between manifolds
        is a manifold morphism.

        .. SEEALSO:: :meth:`Category.additional_structure`

        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: Manifolds().additional_structure()
        """
        return None

    class Connected:
        """
        The category of connected manifolds.
        """
        class ParentMethods:
            @abstract_method
            def dimension(self):
                """
                Return the dimension of ``self``.
                """

