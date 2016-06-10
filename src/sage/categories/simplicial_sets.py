"""
Simplicial Sets
"""
#*****************************************************************************
#  Copyright (C) 2015 John H. Palmieri <palmieri at math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.sets_cat import Sets

class SimplicialSets(Category_singleton):
    r"""
    The category of simplicial sets.

    A simplicial set `X` is a collection of sets `X_i`, indexed by
    the non-negative integers, together with maps

    .. math::

        d_i: X_n \to X_{n-1}, \quad 0 \leq i \leq n \quad \text{(face maps)} \\
        s_j: X_n \to X_{n+1}, \quad 0 \leq j \leq n \quad \text{(degeneracy maps)}

    satisfying the *simplicial identities*:

    .. math::

        d_i d_j &= d_{j-1} d_i \quad \text{if } i<j \\
        d_i s_j &= s_{j-1} d_i \quad \text{if } i<j \\
        d_j s_j &= 1 = d_{j+1} s_j \\
        d_i s_j &= s_{j} d_{i-1} \quad \text{if } i>j+1 \\
        s_i s_j &= s_{j+1} s_{i} \quad \text{if } i \leq j

    Morphisms are sequences of maps `f_i : X_i \to Y_i` which commute
    with the face and degeneracy maps.

    EXAMPLES::

        sage: from sage.categories.simplicial_sets import SimplicialSets
        sage: C = SimplicialSets(); C
        Category of simplicial sets

    TESTS::

        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.simplicial_sets import SimplicialSets
            sage: SimplicialSets().super_categories()
            [Category of sets]
        """
        return [Sets()]

    class Finite(CategoryWithAxiom):
        """
        Category of finite simplicial sets.

        That is, the simplicial sets with finitely many non-degenerate
        simplices.
        """
        class ParentMethods:
            pass

        class ElementMethods:
            pass

        class MorphismMethods:
            pass
