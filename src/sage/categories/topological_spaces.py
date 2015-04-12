r"""
Topological Spaces
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.categories.category import Category
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory
from sage.categories.sets_cat import Sets

class TopologicalSpacesCategory(RegressiveCovariantConstructionCategory):

    _functor_category = "Topological"

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Groups().Topological()  # indirect doctest
            Category of graded algebras with basis over Rational Field
        """
        return "topological {}".format(self.base_category()._repr_object_names())

class TopologicalSpaces(TopologicalSpacesCategory):
    """
    The category of topological spaces.

    EXAMPLES::

        sage: Sets().Topological()
        Category of topological spaces
        sage: Sets().Topological().super_categories()
        [Category of modules over Integer Ring]

    The category of topological spaces defines the topological structure,
    which shall be preserved by morphisms::

        sage: Sets().Topological().additional_structure()
        Category of topological spaces

    TESTS::

        sage: TestSuite(Sets().Topological()).run()
    """
    # We must override the general object because the names don't match
    _base_category_class = (Sets,)

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Sets().Topological()  # indirect doctest
            Category of topological spaces
        """
        return "topological spaces"

