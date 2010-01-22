r"""
Base class for all monoids
"""

from sage.structure.parent_gens import ParentWithGens

def is_Monoid(x):
    return isinstance(x, Monoid_class)

class Monoid_class(ParentWithGens):
    def category(self):
        # Defining a category method is deprecated for parents.
        # Instead, the category should be specified in the constructor.
        # See: http://sagetrac.org/sage_trac/wiki/CategoriesRoadMap
        if self._is_category_initialized():
            return Parent.category(self)

        # Import this locally, rather than at the top level, to avoid
        # a circular import chain
        # (category_types -> sage.algebras.all -> sage.algebras.free_algebra ->
        # sage.monoids.free_monoid -> sage.monoids.monoid ->
        # sage.categories.all -> category_types)
        from sage.categories.monoids import Monoids
        return Monoids()
