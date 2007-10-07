r"""
Base class for all monoids
"""

from sage.structure.parent_gens import ParentWithGens

def is_Monoid(x):
    return isinstance(x, Monoid_class)

class Monoid_class(ParentWithGens):
    def category(self):
        # Import this locally, rather than at the top level, to avoid
        # a circular import chain
        # (category_types -> sage.algebras.all -> sage.algebras.free_algebra ->
        # sage.monoids.free_monoid -> sage.monoids.monoid ->
        # sage.categories.all -> category_types)
        import sage.categories.all
        return sage.categories.all.Monoids()
