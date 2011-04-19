r"""
Base class for all monoids
"""

from sage.structure.parent import Parent
from sage.misc.cachefunc import cached_method

def is_Monoid(x):
    return isinstance(x, Monoid_class)

class Monoid_class(Parent):
    def __init__(self,names):
        from sage.categories.monoids import Monoids
        #self._assign_names(names)
        Parent.__init__(self, base=self,names=names,category=Monoids())
    @cached_method
    def gens(self):
        return tuple(self.gen(i) for i in range(self.ngens()))
#    def category(self):
#        # Defining a category method is deprecated for parents.
#        # Instead, the category should be specified in the constructor.
#        # See: http://sagetrac.org/sage_trac/wiki/CategoriesRoadMap
#        if self._is_category_initialized():
#            return Parent.category(self)
#
#        # Import this locally, rather than at the top level, to avoid
#        # a circular import chain
#        # (category_types -> sage.algebras.all -> sage.algebras.free_algebra ->
#        # sage.monoids.free_monoid -> sage.monoids.monoid ->
#        # sage.categories.all -> category_types)
