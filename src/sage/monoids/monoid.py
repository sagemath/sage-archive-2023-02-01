r"""
Base class for all monoids
"""

from sage.structure.parent import Parent
from sage.misc.cachefunc import cached_method

def is_Monoid(x):
    """
    Return True if x is of type Monoid_class.

    EXAMPLES::

        sage: from sage.monoids.monoid import is_Monoid
        sage: is_Monoid(0)
        False
        sage: is_Monoid(ZZ)   # not technical math meaning of monoid
        False
        sage: is_Monoid(sage.monoids.monoid.Monoid_class(('a','b')))
        True
        sage: F.<a,b,c,d,e> = FreeMonoid(5)
        sage: is_Monoid(F)
        True
    """
    return isinstance(x, Monoid_class)

class Monoid_class(Parent):
    def __init__(self,names):
        """
        EXAMPLES::

            sage: from sage.monoids.monoid import Monoid_class
            sage: Monoid_class(('a','b'))
            <class 'sage.monoids.monoid.Monoid_class_with_category'>
        """
        from sage.categories.monoids import Monoids
        Parent.__init__(self, base=self,names=names,category=Monoids())

    @cached_method
    def gens(self):
        """
        Return generators for self.

        EXAMPLES::

            sage: F.<a,b,c,d,e> = FreeMonoid(5)
            sage: F.gens()
            (a, b, c, d, e)
        """
        return tuple(self.gen(i) for i in range(self.ngens()))
