r"""
Monoids
"""

from sage.structure.parent import Parent
from sage.misc.cachefunc import cached_method

def is_Monoid(x):
    r"""
    Returns True if ``x`` is of type ``Monoid_class``.

    EXAMPLES::

        sage: from sage.monoids.monoid import is_Monoid
        sage: is_Monoid(0)
        False
        sage: is_Monoid(ZZ)   # The technical math meaning of monoid has
        ...                   # no bearing whatsoever on the result: it's
        ...                   # a typecheck which is not satisfied by ZZ
        ...                   # since it does not inherit from Monoid_class.
        False
        sage: is_Monoid(sage.monoids.monoid.Monoid_class(('a','b')))
        True
        sage: F.<a,b,c,d,e> = FreeMonoid(5)
        sage: is_Monoid(F)
        True
    """
    return isinstance(x, Monoid_class)

class Monoid_class(Parent):
    def __init__(self, names):
        r"""
        EXAMPLES::

            sage: from sage.monoids.monoid import Monoid_class
            sage: Monoid_class(('a','b'))
            <class 'sage.monoids.monoid.Monoid_class_with_category'>

        TESTS::

            sage: F.<a,b,c,d,e> = FreeMonoid(5)
            sage: TestSuite(F).run()
        """
        from sage.categories.monoids import Monoids
        category = Monoids().FinitelyGeneratedAsMagma()
        Parent.__init__(self, base=self, names=names, category=category)

    @cached_method
    def gens(self):
        r"""
        Returns the generators for ``self``.

        EXAMPLES::

            sage: F.<a,b,c,d,e> = FreeMonoid(5)
            sage: F.gens()
            (a, b, c, d, e)
        """
        return tuple(self.gen(i) for i in range(self.ngens()))

    def monoid_generators(self):
        r"""
        Returns the generators for ``self``.

        EXAMPLES::

            sage: F.<a,b,c,d,e> = FreeMonoid(5)
            sage: F.monoid_generators()
            Family (a, b, c, d, e)
        """
        from sage.sets.family import Family
        return Family(self.gens())

