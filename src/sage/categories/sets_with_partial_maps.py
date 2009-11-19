r"""
SetsWithPartialMaps
"""
#*****************************************************************************
#  Copyright (C) 2008 David Kohel <kohel@maths.usyd.edu> and
#                     William Stein <wstein@math.ucsd.edu>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.misc.cachefunc import cached_method

class SetsWithPartialMaps(Category):
    """
    The category whose objects are sets and whose morphisms are
    maps that are allowed to raise a ValueError on some inputs.

    This category is equivalent to the category of pointed sets,
    via the equivalence sending an object X to X union {error},
    a morphism f to the morphism of pointed sets that sends x
    to f(x) if f does not raise an error on x, or to error if it
    does.

    EXAMPLES::

        sage: SetsWithPartialMaps()
        Category with objects Sets and morphisms partially defined maps

        sage: SetsWithPartialMaps().super_categories()
        [Category of sets]

    TESTS::

        sage: TestSuite(SetsWithPartialMaps()).run()
    """
    #def __call__(self, X, pt):
    #    import sage.sets.all
    #    return sage.sets.all.Set(X, pt)

    def __repr__(self):
        """
        EXAMPLES::

            sage: SetsWithPartialMaps()
            Category with objects Sets and morphisms partially defined maps

        FIXME: is this the desired result?
        """
        return "Category with objects Sets and morphisms partially defined maps"

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: SetsWithPartialMaps().super_categories()
            [Category of sets]
        """
        from sets_cat import Sets
        return [Sets()]
