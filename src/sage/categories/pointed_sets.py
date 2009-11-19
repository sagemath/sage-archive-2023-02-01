r"""
PointedSets
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

class PointedSets(Category):
    """
    The category of pointed sets.

    EXAMPLES::

        sage: PointedSets()
        Category of pointed sets

    TESTS::

        sage: TestSuite(PointedSets()).run()
    """
    #def __call__(self, X, pt):
    #    import sage.sets.all
    #    return sage.sets.all.Set(X, pt)

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: PointedSets().super_categories()
            [Category of sets]
        """
        from sets_cat import Sets
        return [Sets()] # ???
