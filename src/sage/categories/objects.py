r"""
Objects
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import HomCategory
from sage.categories.category_singleton import Category_singleton

#############################################################
# Generic category (default when requesting category of
# an object using misc.functional.category
#############################################################

class Objects(Category_singleton):
    """
    The category of all objects
    the basic category

    EXAMPLES::

        sage: Objects()
        Category of objects
        sage: Objects().super_categories()
        []

    TESTS::

        sage: TestSuite(Objects()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: Objects().super_categories()
            []
        """
        return []

    def __contains__(self, x):
        """
        Anything is in the category of objects.

        EXAMPLES::

            sage: int(1) in Objects()
            True
            sage: ZZ     in Objects()
            True
            sage: 2/3    in Objects()
            True
        """
        return True

    class ParentMethods:
        pass

    class ElementMethods:
        ## Will eventually contain the basic operations which are no math
        ## latex, hash, ...
        pass

    class HomCategory(HomCategory):
        def extra_super_categories(self):
            """
            This declares that any homset `Hom(A, B)` for `A` and `B`
            in the category of objects is a set.
            This more or less assumes that the category is locally small.
            See http://en.wikipedia.org/wiki/Category_(mathematics)

            EXAMPLES::

                sage: Objects().hom_category().extra_super_categories()
                [Category of sets]
            """
            from sets_cat import Sets
            return [Sets()]

        class ParentMethods:
            pass
