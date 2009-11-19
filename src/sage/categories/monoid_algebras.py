r"""
MonoidAlgebras
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.misc.cachefunc import cached_method

class MonoidAlgebras(Category_over_base_ring):
    """
    The category of all monoid algebras over a given base ring.

    EXAMPLES::

        sage: MonoidAlgebras(GF(2))
        Category of monoid algebras over Finite Field of size 2

    TESTS::

        sage: TestSuite(MonoidAlgebras(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: MonoidAlgebras(QQ).super_categories()
            [Category of algebras over Rational Field]
        """
        from algebras import Algebras
        R = self.base_ring()
        return [Algebras(R)]
