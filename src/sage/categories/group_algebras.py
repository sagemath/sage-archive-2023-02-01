r"""
GroupAlgebras
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

class GroupAlgebras(Category_over_base_ring):
    """
    EXAMPLES::

        sage: GroupAlgebras(IntegerRing())
        Category of group algebras over Integer Ring
        sage: GroupAlgebras(IntegerRing()).super_categories()
        [Category of monoid algebras over Integer Ring]

    TESTS::

        sage: C = GroupAlgebras(ZZ)
        sage: TestSuite(C).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GroupAlgebras(QQ).super_categories()
            [Category of monoid algebras over Rational Field]
        """
        from monoid_algebras import MonoidAlgebras
        R = self.base_ring()
        return [MonoidAlgebras(R)] # TODO: could become Kac algebras / Hopf algebras
