r"""
CommutativeAlgebras
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from category_types import Category_over_base_ring
from algebras import Algebras

class CommutativeAlgebras(Category_over_base_ring):
    """
    The category of commutative algebras with unit over a given base ring.

    EXAMPLES::

        sage: M = CommutativeAlgebras(GF(19))
        sage: M
        Category of commutative algebras over Finite Field of size 19

    TESTS::

        sage: TestSuite(CommutativeAlgebras(ZZ)).run()

    Todo:

     - product   ( = cartesian product)
     - coproduct ( = tensor product over base ring)
    """


    def __contains__(self, A):
        """
        EXAMPLES::

            sage: QQ['a'] in CommutativeAlgebras(QQ)
            True
            sage: QQ['a,b'] in CommutativeAlgebras(QQ) # todo: not implemented
            True
            sage: FreeAlgebra(QQ,2,'a,b') in CommutativeAlgebras(QQ)
            False

        TODO: get rid of this method once all commutative algebras in
        Sage declare themselves in this category
        """
        return super(CommutativeAlgebras, self).__contains__(A) or \
            (A in Algebras(self.base_ring()) and hasattr(A, "is_commutative") and A.is_commutative())

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: CommutativeAlgebras(QQ).super_categories()
            [Category of algebras over Rational Field]
        """
        R = self.base_ring()
        return [Algebras(R)]
