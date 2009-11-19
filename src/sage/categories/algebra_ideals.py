r"""
AlgebraIdeals
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
from category_types import Category_ideal

class AlgebraIdeals(Category_ideal):
    """
    The category of two-sided ideals in a fixed algebra `A`.

    EXAMPLES::

        sage: AlgebraIdeals(FreeAlgebra(QQ,2,'a,b'))
        Category of algebra ideals in Free Algebra on 2 generators (a, b) over Rational Field

    TODO:
     - If useful, implement AlgebraLeftIdeals and AlgebraRightIdeals
       of which AlgebraIdeals would be a subcategory

     - Make AlgebraIdeals(R), return CommutativeAlgebraIdeals(R) when R is
       commutative
    """
    def __init__(self, A):
        """
        EXAMPLES::

            sage: AlgebraIdeals(FreeAlgebra(QQ,2,'a,b'))
            Category of algebra ideals in Free Algebra on 2 generators (a, b) over Rational Field
            sage: AlgebraIdeals(QQ)
            Traceback (most recent call last):
            ...
            TypeError: A (=Rational Field) must be an algebra

        TESTS::

            sage: TestSuite(AlgebraIdeals(FreeAlgebra(QQ,2,'a,b'))).run()
        """
        from sage.algebras.algebra import is_Algebra
        if not is_Algebra(A): # A not in Algebras() ?
            raise TypeError, "A (=%s) must be an algebra"%A
        Category_ideal.__init__(self, A)

    def algebra(self):
        """
        EXAMPLES::

            sage: AlgebraIdeals(QQ[x]).algebra()
            Univariate Polynomial Ring in x over Rational Field
        """
        return self.ambient()

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: AlgebraIdeals(QQ[x]).super_categories()
            [Category of algebra modules over Univariate Polynomial Ring in x over Rational Field]
        """
        from algebra_modules import AlgebraModules
        R = self.algebra()
        return [AlgebraModules(R)]
