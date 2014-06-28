r"""
Algebra modules
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_module
from modules import Modules

class AlgebraModules(Category_module):
    """
    The category of modules over a fixed algebra $A$.

    EXAMPLES::

        sage: AlgebraModules(QQ['a'])
        Category of algebra modules over Univariate Polynomial Ring in a over Rational Field
        sage: AlgebraModules(QQ['a']).super_categories()
        [Category of modules over Univariate Polynomial Ring in a over Rational Field]

    Note: as of now, `A` is required to be commutative, ensuring that
    the categories of left and right modules are isomorphic. Feedback
    and use cases for potential generalizations to the non commutative
    case are welcome.

    """
    def __init__(self, A):
        """
        EXAMPLES::

            sage: AlgebraModules(QQ['a'])
            Category of algebra modules over Univariate Polynomial Ring in a over Rational Field
            sage: AlgebraModules(QQ['a,b']) # todo: not implemented (QQ['a,b'] should be in Algebras(QQ))
            sage: AlgebraModules(FreeAlgebra(QQ,2,'a,b'))
            Traceback (most recent call last):
            ...
            TypeError: A (=Free Algebra on 2 generators (a, b) over Rational Field) must be a commutative algebra
            sage: AlgebraModules(QQ)
            Traceback (most recent call last):
            ...
            TypeError: A (=Rational Field) must be a commutative algebra

        TESTS::

            sage: TestSuite(AlgebraModules(QQ['a'])).run()
        """
        from sage.categories.commutative_algebras import CommutativeAlgebras
        if not hasattr(A, "base_ring") or not A in CommutativeAlgebras(A.base_ring()):
            raise TypeError("A (=%s) must be a commutative algebra"%A)
        Category_module.__init__(self, A)

    @classmethod
    def an_instance(cls):
        """
        Returns an instance of this class

        EXAMPLES::

            sage: AlgebraModules.an_instance()
            Category of algebra modules over Univariate Polynomial Ring in x over Rational Field
        """
        from sage.rings.rational_field import QQ
        return cls(QQ['x'])

    def algebra(self):
        """
        EXAMPLES::

            sage: AlgebraModules(QQ[x]).algebra()
            Univariate Polynomial Ring in x over Rational Field
        """
        return self.base_ring()

    def super_categories(self):
        """
        EXAMPLES::

            sage: AlgebraModules(QQ[x]).super_categories()
            [Category of modules over Univariate Polynomial Ring in x over Rational Field]
        """
        R = self.algebra()
        return [Modules(R)]
