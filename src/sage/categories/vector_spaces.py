r"""
VectorSpaces
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_module
from sage.categories.all import Fields, Modules
from sage.misc.cachefunc import cached_method

class VectorSpaces(Category_module):
    """
    The category of (abstract) vector spaces over a given field

    ??? with an embedding in an ambient vector space ???

    EXAMPLES::

        sage: VectorSpaces(QQ)
        Category of vector spaces over Rational Field
        sage: VectorSpaces(QQ).super_categories()
        [Category of modules over Rational Field]
    """
    def __init__(self, K):
        """
        EXAMPLES::

            sage: VectorSpaces(QQ)
            Category of vector spaces over Rational Field
            sage: VectorSpaces(ZZ)
            Traceback (most recent call last):
            ...
              assert K in Fields()
            AssertionError

        TESTS::

            sage: C = QQ^10      # vector space
            sage: TestSuite(C).run()
            sage: TestSuite(VectorSpaces(QQ)).run()
        """
        assert K in Fields()
        Category_module.__init__(self, K)

    def __call__(self, x):
        """
        Try to coerce ``x`` into an object of this category

        EXAMPLES::

            sage: VectorSpaces(QQ)(ZZ^3)
            Vector space of dimension 3 over Rational Field

        """
        try:
            V = x.vector_space(self.base_field())
            if V.base_field() != self.base_field():
                V = V.change_ring(self.base_field())
        except (TypeError, AttributeError), msg:
            raise TypeError, "%s\nunable to coerce x (=%s) into %s"%(msg,x,self)
        return V

    def base_field(self):
        """
        Returns the base field over which the vector spaces of this
        category are all defined.

        EXAMPLES::

            sage: VectorSpaces(QQ).base_field()
            Rational Field
        """
        return self.base_ring()

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: VectorSpaces(QQ).super_categories()
            [Category of modules over Rational Field]
        """
        R = self.base_field()
        return [Modules(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
