r"""
Vector Spaces
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
from sage.categories.fields import Fields
from sage.categories.dual import DualObjectsCategory
from sage.misc.cachefunc import cached_method
from sage.categories.fields import Fields
_Fields = Fields()

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
    @staticmethod
    def __classcall_private__(cls, K, check=True):
        """
        INPUT:

        - `K` -- a field
        - ``check`` -- a boolean (default: True) whether to check that `K` is a field.

        EXAMPLES::

            sage: VectorSpaces(QQ) is VectorSpaces(QQ, check=False)
            True

        By default, it is checked that ``K`` is a field::

            sage: VectorSpaces(ZZ)
            Traceback (most recent call last):
            ...
            AssertionError: The base ring must be a field.

        With ``check=False``, the check is disabled, possibly enabling
        incorrect inputs::

            sage: VectorSpaces(ZZ, check=False)
            Category of vector spaces over Integer Ring
        """
        if check:
            assert K in _Fields, "The base ring must be a field."
        return super(VectorSpaces, cls).__classcall__(cls, K)

    def __init__(self, K):
        """
        EXAMPLES::

            sage: VectorSpaces(QQ)
            Category of vector spaces over Rational Field
            sage: VectorSpaces(ZZ)
            Traceback (most recent call last):
            ...
            AssertionError: The base ring must be a field.

        TESTS::

            sage: C = QQ^10      # vector space
            sage: TestSuite(C).run()
            sage: TestSuite(VectorSpaces(QQ)).run()
        """
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
        except (TypeError, AttributeError) as msg:
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

    def super_categories(self):
        """
        EXAMPLES::

            sage: VectorSpaces(QQ).super_categories()
            [Category of modules over Rational Field]
        """
        R = self.base_field()
        from sage.categories.modules import Modules
        return [Modules(R, dispatch = False)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass

    class DualObjects(DualObjectsCategory):

        def extra_super_categories(self):
            r"""
            Returns the dual category

            EXAMPLES:

            The category of algebras over the Rational Field is dual
            to the category of coalgebras over the same field::

                sage: C = VectorSpaces(QQ)
                sage: C.dual()
                Category of duals of vector spaces over Rational Field
                sage: C.dual().super_categories() # indirect doctest
                [Category of vector spaces over Rational Field]
            """
            return [self.base_category()]
