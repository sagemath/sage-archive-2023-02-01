r"""
Bialgebras
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_types import Category_over_base_ring
from sage.categories.all import Algebras, Coalgebras

class Bialgebras(Category_over_base_ring):
    """
    The category of bialgebras

    EXAMPLES::

        sage: Bialgebras(ZZ)
        Category of bialgebras over Integer Ring
        sage: Bialgebras(ZZ).super_categories()
        [Category of algebras over Integer Ring, Category of coalgebras over Integer Ring]

    TESTS::

        sage: TestSuite(Bialgebras(ZZ)).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: Bialgebras(QQ).super_categories()
            [Category of algebras over Rational Field, Category of coalgebras over Rational Field]
        """
        R = self.base_ring()
        return [Algebras(R), Coalgebras(R)]

    def is_structure_category(self):
        r"""
        Return whether ``self`` is a structure category.

        .. SEEALSO:: :meth:`Category.is_structure_category`

        The category of bialgebras defines no new structure: a
        morphism of coalgebras and of algebras between two bialgebras
        is a bialgebra morphism.

        .. TODO:: This category should be a :class:`CategoryWithAxiom`

        EXAMPLES::

            sage: Bialgebras(QQ).is_structure_category()
            False
        """
        return False

    class ParentMethods:
        pass

    class ElementMethods:
        pass
