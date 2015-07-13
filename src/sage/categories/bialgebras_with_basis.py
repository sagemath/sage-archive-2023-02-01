r"""
Bialgebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#  Copyright (C) 2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring

class BialgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The category of bialgebras with a distinguished basis.

    EXAMPLES::

        sage: C = BialgebrasWithBasis(QQ); C
        Category of bialgebras with basis over Rational Field

        sage: sorted(C.super_categories(), key=str)
        [Category of algebras with basis over Rational Field,
         Category of bialgebras over Rational Field,
         Category of coalgebras with basis over Rational Field]

    TESTS::

        sage: TestSuite(BialgebrasWithBasis(ZZ)).run()
    """
    class ParentMethods:
        def foo(self):
            return 'bar'

    class ElementMethods:
        def bar(self):
            return 'foo fighters'

