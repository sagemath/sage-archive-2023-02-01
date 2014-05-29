r"""
Commutative rings
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

from sage.categories.category_with_axiom import CategoryWithAxiom

class CommutativeRings(CategoryWithAxiom):
    """
    The category of commutative rings

    commutative rings with unity, i.e. rings with commutative * and
    a multiplicative identity

    EXAMPLES::

         sage: C = CommutativeRings(); C
         Category of commutative rings
         sage: C.super_categories()
         [Category of rings, Category of commutative magmas]

    TESTS::

        sage: TestSuite(C).run()

        sage: QQ['x,y,z'] in CommutativeRings()
        True
        sage: GroupAlgebra(DihedralGroup(3), QQ) in CommutativeRings()
        False
        sage: MatrixSpace(QQ,2,2) in CommutativeRings()
        False

    GroupAlgebra should be fixed::

        sage: GroupAlgebra(CyclicPermutationGroup(3), QQ) in CommutativeRings() # todo: not implemented
        True

    """

    class ElementMethods:
        pass
