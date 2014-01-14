r"""
Commutative rings
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

from sage.categories.category import Category
from sage.categories.category_singleton import Category_singleton
from sage.misc.cachefunc import cached_method

class CommutativeRings(Category_singleton):
    """
    The category of commutative rings

    commutative rings with unity, i.e. rings with commutative * and
    a multiplicative identity

    EXAMPLES::

      sage: CommutativeRings()
      Category of commutative rings
      sage: CommutativeRings().super_categories()
      [Category of rings]

    TESTS::

        sage: TestSuite(CommutativeRings()).run()

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

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: CommutativeRings().super_categories()
            [Category of rings]
        """
        from sage.categories.rings import Rings
        return [Rings()]

    class ParentMethods:
        def is_commutative(self):
            """
            Return True, since commutative rings are commutative.

            EXAMPLES::

                sage: Parent(QQ,category=CommutativeRings()).is_commutative()
                True

            """
            return True

    class ElementMethods:
        pass
