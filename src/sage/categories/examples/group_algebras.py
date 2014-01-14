#*****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.all import GroupAlgebras
from sage.combinat.free_module import CombinatorialFreeModule

class MyGroupAlgebra(CombinatorialFreeModule):
    r"""
    An example of group algebra.

    .. TODO::

        either merge with the example of AlgebrasWithBasis and HopfAlgebrasWithBasis, or make them really different examples.
    """

    def __init__(self, R, G):
        """
        EXAMPLES::

            sage: from sage.categories.examples.group_algebras import MyGroupAlgebra
            sage: MyGroupAlgebra(ZZ,DihedralGroup(4)) # indirect doctest
            The group algebra of the Dihedral group of order 8 as a permutation group over Integer Ring
        """
        self._group = G
        CombinatorialFreeModule.__init__(self, R, G, category = GroupAlgebras(R))

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.examples.group_algebras import MyGroupAlgebra
            sage: MyGroupAlgebra(ZZ,DihedralGroup(4))
            The group algebra of the Dihedral group of order 8 as a permutation group over Integer Ring
        """
        return "The group algebra of the %s over %s"%(self._group, self.base_ring())

    @cached_method
    def one_basis(self):
        """
        Returns the unit

        EXAMPLES::

            sage: from sage.categories.examples.group_algebras import MyGroupAlgebra
            sage: A = MyGroupAlgebra(ZZ,DihedralGroup(4))
            sage: A.one_basis()
            ()
        """
        return self.group().one()

    def product_on_basis(self, g1, g2):
        """
        Return the product of two elements of the group

        EXAMPLES::

            sage: from sage.categories.examples.group_algebras import MyGroupAlgebra
            sage: A=MyGroupAlgebra(ZZ,PermutationGroup([[(1,2,3),(4,5)],[(3,4)]]))
            sage: x,y=A.group().gens()
            sage: A.product_on_basis(x,y)
            B[(1,2,3,5,4)]
        """
        return self.basis()[g1 * g2]

    def group(self):
        r"""
        Returns the underlying group of the group algebra

        EXAMPLES::

            sage: from sage.categories.examples.group_algebras import MyGroupAlgebra
            sage: A=MyGroupAlgebra(ZZ,PermutationGroup([[(1,2,3),(4,5)],[(3,4)]]))
            sage: A.group()
            Permutation Group with generators [(3,4), (1,2,3)(4,5)]
        """
        return self._group
