"""
Root system data for (untwisted) type C affine
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from cartan_type import CartanType_standard_untwisted_affine
class CartanType(CartanType_standard_untwisted_affine):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['C',4,1])
            sage: ct
            ['C', 4, 1]
            sage: ct._repr_(compact = True)
            'C4~'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            False
            sage: ct.is_affine()
            True
            sage: ct.is_untwisted_affine()
            True
            sage: ct.is_crystalographic()
            True
            sage: ct.is_simply_laced()
            False
            sage: ct.classical()
            ['C', 4]
            sage: ct.dual()
            ['C', 4, 1]^*
            sage: ct.dual().is_untwisted_affine()
            False

        TESTS::
            sage: ct == loads(dumps(ct))
            True
        """
        assert n >= 1
        CartanType_standard_untwisted_affine.__init__(self, "C", n)

    def dynkin_diagram(self):
        """
        Returns the extended Dynkin diagram for affine type C.

        EXAMPLES::

            sage: c = CartanType(['C',3,1]).dynkin_diagram()
            sage: c
             O=>=O---O=<=O
             0   1   2   3
             C3~
            sage: sorted(c.edges())
            [(0, 1, 2), (1, 0, 1), (1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 2)]

        """
        n = self.n
        if n == 1:
            import cartan_type
            res = cartan_type.CartanType(["A",1,1]).dynkin_diagram()
            res._cartan_type = self
            return res
        from dynkin_diagram import DynkinDiagram_class
        g = DynkinDiagram_class(self)
        for i in range(1, n):
            g.add_edge(i, i+1)
        g.set_edge_label(n,n-1,2)
        g.add_edge(0,1,2)
        return g

    def ascii_art(self, label = lambda x: x):
        """
        Returns a ascii art representation of the extended Dynkin diagram

        EXAMPLES::

            sage: print CartanType(['C',5,1]).ascii_art(label = lambda x: x+2)
            O=>=O---O---O---O=<=O
            2   3   4   5   6   7

            sage: print CartanType(['C',3,1]).ascii_art()
            O=>=O---O=<=O
            0   1   2   3

            sage: print CartanType(['C',2,1]).ascii_art()
            O=>=O=<=O
            0   1   2

            sage: print CartanType(['C',1,1]).ascii_art()
            O<=>O
            0   1
        """
        n = self.n
        from cartan_type import CartanType
        if n == 1:
            return CartanType(["A",1,1]).ascii_art(label)
        ret = "O=>=O"+(n-2)*"---O"+"=<=O\n%s   "%label(0)
        ret += "   ".join("%s"%label(i) for i in range(1,n+1))
        return ret
