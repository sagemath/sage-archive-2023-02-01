"""
Root system data for (untwisted) type B affine
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

            sage: ct = CartanType(['B',4,1])
            sage: ct
            ['B', 4, 1]
            sage: ct._repr_(compact = True)
            'B4~'

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
            ['B', 4]
            sage: ct.dual()
            ['B', 4, 1]^*
            sage: ct.dual().is_untwisted_affine()
            False

        TESTS:
            sage: ct == loads(dumps(ct))
            True
        """
        assert n >= 1
        CartanType_standard_untwisted_affine.__init__(self, "B", n)

    def dynkin_diagram(self):
        """
        Returns the extended Dynkin diagram for affine type B.

        EXAMPLES::

            sage: b = CartanType(['B',3,1]).dynkin_diagram()
            sage: b
                O 0
                |
                |
            O---O=>=O
            1   2   3
            B3~
            sage: sorted(b.edges())
            [(0, 2, 1), (1, 2, 1), (2, 0, 1), (2, 1, 1), (2, 3, 2), (3, 2, 1)]

            sage: b = CartanType(['B',2,1]).dynkin_diagram(); b
            O=>=O=<=O
            0   2   1
            B2~
            sage: sorted(b.edges())
            [(0, 2, 2), (1, 2, 2), (2, 0, 1), (2, 1, 1)]

            sage: b = CartanType(['B',1,1]).dynkin_diagram(); b
            O<=>O
            0   1
            B1~
            sage: sorted(b.edges())
            [(0, 1, 2), (1, 0, 2)]

        """
        import cartan_type
        n = self.n
        if n == 1:
            res = cartan_type.CartanType(["A",1,1]).dynkin_diagram()
            res._cartan_type = self
            return res
        if n == 2:
            res = cartan_type.CartanType(["C",2,1]).relabel({0:0, 1:2, 2:1}).dynkin_diagram()
            res._cartan_type = self
            return res
        from dynkin_diagram import DynkinDiagram_class
        g = DynkinDiagram_class(self)
        for i in range(1, n):
            g.add_edge(i, i+1)
        g.set_edge_label(n-1, n, 2)
        g.add_edge(0,2)
        return g

    def ascii_art(self, label = lambda x: x):
        """
        Returns a ascii art representation of the extended Dynkin diagram

        EXAMPLES::

            sage: print CartanType(['B',3,1]).ascii_art()
                O 0
                |
                |
            O---O=>=O
            1   2   3

            sage: print CartanType(['B',5,1]).ascii_art(label = lambda x: x+2)
                O 2
                |
                |
            O---O---O---O=>=O
            3   4   5   6   7

            sage: print CartanType(['B',2,1]).ascii_art(label = lambda x: x+2)
            O=>=O=<=O
            2   4   3
            sage: print CartanType(['B',1,1]).ascii_art(label = lambda x: x+2)
            O<=>O
            2   3
        """
        n = self.n
        from cartan_type import CartanType
        if n == 1:
            return CartanType(["A",1,1]).ascii_art(label)
        if n == 2:
            return CartanType(["C",2,1]).relabel({0:0, 1:2, 2:1}).ascii_art(label)
        ret  = "    O %s\n    |\n    |\n"%label(0)
        ret += (n-2)*"O---" + "O=>=O\n"
        ret += "   ".join("%s"%label(i) for i in range(1,n+1))
        return ret
