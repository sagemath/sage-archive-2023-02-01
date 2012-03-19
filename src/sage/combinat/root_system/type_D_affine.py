"""
Root system data for (untwisted) type D affine
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cartan_type import CartanType_standard_untwisted_affine, CartanType_simply_laced
class CartanType(CartanType_standard_untwisted_affine, CartanType_simply_laced):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['D',4,1])
            sage: ct
            ['D', 4, 1]
            sage: ct._repr_(compact = True)
            'D4~'

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
            True
            sage: ct.classical()
            ['D', 4]
            sage: ct.dual()
            ['D', 4, 1]

        TESTS::

            sage: ct == loads(dumps(ct))
            True
        """
        assert n >= 3
        CartanType_standard_untwisted_affine.__init__(self, "D", n)

    def dynkin_diagram(self):
        """
        Returns the extended Dynkin diagram for affine type D.

        EXAMPLES::

           sage: d = CartanType(['D', 6, 1]).dynkin_diagram()
           sage: d
              0 O       O 6
                |       |
                |       |
            O---O---O---O---O
            1   2   3   4   5
            D6~
           sage: sorted(d.edges())
           [(0, 2, 1), (1, 2, 1), (2, 0, 1), (2, 1, 1), (2, 3, 1),
            (3, 2, 1), (3, 4, 1), (4, 3, 1), (4, 5, 1), (4, 6, 1), (5, 4, 1), (6, 4, 1)]

           sage: d = CartanType(['D', 4, 1]).dynkin_diagram()
           sage: d
               O 4
               |
               |
           O---O---O
           1   |2  3
               |
               O 0
           D4~
           sage: sorted(d.edges())
           [(0, 2, 1),
            (1, 2, 1),
            (2, 0, 1),
            (2, 1, 1),
            (2, 3, 1),
            (2, 4, 1),
            (3, 2, 1),
            (4, 2, 1)]

           sage: d = CartanType(['D', 3, 1]).dynkin_diagram()
           sage: d
           0
           O-------+
           |       |
           |       |
           O---O---O
           3   1   2
           D3~
           sage: sorted(d.edges())
           [(0, 2, 1), (0, 3, 1), (1, 2, 1), (1, 3, 1), (2, 0, 1), (2, 1, 1), (3, 0, 1), (3, 1, 1)]

        """
        from dynkin_diagram import DynkinDiagram_class
        n = self.n
        if n == 3:
            import cartan_type
            res = cartan_type.CartanType(["A",3,1]).relabel({0:0, 1:3, 2:1, 3: 2}).dynkin_diagram()
            res._cartan_type = self
            return res
        g = DynkinDiagram_class(self)
        for i in range(1, n-1):
            g.add_edge(i, i+1)
        g.add_edge(n-2,n)
        g.add_edge(0,2)
        return g

    def ascii_art(self, label = lambda x: x):
        """
        Returns a ascii art representation of the extended Dynkin diagram

        TESTS::

            sage: print CartanType(['D',6,1]).ascii_art(label = lambda x: x+2)
              2 O       O 8
                |       |
                |       |
            O---O---O---O---O
            3   4   5   6   7

            sage: print CartanType(['D',4,1]).ascii_art(label = lambda x: x+2)
                O 6
                |
                |
            O---O---O
            3   |4  5
                |
                O 2

            sage: print CartanType(['D',3,1]).ascii_art(label = lambda x: x+2)
            2
            O-------+
            |       |
            |       |
            O---O---O
            5   3   4
        """
        n = self.n
        if n == 3:
            import cartan_type
            return cartan_type.CartanType(["A",3,1]).relabel({0:0, 1:3, 2:1, 3: 2}).ascii_art(label)
        if n == 4:
            return "    O %s\n    |\n    |\nO---O---O\n%s   |%s  %s\n    |\n    O %s"%tuple(label(i) for i in (4,1,2,3,0))
        ret = "  %s O"%label(0) +(4*(n-4)-1)*" "+"O %s"%label(n)+ "\n"
        ret += "    |"          +(4*(n-4)-1)*" "                  +"|\n"
        ret += "    |"          +(4*(n-4)-1)*" "                  +"|\n"
        ret += (n-2)*"O---"+"O\n"
        ret += "   ".join("%s"%label(i) for i in range(1,n))
        return ret
