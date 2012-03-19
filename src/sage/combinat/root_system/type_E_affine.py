"""
Root system data for (untwisted) type E affine
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

            sage: ct = CartanType(['E',6,1])
            sage: ct
            ['E', 6, 1]
            sage: ct._repr_(compact = True)
            'E6~'

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
            ['E', 6]
            sage: ct.dual()
            ['E', 6, 1]

        TESTS::

            sage: ct == loads(dumps(ct))
            True
        """
        assert n >= 6 and n <= 8
        CartanType_standard_untwisted_affine.__init__(self, "E", n)

    def dynkin_diagram(self):
        """
        Returns the extended Dynkin diagram for affine type E.

        EXAMPLES::

            sage: e = CartanType(['E', 6, 1]).dynkin_diagram()
            sage: e
                    O 0
                    |
                    |
                    O 2
                    |
                    |
            O---O---O---O---O
            1   3   4   5   6
            E6~
            sage: sorted(e.edges())
            [(0, 2, 1),
             (1, 3, 1),
             (2, 0, 1),
             (2, 4, 1),
             (3, 1, 1),
             (3, 4, 1),
             (4, 2, 1),
             (4, 3, 1),
             (4, 5, 1),
             (5, 4, 1),
             (5, 6, 1),
             (6, 5, 1)]

            sage: e = CartanType(['E', 7, 1]).dynkin_diagram()
            sage: e
                        O 2
                        |
                        |
            O---O---O---O---O---O---O
            0   1   3   4   5   6   7
            E7~
            sage: sorted(e.edges())
            [(0, 1, 1), (1, 0, 1), (1, 3, 1), (2, 4, 1), (3, 1, 1), (3, 4, 1),
             (4, 2, 1), (4, 3, 1), (4, 5, 1), (5, 4, 1), (5, 6, 1),
             (6, 5, 1), (6, 7, 1), (7, 6, 1)]
            sage: e = CartanType(['E', 8, 1]).dynkin_diagram()
            sage: e
                    O 2
                    |
                    |
            O---O---O---O---O---O---O---O
            1   3   4   5   6   7   8   0
            E8~
            sage: sorted(e.edges())
            [(0, 8, 1), (1, 3, 1), (2, 4, 1), (3, 1, 1), (3, 4, 1),
             (4, 2, 1), (4, 3, 1), (4, 5, 1), (5, 4, 1), (5, 6, 1),
             (6, 5, 1), (6, 7, 1), (7, 6, 1), (7, 8, 1), (8, 0, 1), (8, 7, 1)]

        """
        from dynkin_diagram import DynkinDiagram_class
        n = self.n
        g = DynkinDiagram_class(self)
        g.add_edge(1,3)
        g.add_edge(2,4)
        for i in range(3,n):
            g.add_edge(i, i+1)
        if n == 6:
            g.add_edge(0, 2)
        elif n == 7:
            g.add_edge(0, 1)
        elif n == 8:
            g.add_edge(0, 8)
        else:
            assert False, "Invalid Cartan Type for Type E affine"
        return g



    def ascii_art(self, label = lambda x: x):
        """
        Returns a ascii art representation of the extended Dynkin diagram

        EXAMPLES::

            sage: print CartanType(['E',6,1]).ascii_art(label = lambda x: x+2)
                    O 2
                    |
                    |
                    O 4
                    |
                    |
            O---O---O---O---O
            3   5   6   7   8
            sage: print CartanType(['E',7,1]).ascii_art(label = lambda x: x+2)
                        O 4
                        |
                        |
            O---O---O---O---O---O---O
            2   3   5   6   7   8   9
            sage: print CartanType(['E',8,1]).ascii_art(label = lambda x: x+1)
                    O 3
                    |
                    |
            O---O---O---O---O---O---O---O
            2   4   5   6   7   8   9   1
        """
        n = self.n
        if n == 6:
            return "        O %s\n        |\n        |\n        O %s\n        |\n        |\nO---O---O---O---O\n%s   %s   %s   %s   %s"\
                %tuple(label(i) for i in (0,2,1,3,4,5,6))
        elif n == 7:
            return "            O %s\n            |\n            |\nO---O---O---O---O---O---O\n%s   %s   %s   %s   %s   %s   %s"\
                %tuple(label(i) for i in (2,0,1,3,4,5,6,7))
        elif n == 8:
            return "        O %s\n        |\n        |\nO---O---O---O---O---O---O---O\n%s   %s   %s   %s   %s   %s   %s   %s"\
                %tuple(label(i) for i in (2,1,3,4,5,6,7,8,0))
