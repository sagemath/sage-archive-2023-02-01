"""
Root system data for (untwisted) type A affine
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cartan_type import CartanType_standard_untwisted_affine, CartanType_simply_laced
class CartanType(CartanType_standard_untwisted_affine):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['A',4,1])
            sage: ct
            ['A', 4, 1]
            sage: ct._repr_(compact = True)
            'A4~'

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
            ['A', 4]
            sage: ct.dual()
            ['A', 4, 1]

            sage: ct = CartanType(['A', 1, 1])
            sage: ct.is_simply_laced()
            False
            sage: ct.dual()
            ['A', 1, 1]

        TESTS::

            sage: ct == loads(dumps(ct))
            True
        """
        assert n >= 1
        CartanType_standard_untwisted_affine.__init__(self, "A", n)
        if n >= 2:
            self._add_abstract_superclass(CartanType_simply_laced)

    def dynkin_diagram(self):
        """
        Returns the extended Dynkin diagram for affine type A.

        EXAMPLES::

            sage: a = CartanType(['A',3,1]).dynkin_diagram()
            sage: a
             0
             O-------+
             |       |
             |       |
             O---O---O
             1   2   3
             A3~
            sage: sorted(a.edges())
            [(0, 1, 1),
             (0, 3, 1),
             (1, 0, 1),
             (1, 2, 1),
             (2, 1, 1),
             (2, 3, 1),
             (3, 0, 1),
             (3, 2, 1)]

            sage: a = DynkinDiagram(['A',1,1])
            sage: a
            O<=>O
            0   1
            A1~
            sage: sorted(a.edges())
            [(0, 1, 2), (1, 0, 2)]
        """
        from dynkin_diagram import DynkinDiagram_class
        n = self.n
        g = DynkinDiagram_class(self)

        if n == 1:
            g.add_edge(0, 1, 2)
            g.add_edge(1, 0, 2)
        else:
            for i in range(1, n):
                g.add_edge(i, i+1)
            g.add_edge(0, 1)
            g.add_edge(0, n)
        return g


    def ascii_art(self, label = lambda x: x):
        """
        Returns a ascii art representation of the extended Dynkin diagram

        EXAMPLES::

            sage: print CartanType(['A',3,1]).ascii_art()
            0
            O-------+
            |       |
            |       |
            O---O---O
            1   2   3

            sage: print CartanType(['A',5,1]).ascii_art(label = lambda x: x+2)
            2
            O---------------+
            |               |
            |               |
            O---O---O---O---O
            3   4   5   6   7

            sage: print CartanType(['A',1,1]).ascii_art()
            O<=>O
            0   1

            sage: print CartanType(['A',1,1]).ascii_art(label = lambda x: x+2)
            O<=>O
            2   3
        """
        n = self.n
        if n == 1:
            return "O<=>O\n%s   %s"%(label(0), label(1))
        ret  = "%s\nO"%label(0)+(n-2)*"----"+"---+\n|"+(n-2)*"    "+"   |\n|"+(n-2)*"    "+"   |\n"
        ret += "---".join("O"           for i in range(1,n+1)) + "\n"
        ret += "   ".join("%s"%label(i) for i in range(1,n+1))
        return ret

    def dual(self):
        """
        Type `A_1^1` is self dual despite not being simply laced

        EXAMPLES::

            sage: CartanType(['A',1,1]).dual()
            ['A', 1, 1]

        """
        return self
