"""
Root system data for type BC affine
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cartan_type import CartanType_standard_affine
from sage.rings.integer_ring import ZZ
class CartanType(CartanType_standard_affine):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['BC',4,2])
            sage: ct
            ['BC', 4, 2]
            sage: ct._repr_(compact = True)
            'BC4~'
            sage: ct.dynkin_diagram()
            O=<=O---O---O=<=O
            0   1   2   3   4
            BC4~

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            False
            sage: ct.is_affine()
            True
            sage: ct.is_crystalographic()
            True
            sage: ct.is_simply_laced()
            False
            sage: ct.classical()
            ['C', 4]

            sage: dual = ct.dual()
            sage: dual.dynkin_diagram()
            O=<=O---O---O=<=O
            4   3   2   1   0
            BC4~ relabelled by {0: 4, 1: 3, 2: 2, 3: 1, 4: 0}

            sage: dual.special_node()
            4
            sage: dual.classical().dynkin_diagram()
            O---O---O=<=O
            3   2   1   0
            C4 relabelled by {1: 3, 2: 2, 3: 1, 4: 0}

            sage: CartanType(['BC',1,2]).dynkin_diagram()
              4
            O=<=O
            0   1
            BC1~

        TESTS::

            sage: TestSuite(ct).run()
        """
        assert n in ZZ and n >= 1
        CartanType_standard_affine.__init__(self, "BC", n, 2)

    def dynkin_diagram(self):
        """
        Returns the extended Dynkin diagram for affine type BC.

        EXAMPLES::

            sage: c = DynkinDiagram(['BC',3,2])
            sage: c
            O=<=O---O=<=O
            0   1   2   3
            BC3~
            sage: sorted(c.edges())
            [(0, 1, 1), (1, 0, 2), (1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 2)]

            sage: c = DynkinDiagram(["A", 6, 2]) # should be the same as above; did fail at some point!
            sage: c
            O=<=O---O=<=O
            0   1   2   3
            BC3~
            sage: sorted(c.edges())
            [(0, 1, 1), (1, 0, 2), (1, 2, 1), (2, 1, 1), (2, 3, 1), (3, 2, 2)]

            sage: c = DynkinDiagram(['BC',2,2])
            sage: c
            O=<=O=<=O
            0   1   2
            BC2~
            sage: sorted(c.edges())
            [(0, 1, 1), (1, 0, 2), (1, 2, 1), (2, 1, 2)]

            sage: c = DynkinDiagram(['BC',1,2])
            sage: c
              4
            O=<=O
            0   1
            BC1~
            sage: sorted(c.edges())
            [(0, 1, 1), (1, 0, 4)]

        """
        from dynkin_diagram import DynkinDiagram_class
        n = self.n
        g = DynkinDiagram_class(self)
        if n == 1:
            g.add_edge(1,0,4)
            return g
        g.add_edge(1,0,2)
        for i in range(1, n-1):
            g.add_edge(i, i+1)
        g.add_edge(n,n-1,2)
        return g

    def ascii_art(self, label = lambda x: x):
        """
        Returns a ascii art representation of the extended Dynkin diagram

        EXAMPLES::

            sage: print CartanType(['BC',2,2]).ascii_art()
            O=<=O=<=O
            0   1   2
            sage: print CartanType(['BC',3,2]).ascii_art()
            O=<=O---O=<=O
            0   1   2   3
            sage: print CartanType(['BC',5,2]).ascii_art(label = lambda x: x+2)
            O=<=O---O---O---O=<=O
            2   3   4   5   6   7

            sage: print CartanType(['BC',1,2]).ascii_art(label = lambda x: x+2)
              4
            O=<=O
            2   3
        """
        n = self.n
        if n == 1:
            return "  4\nO=<=O\n%s   %s"%(label(0), label(1))
        ret = "O=<=O"+(n-2)*"---O"+"=<=O\n%s   "%label(0)
        ret += "   ".join("%s"%label(i) for i in range(1,n+1))
        return ret

    def classical(self):
        """
        Returns the classical Cartan type associated with self

            sage: CartanType(["BC", 3, 2]).classical()
            ['C', 3]
        """
        import cartan_type
        return cartan_type.CartanType(["C", self.n])

    def dual(self):
        """
        Implements :meth:`sage.combinat.root_system.cartan_type.CartanType_abstract.dual`.

        EXAMPLES::

            sage: CartanType(["BC", 3, 2]).dual()
            ['BC', 3, 2] relabelled by {0: 3, 1: 2, 2: 1, 3: 0}

        """
        return self.relabel(lambda x: self.n-x)
