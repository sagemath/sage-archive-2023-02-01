"""
Root system data for (untwisted) type G affine
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cartan_type import CartanType_standard_untwisted_affine
class CartanType(CartanType_standard_untwisted_affine):
    def __init__(self):
        """
        EXAMPLES::

            sage: ct = CartanType(['G',2,1])
            sage: ct
            ['G', 2, 1]
            sage: ct._repr_(compact = True)
            'G2~'

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
            ['G', 2]
            sage: ct.dual()
            ['G', 2, 1]^*
            sage: ct.dual().is_untwisted_affine()
            False

        TESTS::

            sage: ct == loads(dumps(ct))
            True
        """
        CartanType_standard_untwisted_affine.__init__(self, "G",2)

    def dynkin_diagram(self):
        """
        Returns the extended Dynkin diagram for type G.

        EXAMPLES::

            sage: g = CartanType(['G',2,1]).dynkin_diagram()
            sage: g
              3
            O=<=O---O
            1   2   0
            G2~
            sage: sorted(g.edges())
            [(0, 2, 1), (1, 2, 1), (2, 0, 1), (2, 1, 3)]
        """
        from dynkin_diagram import DynkinDiagram_class
        g = DynkinDiagram_class(self)
        g.add_edge(1, 2)
        g.set_edge_label(2,1,3)
        g.add_edge(0, 2)
        return g

    def ascii_art(self, label = lambda x: x):
        """
        Returns an ascii art representation of the Dynkin diagram

        EXAMPLES::

            sage: print CartanType(['G',2,1]).ascii_art(label = lambda x: x+2)
              3
            O=<=O---O
            3   4   2
        """
        return "  3\nO=<=O---O\n%s   %s   %s"%tuple(label(i) for i in (1,2,0))
