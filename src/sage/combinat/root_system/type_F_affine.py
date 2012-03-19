"""
Root system data for (untwisted) type F affine
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

            sage: ct = CartanType(['F',4,1])
            sage: ct
            ['F', 4, 1]
            sage: ct._repr_(compact = True)
            'F4~'

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
            ['F', 4]
            sage: ct.dual()
            ['F', 4, 1]^*
            sage: ct.dual().is_untwisted_affine()
            False

        TESTS::

            sage: ct == loads(dumps(ct))
            True
        """
        CartanType_standard_untwisted_affine.__init__(self, "F", 4)

    def dynkin_diagram(self):
        """
        Returns the extended Dynkin diagram for affine type F.

        EXAMPLES::

            sage: f = CartanType(['F', 4, 1]).dynkin_diagram()
            sage: f
            O---O---O=>=O---O
            0   1   2   3   4
            F4~
            sage: sorted(f.edges())
            [(0, 1, 1), (1, 0, 1), (1, 2, 1), (2, 1, 1), (2, 3, 2), (3, 2, 1), (3, 4, 1), (4, 3, 1)]

        """
        from dynkin_diagram import DynkinDiagram_class
        g = DynkinDiagram_class(self)
        for i in range(1, 4):
            g.add_edge(i, i+1)
        g.set_edge_label(2,3,2)
        g.add_edge(0, 1)
        return g

    def ascii_art(self, label = lambda x: x):
        """
        Returns a ascii art representation of the extended Dynkin diagram

        EXAMPLES::

            sage: print CartanType(['F',4,1]).ascii_art(label = lambda x: x+2)
            O---O---O=>=O---O
            2   3   4   5   6
        """
        return "O---O---O=>=O---O\n%s   %s   %s   %s   %s"%tuple(label(i) for i in (0,1,2,3,4))
