"""
Root system data for type H
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cartan_type import CartanType_standard_finite, CartanType_simple
class CartanType(CartanType_standard_finite, CartanType_simple):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['H',3])
            sage: ct
            ['H', 3]
            sage: ct._repr_(compact = True)
            'H3'
            sage: ct.rank()
            3

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_affine()
            False
            sage: ct.is_crystalographic()
            False
            sage: ct.is_simply_laced()
            False

        TESTS:
            sage: ct == loads(dumps(ct))
            True
        """
        assert n in [3, 4]
        CartanType_standard_finite.__init__(self, "H", n)
