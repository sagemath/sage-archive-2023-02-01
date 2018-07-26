"""
Root system data for type Q
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Daniel Bump
#       Copyright (C) 2008-2009 Justin Walker
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from __future__ import absolute_import

from .type_A import AmbientSpace
from .type_A import CartanType
from .cartan_type import CartanType_standard_finite

class CartanType(CartanType):
    """
    Cartan Type `Q_n`

    .. SEEALSO:: :func:`~sage.combinat.root_systems.cartan_type.CartanType`
    """
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['Q',4])
            sage: ct
            ['Q', 4]
            sage: ct._repr_(compact = True)
            'Q4'

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_affine()
            False
            sage: ct.is_crystallographic()
            True
            sage: ct.is_simply_laced()
            True
            sage: ct.affine()
            ['Q', 4, 1]
            sage: ct.dual()
            ['Q', 4]

        TESTS::

            sage: TestSuite(ct).run()
        """
        assert n >= 0
        CartanType_standard_finite.__init__(self, "Q", n)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CartanType(['Q',4]))
            Q_{4}
        """
        return "Q_{%s}"%self.n

    AmbientSpace = AmbientSpace

# For unpickling backward compatibility (Sage <= 4.1)
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.root_system.type_Q', 'ambient_space',  AmbientSpace)
