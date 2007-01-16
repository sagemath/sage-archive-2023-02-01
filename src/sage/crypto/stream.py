#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.monoids.string_monoid import BinaryStrings
from sage.rings.finite_field import FiniteField
from sage.rings.polynomial_ring import PolynomialRing

from cryptosystem import SymmetricKeyCryptosystem

class LFSRCryptosystem(SymmetricKeyCryptosystem):
    """
    Linear feedback shift register cryptosystem class
    """
    def __init__(self):
        """
        Create a linear feedback shift cryptosystem.

        INPUT:
            A string monoid over a binary alphabet.

        EXAMPLES:
            sage: E = LFSRCryptosystem()
        """
	S = BinaryStrings()
	P = PolynomialRing(FiniteField(2),'x')
        SymmetricKeyCryptosystem.__init__(self, S, S, None)


