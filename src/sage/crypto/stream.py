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
from sage.rings.polynomial_element import is_Polynomial
from stream_cipher import LFSRCipher, ShrinkingGeneratorCipher

from cryptosystem import SymmetricKeyCryptosystem

class LFSRCryptosystem(SymmetricKeyCryptosystem):
    """
    Linear feedback shift register cryptosystem class
    """
    def __init__(self, field = None):
        """
        Create a linear feedback shift cryptosystem.

        INPUT:
            A string monoid over a binary alphabet.

        OUTPUT:

        EXAMPLES:
            sage: E = LFSRCryptosystem(FiniteField(2))
	    sage: E
	    LFSR cryptostem over Finite Field of size 2

	TODO: Implement LFSR cryptosytem for arbitrary rings. The current
	implementation if limitated to the finite field of 2 elements only
        because of the dependence on binary strings.
        """
	if field is None:
	   field = FiniteField(2)
	if field.cardinality() != 2:
	    raise NotImplementedError, "Not yet implemented."
	S = BinaryStrings()
	P = PolynomialRing(FiniteField(2),'x')
        SymmetricKeyCryptosystem.__init__(self, S, S, None)
	self._field = field

    def __call__(self, key):
        """
        Create a LFSR cipher.

        INPUT:
            A polynomial and inital state of the LFSR.
        """
	if not isinstance(key, (list,tuple)) and len(key) == 2:
	    raise TypeError, "Argument key (= %s) must be a list of tuple of length 2" % key
	poly = key[0]; IS = key[1]
        if not is_Polynomial(poly):
            raise TypeError, "poly (= %s) must be a polynomial." % poly
        if not isinstance(IS, (list,tuple)):
            raise TypeError, "IS (= %s) must be an initial in the key space."%K
        if len(IS) != poly.degree():
            raise TypeError, \
  	        "The length of IS (= %s) must equal the degree of poly (= %s)" % (IS, poly)
        return LFSRCipher(self, poly, IS)

    def __repr__(self):
        return "LFSR cryptosystem over %s" % self._field


class ShrinkingGeneratorCryptosystem(SymmetricKeyCryptosystem):
    """
    Shrinking generator cryptosystem class
    """
    def __init__(self):
        """
        Create a shrinking generator cryptosystem.

        INPUT:
            A string monoid over a binary alphabet.

        OUTPUT:


        EXAMPLES:
            sage: E = LFSRCryptosystem()
        """
	S = BinaryStrings()
	P = PolynomialRing(FiniteField(2),'x')
        SymmetricKeyCryptosystem.__init__(self, S, S, None)
