#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from lfsr import lfsr_sequence
from cipher import SymmetricKeyCipher
from sage.monoids.string_monoid_element import StringMonoidElement

class LFSRCipher(SymmetricKeyCipher):
    """
    LFSR cipher class
    """
    def __init__(self, parent, poly, IS):
        """
        Create a LFSR cipher.

        INPUT: Parent, connection polynomial, and initial state

        EXAMPLES:
            sage: FF = FiniteField(2)
	    sage: P.<x> = PolynomialRing(FF)
            sage: E = LFSRCryptosystem(FF)
	    sage: E
            LFSR cryptosystem on Finite Field of 2 elements
	    sage: IS = [ FF.random_element() for i in range(7) ]
	    sage: g = x^7 + x + 1
            sage: e = E((g,IS))
	    sage: B = BinaryStrings()
	    sage: m = B.encoding("THECATINTHEHAT")
	    sage: e(m)
	    sage: FF = FiniteField(2)
	    sage: P.<x> = PolynomialRing(FF)
	    sage: LFSR = LFSRCryptosystem(FF)
	    sage: e = LFSR((x^2+x+1,[FF(0),FF(1)]))
	    sage: B = e.domain()
	    sage: m = B.encoding("The cat in the hat.")
	    sage: e(m)
	    00111001110111101011111001001101110101011011101000011001100101101011001000000011100101101010111100000101110100111111101100000101110101111010111101000011
	    sage: m == e(e(m))
	    True
        """
	SymmetricKeyCipher.__init__(self, parent, key = (poly, IS))

    def __call__(self, M, mode = "ECB"):
	B = self.domain() # = plaintext_space = ciphertext_space
	if not isinstance(M, StringMonoidElement) and M.parent() == B:
	    raise TypeError, "Argument M (= %s) must be a string in the plaintext space." % M
	(poly, IS) = self.key()
	n = B.ngens() # two for binary strings
	N = len(M)
	Melt = M._element_list
	Kelt = lfsr_sequence(poly.list(), IS, N)
	return B([ (Melt[i]+int(Kelt[i]))%n for i in range(N) ])

    def connection_polynomial(self):
        """
	The connection polynomial defining the LFSR of the cipher.
        """
	return self.key()[0]

    def initial_state(self):
        """
	The initial state of the LFSR cipher.
        """
	return self.key()[1]

class ShrinkingGeneratorCipher(SymmetricKeyCipher):
    """
    A shrinking generator cipher.
    """

    def __init__(self, parent, e1, e2):
        """
        Create a shrinking generator cipher.

        INPUT: Parent, key stream cipher, decimation cipher

        EXAMPLES:
            sage: FF = FiniteField(2)
	    sage: P.<x> = PolynomialRing(FF)
            sage: LFSR = LFSRCryptosystem(FF)
	    sage: IS_1 = [ FF.random_element() for i in range(7) ]
            sage: e1 = LFSR((x^7 + x + 1,IS_1))
	    sage: IS_2 = [ FF.random_element() for i in range(9) ]
            sage: e2 = LFSR((x^9 + x^3 + 1,IS_2))
	    sage: E = ShrinkingGeneratorCryptosystem(FF)
	    sage: e = E((e1,e2))
	    sage: B = BinaryStrings()
	    sage: m = B.encoding("THECATINTHEHAT")
	    sage: e(m)
        """
	if not isinstance(e1, LFSRCipher):
	    raise TypeError, "Argument e1 (= %s) must be a LFSR cipher." % e1
	if not isinstance(e2, LFSRCipher):
	    raise TypeError, "Argument e2 (= %s) must be a LFSR cipher." % e2
	SymmetricKeyCipher.__init__(self, parent, key = (e1, e2))

    def keystream_cipher(self):
        """
	The LFSR cipher generating the output key stream.
        """
	return self.key()[0]

    def decimating_cipher(self):
        """
	The LFSR cipher generating the decimating key stream.
        """
	return self.key()[1]

    def __call__(self, M, mode = "ECB"):
        raise NotImplementedError, "Not yet implemented."

