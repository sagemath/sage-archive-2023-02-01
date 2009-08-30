"""
Stream Cryptosystems
"""

#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.monoids.string_monoid import BinaryStrings
from sage.rings.finite_field import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_element import is_Polynomial
from stream_cipher import LFSRCipher, ShrinkingGeneratorCipher

from cryptosystem import SymmetricKeyCryptosystem

class LFSRCryptosystem(SymmetricKeyCryptosystem):
    """
    Linear feedback shift register cryptosystem class
    """
    def __init__(self, field = None):
        """
        Create a linear feedback shift cryptosystem.

        INPUT: A string monoid over a binary alphabet.

        OUTPUT:

        EXAMPLES::

            sage: E = LFSRCryptosystem(FiniteField(2))
            sage: E
            LFSR cryptosystem over Finite Field of size 2

        TESTS::

            sage: E = LFSRCryptosystem(FiniteField(2))
            sage: E == loads(dumps(E))
            True

        TODO: Implement LFSR cryptosystem for arbitrary rings. The current
        implementation is limited to the finite field of 2 elements only
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

    def __eq__(self,right):
        return type(self) == type(right) and self._field == right._field

    def __call__(self, key):
        """
        Create a LFSR cipher.

        INPUT: A polynomial and initial state of the LFSR.
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

    def _repr_(self):
        r"""
        Return the string representation of this LFSR cryptosystem.

        EXAMPLES::

            sage: LFSRCryptosystem(FiniteField(2))
            LFSR cryptosystem over Finite Field of size 2
        """
        return "LFSR cryptosystem over %s" % self._field

    def encoding(self,M):
        S = self.cipher_domain()
        try:
            return S.encoding(M)
        except:
            raise TypeError, "Argument M = %s does not encode in the cipher domain" % M

class ShrinkingGeneratorCryptosystem(SymmetricKeyCryptosystem):
    """
    Shrinking generator cryptosystem class
    """
    def __init__(self, field = None):
        """
        Create a shrinking generator cryptosystem.

        INPUT: A string monoid over a binary alphabet.

        OUTPUT:

        EXAMPLES::

            sage: E = ShrinkingGeneratorCryptosystem()
            sage: E
            Shrinking generator cryptosystem over Finite Field of size 2
        """
        if field is None:
           field = FiniteField(2)
        if field.cardinality() != 2:
            raise NotImplementedError, "Not yet implemented."
        S = BinaryStrings()
        P = PolynomialRing(field, 'x')
        SymmetricKeyCryptosystem.__init__(self, S, S, None)
        self._field = field

    def __call__(self, key):
        """
        Create a Shrinking generator cipher.

        INPUT: A list or tuple consisting of two LFSR ciphers (e1,e2).

        OUTPUT: The shrinking generator cipher with key stream generator e1
        and decimating cipher e2.
        """
        if not isinstance(key, (list,tuple)) and len(key) == 2:
            raise TypeError, "Argument key (= %s) must be a list of tuple of length 2" % key
        e1 = key[0]; e2 = key[1]
        if not isinstance(e1, LFSRCipher) or not isinstance(e2, LFSRCipher):
            raise TypeError, "The key (= (%s,%s)) must be a tuple of two LFSR ciphers." % key
        return ShrinkingGeneratorCipher(self, e1, e2)

    def _repr_(self):
        r"""
        Return the string representation of this shrinking generator
        cryptosystem.

        EXAMPLES::

            sage: ShrinkingGeneratorCryptosystem()
            Shrinking generator cryptosystem over Finite Field of size 2
        """
        return "Shrinking generator cryptosystem over %s" % self._field

    def encoding(self,M):
        S = self.cipher_domain()
        try:
            return S.encoding(M)
        except:
            raise TypeError, "Argument M = %s does not encode in the cipher domain" % M
