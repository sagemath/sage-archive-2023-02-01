"""
Classical Ciphers.
"""

#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cipher import SymmetricKeyCipher
from sage.monoids.string_monoid_element import StringMonoidElement
from sage.modules.free_module import FreeModule

class HillCipher(SymmetricKeyCipher):
    """
    Hill cipher class
    """
    def __init__(self, parent, key):
        """
        Create a Hill cipher.

        INPUT: Parent and key

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = HillCryptosystem(S,3)
            sage: E
            Hill cryptosystem on Free alphabetic string monoid on A-Z of block length 3
            sage: M = E.key_space()
            sage: A = M([[1,0,1],[0,1,1],[2,2,3]])
            sage: A
            [1 0 1]
            [0 1 1]
            [2 2 3]
            sage: e = E(A)
            sage: e
            [1 0 1]
            [0 1 1]
            [2 2 3]
            sage: e(S("LAMAISONBLANCHE"))
            JYVKSKQPELAYKPV

        TESTS::

            sage: S = AlphabeticStrings()
            sage: E = HillCryptosystem(S,3)
            sage: E == loads(dumps(E))
            True
        """
        # TODO: some type checking that the key is an invertible matrix?
        SymmetricKeyCipher.__init__(self, parent, key)

    def __eq__(self, right):
        return type(self) == type(right) and self.parent() == right.parent() and self.key() == right.key()

    def __call__(self, M):
        S = self.domain() # = plaintext_space = ciphertext_space
        if not isinstance(M, StringMonoidElement) and M.parent() == S:
            raise TypeError, "Argument M (= %s) must be a string in the plaintext space." % M
        m = self.parent().block_length()
        if len(M) % m != 0:
            raise TypeError, "The length of M (= %s) must be a multiple of %s." % (M, m )
        Alph = list(S.alphabet())
        A = self.key() # A is an m x m matrix
        R = A.parent().base_ring()
        V = FreeModule(R,m)
        Mstr = str(M)
        C = []
        for i in range(len(M)//m):
            v = V([ Alph.index(Mstr[m*i+j]) for j in range(m) ])
            C += (v * A).list()
        return S([ k.lift() for k in C ])

    def inverse(self):
        E = self.parent()
        try:
            B = E.inverse_key(self.key())
        except:
            raise ValueError, "Argument\n\n%s\n\nmust be an invertible cipher." % self
        return E(B)

class SubstitutionCipher(SymmetricKeyCipher):
    """
    Substitution cipher class
    """
    def __init__(self, parent, key):
        """
        Create a substitution cipher.

        INPUT: Parent and key

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = SubstitutionCryptosystem(S)
            sage: E
            Substitution cryptosystem on Free alphabetic string monoid on A-Z
            sage: K = S([ 25-i for i in range(26) ])
            sage: K
            ZYXWVUTSRQPONMLKJIHGFEDCBA
            sage: e = E(K)
            sage: m = S("THECATINTHEHAT")
            sage: e(m)
            GSVXZGRMGSVSZG

        TESTS::

            sage: S = AlphabeticStrings()
            sage: E = SubstitutionCryptosystem(S)
            sage: E == loads(dumps(E))
            True
        """
        SymmetricKeyCipher.__init__(self, parent, key)

    def __eq__(self, right):
        return type(self) == type(right) and self.parent() == right.parent() and self.key() == right.key()

    def __call__(self, M):
        S = self.domain() # = plaintext_space = ciphertext_space
        if not isinstance(M, StringMonoidElement) and M.parent() == S:
            raise TypeError, "Argument M (= %s) must be a string in the plaintext space." % M
        A = list(S.alphabet())
        K = str(self.key()) # K is a string, while we want the indices:
        I = [ A.index(K[i]) for i in range(len(K)) ]
        Mstr = str(M)
        return S([ I[A.index(Mstr[i])] for i in range(len(Mstr)) ])

    def inverse(self):
        E = self.parent()
        K = E.inverse_key(self.key())
        return E(K)

class TranspositionCipher(SymmetricKeyCipher):
    """
    Transition cipher class
    """
    def __init__(self, parent, key):
        """
        Create a transposition cipher.

        INPUT: Parent and key

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = TranspositionCryptosystem(S,14)
            sage: E
            Transposition cryptosystem on Free alphabetic string monoid on A-Z of block length 14
            sage: K = [ 14-i for i in range(14) ]
            sage: K
            [14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
            sage: e = E(K)
            sage: m = S("THECATINTHEHAT")
            sage: e(m)
            TAHEHTNITACEHT

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = TranspositionCryptosystem(S,15);
            sage: m = S("THECATANDTHEHAT")
            sage: G = E.key_space()
            sage: G
            Symmetric group of order 15! as a permutation group
            sage: g = G([ 3, 2, 1, 6, 5, 4, 9, 8, 7, 12, 11, 10, 15, 14, 13 ])
            sage: e = E(g)
            sage: e(m)
            EHTTACDNAEHTTAH

        TESTS::

            sage: S = AlphabeticStrings()
            sage: E = TranspositionCryptosystem(S,14)
            sage: E == loads(dumps(E))
            True
        """
        n = parent.block_length()
        if isinstance(key, list) and not len(key) == n:
            raise ValueError, "key (= %s) must have block length %s" % (key, n)
        SymmetricKeyCipher.__init__(self, parent, key)

    def __call__(self, M, mode = "ECB"):
        S = self.domain() # = plaintext_space = ciphertext_space
        if not isinstance(M, StringMonoidElement) and M.parent() == S:
            raise TypeError, "Argument M (= %s) must be a string in the plaintext space." % M
        if not mode == "ECB":
            raise NotImplementedError, "Enciphering not implemented for mode (= %s) other than 'ECB'." % mode
        g = self.key()
        N = len(M)
        m = self.parent().block_length()
        if not N%m == 0:
            raise TypeError, "Argument M (= %s) must be a string of length k*%s." % (M, m)
        Melt = M._element_list # this uses the internal structure of string monoids
        # Caution: this is parsed as an outer loop in k and an inner loop in i:
        #     for k in range(N//m):
        #         for i in range(m):
        #             S([ Melt[g(i+1)-1+k*m]
        return S([ Melt[g(i+1)-1+k*m] for k in range(N//m) for i in range(m) ])

    def inverse(self):
        E = self.parent()
        K = E.inverse_key(self.key())
        return E(K)

class VigenereCipher(SymmetricKeyCipher):
    """
    Vigenere cipher class
    """
    def __init__(self, parent, key):
        """
        Create a Vigenere cipher.

        INPUT: Parent and key

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = VigenereCryptosystem(S,11)
            sage: K = S("SHAKESPEARE")
            sage: e = E(K)
            sage: m = S("THECATINTHEHAT")
            sage: e(m)
            LOEMELXRTYIZHT

        TESTS::

            sage: S = AlphabeticStrings()
            sage: E = VigenereCryptosystem(S,11)
            sage: E == loads(dumps(E))
            True
        """
        SymmetricKeyCipher.__init__(self, parent, key)

    def __call__(self, M, mode = "ECB"):
        S = self.domain() # = plaintext_space = ciphertext_space
        if not isinstance(M, StringMonoidElement) and M.parent() == S:
            raise TypeError, "Argument M (= %s) must be a string in the plaintext space." % M
        if not mode == "ECB":
            raise NotImplementedError, "Enciphering not implemented for mode (= %s) other than 'ECB'." % mode
        K = self.key()
        m = self.parent().period()
        n = S.ngens()
        # This uses the internal structure of string monoids
        Melt = M._element_list
        Kelt = K._element_list
        return S([ (Melt[i]+Kelt[i%m])%n for i in range(len(M)) ])

    def inverse(self):
        E = self.parent()
        K = E.inverse_key(self.key())
        return E(K)




