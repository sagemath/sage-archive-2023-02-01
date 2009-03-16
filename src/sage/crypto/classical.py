"""
Classical Cryptosystems
"""

#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.monoids.string_monoid import StringMonoid_class, AlphabeticStringMonoid
from sage.monoids.string_monoid_element import StringMonoidElement
from sage.monoids.string_ops import strip_encoding
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.integer_mod_ring import IntegerModRing
from sage.rings.arith import xgcd
from random import randint
from sage.matrix.matrix_space import MatrixSpace

from cryptosystem import SymmetricKeyCryptosystem
from classical_cipher import (
     HillCipher,
     SubstitutionCipher,
     TranspositionCipher,
     VigenereCipher)

class HillCryptosystem(SymmetricKeyCryptosystem):
    """
    Hill cryptosystem class
    """

    def __init__(self, S, m):
        """
        Create a Hill cryptosystem defined by the m x m matrix space
        over `\mathbb{Z} / N \mathbb{Z}`, where `N` is the alphabet size of
        the string monoid ``S``.

        INPUT:

        - ``S`` - a string monoid over some alphabet

        - ``m`` - integer `> 0`; the block length of matrices that specify
          block permutations

        OUTPUT:

        - A Hill cryptosystem of block length ``m`` over the alphabet ``S``.

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = HillCryptosystem(S,3)
            sage: E
            Hill cryptosystem on Free alphabetic string monoid on A-Z of block length 3
            sage: R = IntegerModRing(26)
            sage: M = MatrixSpace(R,3,3)
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
        if not isinstance(S, StringMonoid_class):
            raise TypeError, "S (= %s) must be a string monoid."%S
	R = IntegerModRing(S.ngens())
	M = MatrixSpace(R,m,m)
        SymmetricKeyCryptosystem.__init__(self, S, S, M, block_length = m)

    def __call__(self, A):
        """
        Create a Hill cipher.

        INPUT:

        - ``A`` - a matrix which specifies a block permutation

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
            sage: m = S("LAMAISONBLANCHE")
            sage: e(m)
            JYVKSKQPELAYKPV
            sage: c = e.inverse()
            sage: c(e(m))
            LAMAISONBLANCHE
        """
        M = self.key_space()
        m = self.block_length()
        if isinstance(A, list):
            try:
                A = M(A)
            except:
                raise TypeError, "A (= %s) must specify a square matrix of degree %s." % (A, m)
        return HillCipher(self, A)

    def _repr_(self):
        """
        Return a string representation of self.

        EXAMPLE::

            sage: A = AlphabeticStrings()
            sage: H = HillCryptosystem(A, 3)
            sage: H
            Hill cryptosystem on Free alphabetic string monoid on A-Z of block length 3
            sage: H._repr_()
            'Hill cryptosystem on Free alphabetic string monoid on A-Z of block length 3'
        """
        return "Hill cryptosystem on %s of block length %s" % (
            self.cipher_domain(), self.block_length())

    def block_length(self):
        """
        The row or column dimension of a matrix specifying a block
        permutation. Encryption and decryption keys of a Hill cipher are
        square matrices, i.e. the row and column dimensions of an encryption
        or decryption key are the same. This row/column dimension is referred
        to as the block length.

        OUTPUT:

        - The block length of an encryption/decryption key.

        EXAMPLES::

            sage: A = AlphabeticStrings()
            sage: n = randint(1, A.ngens() - 1)
            sage: H = HillCryptosystem(A, n)
            sage: H.block_length() == n
            True
        """
        return self.key_space().nrows()

    def random_key(self):
        """
        A random key within the key space of this Hill cipher. That is,
        generate a random m x m matrix to be used as a block
        permutation, where `m` is the block length of this Hill cipher. If
        `n` is the size of the cryptosystem alphabet, then there is a total
        of `n^{m^2}` possible keys. However the number of valid keys,
        i.e. invertible m x m square matrices, is smaller than
        `n^{m^2}`.

        OUTPUT:

        - A random key within the key space of this Hill cipher.

        EXAMPLE::

            sage: A = AlphabeticStrings()
            sage: n = 3
            sage: H = HillCryptosystem(A, n)
            sage: K = H.random_key()
            sage: Ki = H.inverse_key(K)
            sage: M = "LAMAISONBLANCHE"
            sage: e = H(K)
            sage: d = H(Ki)
            sage: d(e(A(M))) == A(M)
            True
        """
        M = self.key_space()
	R = M.base_ring()
        m = M.nrows()
        N = Integer(self.cipher_domain().ngens())
	while True:
	    A = M([ randint(0,N-1) for i in range(m**2) ])
	    if N.gcd(A.det()) == 1:
	        break
        return A

    def inverse_key(self, A):
        """
        The inverse key corresponding to the key ``A``.

        INPUT:

        - ``A`` - an invertible matrix of the key space of this Hill cipher

        OUTPUT:

        - The inverse matrix of ``A``.

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = HillCryptosystem(S,3)
            sage: A = E.random_key()
            sage: B = E.inverse_key(A)
            sage: M = S("LAMAISONBLANCHE")
            sage: e = E(A)
            sage: c = E(B)
            sage: c(e(M))
            LAMAISONBLANCHE
        """
	S = self.plaintext_space()
	M = self.key_space()
	if not A in M:
	    raise TypeError, "A (= %s) must be a matrix in the key space of %s." % (A, self)
	m = self.block_length()
	MatZZ = MatrixSpace(ZZ,m)
	AZ = MatZZ([ [ A[i,j].lift() for j in range(m) ] for i in range(m) ])
	AZ_adj = AZ.adjoint()
	u, r, s = xgcd(A.det().lift(),S.ngens())
	if u != 1:
	    raise ValueError, "Argument:\n\n%s\n\nis not invertible."%(A)
        return r * A.parent()(AZ_adj)

    def encoding(self, M):
        """
        The encoding of the string ``M`` over the string monoid of this
        Hill cipher. For example, if the string monoid of this Hill cipher
        is :class:`AlphabeticStringMonoid`, then the encoding of ``M`` would
        be its upper-case equivalent stripped of all non-alphabetic
        characters.

        INPUT:

        - ``M`` - a string, possibly empty

        OUTPUT:

        - The encoding of ``M`` over the string monoid of this Hill cipher.

        EXAMPLES::

            sage: M = "The matrix cipher by Lester S. Hill."
            sage: A = AlphabeticStrings()
            sage: H = HillCryptosystem(A, 7)
            sage: H.encoding(M) == A.encoding(M)
            True
        """
        S = self.cipher_domain()
        if isinstance(S,AlphabeticStringMonoid):
            return S(strip_encoding(M))
        try:
            return S.encoding(M)
        except:
            raise TypeError, "Argument M = %s does not encode in the cipher domain" % M

    def deciphering(self, A, C):
        """
        Decrypt the ciphertext ``C`` using the key ``A``.

        INPUT:

        - ``A`` - a key within the key space of this Hill cipher

        - ``C`` - a string (possibly empty) over the string monoid of this
          Hill cipher

        OUTPUT:

        - The plaintext corresponding to the ciphertext ``C``.

        EXAMPLES::

            sage: H = HillCryptosystem(AlphabeticStrings(), 3)
            sage: K = H.random_key()
            sage: M = H.encoding("Good day, mate! How ya going?")
            sage: H.deciphering(K, H.enciphering(K, M)) == M
            True
        """
        # TODO: some type checking that A is invertible hence a valid key
        i = self(self.inverse_key(A))
        return i(C)

    def enciphering(self, A, M):
        """
        Encrypt the plaintext ``M`` using the key ``A``.

        INPUT:

        - ``A`` - a key within the key space of this Hill cipher

        - ``M`` - a string (possibly empty) over the string monoid of this
          Hill cipher.

        OUTPUT:

        - The ciphertext corresponding to the plaintext ``M``.

        EXAMPLES::

            sage: H = HillCryptosystem(AlphabeticStrings(), 3)
            sage: K = H.random_key()
            sage: M = H.encoding("Good day, mate! How ya going?")
            sage: H.deciphering(K, H.enciphering(K, M)) == M
            True
        """
        # TODO: some type checking that A is invertible hence a valid key
        e = self(A)
        return e(M)

class SubstitutionCryptosystem(SymmetricKeyCryptosystem):
    """
    Substitution cryptosystem class
    """
    def __init__(self, S):
        """
        Create a substitution cryptosystem.

        INPUT:

        - ``S`` - a string monoid over some alphabet

        OUTPUT:

        - A substitution cryptosystem over the alphabet ``S``.

        EXAMPLES::

            sage: M = AlphabeticStrings()
            sage: E = SubstitutionCryptosystem(M)
            sage: E
            Substitution cryptosystem on Free alphabetic string monoid on A-Z
            sage: K = M([ 25-i for i in range(26) ])
            sage: K
            ZYXWVUTSRQPONMLKJIHGFEDCBA
            sage: e = E(K)
            sage: m = M("THECATINTHEHAT")
            sage: e(m)
            GSVXZGRMGSVSZG

        TESTS::

            sage: M = AlphabeticStrings()
            sage: E = SubstitutionCryptosystem(M)
            sage: E == loads(dumps(E))
            True
        """
        if not isinstance(S, StringMonoid_class):
            raise TypeError, "S (= %s) must be a string monoid."%S
        SymmetricKeyCryptosystem.__init__(self, S, S, S)

    def __call__(self, K):
        """
        Create a substitution cipher.

        INPUT:

        - ``K`` - a key which is a permutation of the cryptosystem alphabet

        EXAMPLES::

            sage: M = AlphabeticStrings()
            sage: E = SubstitutionCryptosystem(M)
            sage: E
            Substitution cryptosystem on Free alphabetic string monoid on A-Z
            sage: K = M([ 25-i for i in range(26) ])
            sage: K
            ZYXWVUTSRQPONMLKJIHGFEDCBA
            sage: e = E(K)
            sage: m = M("THECATINTHEHAT")
            sage: e(m)
            GSVXZGRMGSVSZG
        """
        if not isinstance(K, StringMonoidElement):
            raise TypeError, "K (= %s) must be a string."%K
        if K.parent() != self.key_space():
            raise TypeError, "K (= %s) must be a string in the key space."%K
        return SubstitutionCipher(self, K)

    def _repr_(self):
        """
        Return a string representation of self.

        EXAMPLE::

            sage: A = AlphabeticStrings()
            sage: S = SubstitutionCryptosystem(A)
            sage: S
            Substitution cryptosystem on Free alphabetic string monoid on A-Z
            sage: S._repr_()
            'Substitution cryptosystem on Free alphabetic string monoid on A-Z'
        """
        return "Substitution cryptosystem on %s" % self.cipher_domain()

    def random_key(self):
        """
        Generate a random key within the key space of this substitution
        cipher. The generated key is a permutation of the cryptosystem
        alphabet. Let `n` be the length of the alphabet. Then there is a
        total of `n!` possible keys in the key space.

        OUTPUT:

        - A random key within the key space of this cryptosystem.

        EXAMPLE::

            sage: A = AlphabeticStrings()
            sage: S = SubstitutionCryptosystem(A)
            sage: K = S.random_key()
            sage: Ki = S.inverse_key(K)
            sage: M = "THECATINTHEHAT"
            sage: e = S(K)
            sage: d = S(Ki)
            sage: d(e(A(M))) == A(M)
            True
        """
        S = self.cipher_domain()
        n = S.ngens()
        I = SymmetricGroup(n).random_element().list()
        return S([ i-1 for i in I ])

    def inverse_key(self, K):
        """
        The inverse key corresponding to the key ``K``. The specified key is a
        permutation of the cryptosystem alphabet.

        INPUT:

        - ``K`` - a key belonging to the key space of this cryptosystem

        OUTPUT:

        - The inverse key of ``K``.

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = SubstitutionCryptosystem(S)
            sage: K = E.random_key()
            sage: L = E.inverse_key(K)
            sage: M = S("THECATINTHEHAT")
            sage: e = E(K)
            sage: c = E(L)
            sage: c(e(M))
            THECATINTHEHAT
        """
        I = K._element_list
        S = self.cipher_domain()
        n = S.ngens()
        return S([ I.index(i) for i in range(n) ])

    def encoding(self, M):
        """
        The encoding of the string ``M`` over the string monoid of this
        substitution cipher. For example, if the string monoid of this
        cryptosystem is :class:`AlphabeticStringMonoid`, then the encoding
        of ``M`` would be its upper-case equivalent stripped of all
        non-alphabetic characters.

        INPUT:

        - ``M`` - a string, possibly empty

        OUTPUT:

        - The encoding of ``M`` over the string monoid of this cryptosystem.

        EXAMPLES::

            sage: M = "Peter Pan(ning) for gold."
            sage: A = AlphabeticStrings()
            sage: S = SubstitutionCryptosystem(A)
            sage: S.encoding(M) == A.encoding(M)
            True
        """
        S = self.cipher_domain()
        if isinstance(S,AlphabeticStringMonoid):
            return S(strip_encoding(M))
        try:
            return S.encoding(M)
        except:
            raise TypeError, "Argument M = %s does not encode in the cipher domain" % M

    def deciphering(self, K, C):
        """
        Decrypt the ciphertext ``C`` using the key ``K``.

        INPUT:

        - ``K`` - a key belonging to the key space of this substitution cipher

        - ``C`` - a string (possibly empty) over the string monoid of this
          cryptosystem.

        OUTPUT:

        - The plaintext corresponding to the ciphertext ``C``.

        EXAMPLES::

            sage: S = SubstitutionCryptosystem(AlphabeticStrings())
            sage: K = S.random_key()
            sage: M = S.encoding("Don't substitute me!")
            sage: S.deciphering(K, S.enciphering(K, M)) == M
            True
        """
        i = self(self.inverse_key(K))
        return i(C)

    def enciphering(self, K, M):
        """
        Encrypt the plaintext ``M`` using the key ``K``.

        INPUT:

        - ``K`` - a key belonging to the key space of this substitution cipher

        - ``M`` - a string (possibly empty) over the string monoid of this
          cryptosystem.

        OUTPUT:

        - The ciphertext corresponding to the plaintext ``M``.

        EXAMPLES::

            sage: S = SubstitutionCryptosystem(AlphabeticStrings())
            sage: K = S.random_key()
            sage: M = S.encoding("Don't substitute me.")
            sage: S.deciphering(K, S.enciphering(K, M)) == M
            True
        """
        e = self(K)
        return e(M)

class TranspositionCryptosystem(SymmetricKeyCryptosystem):
    """
    Transposition cryptosystem class
    """
    def __init__(self, S, n):
        """
        Create a transposition cryptosystem of block length ``n``.

        INPUT:

        - ``S`` - a string monoid over some alphabet

        - ``n`` - integer `> 0`; a block length of a block permutation

        OUTPUT:

        - A transposition cryptosystem of block length ``n`` over the
        alphabet ``S``.

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = TranspositionCryptosystem(S,14)
            sage: E
            Transposition cryptosystem on Free alphabetic string monoid on A-Z of block length 14
            sage: K = [ 14-i for i in range(14) ]
            sage: K
            [14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
            sage: e = E(K)
            sage: e(S("THECATINTHEHAT"))
            TAHEHTNITACEHT

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = TranspositionCryptosystem(S,14)
            sage: E == loads(dumps(E))
            True
        """
        if not isinstance(S, StringMonoid_class):
            raise TypeError, "S (= %s) must be a string monoid."%S
        key_space = SymmetricGroup(n)
        SymmetricKeyCryptosystem.__init__(self, S, S, key_space, block_length = n)

    def __call__(self, K):
        """
        Create a transposition cipher.

        INPUT:

        - ``K`` - a key which specifies a block permutation

        EXAMPLES::

            sage: M = AlphabeticStrings()
            sage: E = TranspositionCryptosystem(M,14)
            sage: E
            Transposition cryptosystem on Free alphabetic string monoid on A-Z of block length 14
            sage: K = [ 14-i for i in range(14) ]
            sage: K
            [14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
            sage: e = E(K)
            sage: m = M("THECATINTHEHAT")
            sage: e(m)
            TAHEHTNITACEHT
        """
        G = self.key_space()
        if isinstance(K, list):
            try:
                K = G(K)
            except:
                raise TypeError, "K (= %s) must specify a permutation."%K
        if not isinstance(K, PermutationGroupElement) and K.parent() == G:
            raise TypeError, "K (= %s) must be a permutation or list specifying a permutation."%K
        return TranspositionCipher(self, K)

    def _repr_(self):
        """
        Return a string representation of self.

        EXAMPLE::

            sage: A = AlphabeticStrings()
            sage: T = TranspositionCryptosystem(A, 14)
            sage: T
            Transposition cryptosystem on Free alphabetic string monoid on A-Z of block length 14
            sage: T._repr_()
            'Transposition cryptosystem on Free alphabetic string monoid on A-Z of block length 14'
        """
        return "Transposition cryptosystem on %s of block length %s" % (
            self.cipher_domain(), self.block_length())

    def random_key(self):
        """
        Generate a random key within the key space of this transposition
        cryptosystem. Let `n > 0` be the block length of this cryptosystem.
        Then there is a total of `n!` possible keys.

        OUTPUT:

        - A random key within the key space of this cryptosystem.

        EXAMPLE::

            sage: S = AlphabeticStrings()
            sage: E = TranspositionCryptosystem(S, 14)
            sage: K = E.random_key()
            sage: Ki = E.inverse_key(K)
            sage: e = E(K)
            sage: d = E(Ki)
            sage: M = "THECATINTHEHAT"
            sage: C = e(S(M))
            sage: d(S(C)) == S(M)
            True
        """
        n = self.block_length()
        return SymmetricGroup(n).random_element()

    def inverse_key(self, K, check=True):
        """
        The inverse key corresponding to the key ``K``.

        INPUT:

        - ``K`` - a key belonging to the key space of this transposition
          cipher

        - ``check`` - bool (default: ``True``); check that ``K`` belongs to
          the key space of this cryptosystem.

        OUTPUT:

        - The inverse key corresponding to ``K``.

        EXAMPLE::

            sage: S = AlphabeticStrings()
            sage: E = TranspositionCryptosystem(S, 14)
            sage: K = E.random_key()
            sage: Ki = E.inverse_key(K)
            sage: e = E(K)
            sage: d = E(Ki)
            sage: M = "THECATINTHEHAT"
            sage: C = e(S(M))
            sage: d(S(C)) == S(M)
            True
        """
        if check:
            if not K in self.key_space():
                raise TypeError, "Argument K (= %s) is not in the key space." % K
        return K**-1

    def encoding(self, M):
        """
        The encoding of the string ``M`` over the string monoid of this
        transposition cipher. For example, if the string monoid of this
        cryptosystem is :class:`AlphabeticStringMonoid`, then the encoding
        of ``M`` would be its upper-case equivalent stripped of all
        non-alphabetic characters.

        INPUT:

        - ``M`` - a string, possibly empty

        OUTPUT:

        - The encoding of ``M`` over the string monoid of this cryptosystem.

        EXAMPLES::

            sage: M = "Transposition cipher is not about matrix transpose."
            sage: A = AlphabeticStrings()
            sage: T = TranspositionCryptosystem(A, 11)
            sage: T.encoding(M) == A.encoding(M)
            True
        """
        S = self.cipher_domain()
        if isinstance(S,AlphabeticStringMonoid):
            return S(strip_encoding(M))
        try:
            return S.encoding(M)
        except:
            raise TypeError, "Argument M = %s does not encode in the cipher domain" % M

    def deciphering(self, K, C):
        """
        Decrypt the ciphertext ``C`` using the key ``K``.

        INPUT:

        - ``K`` - a key belonging to the key space of this transposition
          cipher

       - ``C`` - a string (possibly empty) over the string monoid of this
          cryptosystem.

        OUTPUT:

        - The plaintext corresponding to the ciphertext ``C``.

        EXAMPLES::

            sage: T = TranspositionCryptosystem(AlphabeticStrings(), 14)
            sage: K = T.random_key()
            sage: M = T.encoding("The cat in the hat.")
            sage: T.deciphering(K, T.enciphering(K, M)) == M
            True
        """
        i = self(self.inverse_key(K))
        return i(C)

    def enciphering(self, K, M):
        """
        Encrypt the plaintext ``M`` using the key ``K``.

        INPUT:

        - ``K`` - a key belonging to the key space of this transposition
          cipher

        - ``M`` - a string (possibly empty) over the string monoid of this
          cryptosystem

        OUTPUT:

        - The ciphertext corresponding to the plaintext ``M``.

        EXAMPLES::

            sage: T = TranspositionCryptosystem(AlphabeticStrings(), 14)
            sage: K = T.random_key()
            sage: M = T.encoding("The cat in the hat.")
            sage: T.deciphering(K, T.enciphering(K, M)) == M
            True
        """
        e = self(K)
        return e(M)

class VigenereCryptosystem(SymmetricKeyCryptosystem):
    """
    Vigenere cryptosystem class
    """

    def __init__(self, S, n):
        """
        Create a Vigenere cryptosystem of block length ``n``.

        INPUT:

        - ``S``-- a string monoid over some alphabet

        - ``n`` - integer `> 0`; block length of an encryption/decryption key

        OUTPUT:

        - A Vigenere cryptosystem of block length ``n`` over the alphabet
          ``S``.

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = VigenereCryptosystem(S,14)
            sage: E
            Vigenere cryptosystem on Free alphabetic string monoid on A-Z of period 14
            sage: K = S('ABCDEFGHIJKLMN')
            sage: K
            ABCDEFGHIJKLMN
            sage: e = E(K)
            sage: e
            ABCDEFGHIJKLMN
            sage: e(S("THECATINTHEHAT"))
            TIGFEYOUBQOSMG

        TESTS::

            sage: S = AlphabeticStrings()
            sage: E = VigenereCryptosystem(S,14)
            sage: E == loads(dumps(E))
            True
        """
        if not isinstance(S, StringMonoid_class):
            raise TypeError, "S (= %s) must be a string monoid."%S
        SymmetricKeyCryptosystem.__init__(self, S, S, S, block_length = 1, period = n)

    def __call__(self, K):
        """
        Create a Vigenere cipher.

        INPUT: A key which specifies a block permutation.

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = VigenereCryptosystem(S,14)
            sage: E
            Vigenere cryptosystem on Free alphabetic string monoid on A-Z of period 14
            sage: K = S('ABCDEFGHIJKLMN')
            sage: K
            ABCDEFGHIJKLMN
            sage: e = E(K)
            sage: e
            ABCDEFGHIJKLMN
            sage: e(S("THECATINTHEHAT"))
            TIGFEYOUBQOSMG
        """
        S = self.key_space()
        m = self.period()
        if isinstance(K, list):
            try:
                K = S(K)
            except:
                raise TypeError, "K (= %s) must specify a string of length %s." % (K, m)
        if not len(K) == m:
            raise TypeError, "K (= %s) must specify a string of length %s." % (K, m)
        return VigenereCipher(self, K)

    def _repr_(self):
        """
        Return a string representation of self.

        EXAMPLE::

            sage: A = AlphabeticStrings()
            sage: V = VigenereCryptosystem(A, 14)
            sage: V
            Vigenere cryptosystem on Free alphabetic string monoid on A-Z of period 14
            sage: V._repr_()
            'Vigenere cryptosystem on Free alphabetic string monoid on A-Z of period 14'
        """
        return "Vigenere cryptosystem on %s of period %s" % (
            self.cipher_domain(), self.period())

    def random_key(self):
        """
        Generate a random key within the key space of this Vigenere
        cryptosystem. Let `n > 0` be the length of the cryptosystem alphabet
        and let `m > 0` be the block length of this cryptosystem. Then there
        is a total of `n^m` possible keys.

        OUTPUT:

        - A random key within the key space of this cryptosystem.

        EXAMPLE::

            sage: A = AlphabeticStrings()
            sage: V = VigenereCryptosystem(A, 14)
            sage: M = "THECATINTHEHAT"
            sage: K = V.random_key()
            sage: Ki = V.inverse_key(K)
            sage: e = V(K)
            sage: d = V(Ki)
            sage: d(e(A(M))) == A(M)
            True
        """
        S = self.key_space()
        n = S.ngens()
        m = self.period()
        return S([ randint(0,n-1) for i in range(m) ])

    def inverse_key(self, K):
        """
        The inverse key corresponding to the key ``K``.

        INPUT:

        - ``K`` - a key within the key space of this Vigenere cryptosystem

        OUTPUT:

        - The inverse key corresponding to ``K``.

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = VigenereCryptosystem(S,14)
            sage: K = E.random_key()
            sage: L = E.inverse_key(K)
            sage: M = S("THECATINTHEHAT")
            sage: e = E(K)
            sage: c = E(L)
            sage: c(e(M))
            THECATINTHEHAT
        """
        S = self.key_space()
        n = S.ngens()
        return S([ (-i)%(n) for i in K._element_list ])

    def encoding(self, M):
        """
        The encoding of the string ``M`` over the string monoid of this
        Vigenere cipher. For example, if the string monoid of this
        cryptosystem is :class:`AlphabeticStringMonoid`, then the encoding
        of ``M`` would be its upper-case equivalent stripped of all
        non-alphabetic characters.

        INPUT:

        - ``M`` - a string, possibly empty

        OUTPUT:

        - The encoding of ``M`` over the string monoid of this cryptosystem.

        EXAMPLES::

            sage: A = AlphabeticStrings()
            sage: V = VigenereCryptosystem(A, 24)
            sage: M = "Jack and Jill went up the hill."
            sage: V.encoding(M) == A.encoding(M)
            True
        """
        S = self.cipher_domain()
        if isinstance(S,AlphabeticStringMonoid):
            return S(strip_encoding(M))
        try:
            return S.encoding(M)
        except:
            raise TypeError, "Argument M = %s does not encode in the cipher domain" % M

    def deciphering(self, K, C):
        """
        Decrypt the ciphertext ``C`` using the key ``K``.

        INPUT:

        - ``K`` - a key belonging to the key space of this Vigenere cipher

        - ``C`` - a string (possibly empty) over the string monoid of this
          cryptosystem

        OUTPUT:

        - The plaintext corresponding to the ciphertext ``C``.

        EXAMPLES::

            sage: V = VigenereCryptosystem(AlphabeticStrings(), 24)
            sage: K = V.random_key()
            sage: M = V.encoding("Jack and Jill went up the hill.")
            sage: V.deciphering(K, V.enciphering(K, M)) == M
            True
        """
        i = self(self.inverse_key(K))
        return i(C)

    def enciphering(self, K, M):
        """
        Encrypt the plaintext ``M`` using the key ``K``.

        INPUT:

        - ``K`` - a key belonging to the key space of this Vigenere cipher

        - ``M`` - a string (possibly empty) over the string monoid of this
          cryptosystem

        OUTPUT:

        - The ciphertext corresponding to the plaintext ``M``.

        EXAMPLES::

            sage: V = VigenereCryptosystem(AlphabeticStrings(), 24)
            sage: K = V.random_key()
            sage: M = V.encoding("Jack and Jill went up the hill.")
            sage: V.deciphering(K, V.enciphering(K, M)) == M
            True
        """
        e = self(K)
        return e(M)

