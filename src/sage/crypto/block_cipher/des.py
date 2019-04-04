r"""
DES

AUTHORS:

- Lukas Stennes (2019-03-29): initial version
"""

# ****************************************************************************
#       Copyright (C) 2013 Lukas Stennes <lukas.stennes@rub.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF


class DES(SageObject):
    r"""
    This class implements DES described in [TODO: ADD REF]_.
    """

    def __init__(self, rounds=None, keySchedule=False):
        r"""
        Construct an instance of DES.

        INPUT:

        - ``rounds``  -- integer (default: ``None``); the number of rounds. If
          ``None`` the number of rounds of the key schedule is used.

        - ``keySchedule`` -- (default: ``None``); the key schedule that will be
          used for encryption and decryption. If ``None`` the default DES key
          schedule is used.

        """
        if keySchedule is None:
            self._keySchedule = DES_KS()
        else:
            self._keySchedule = keySchedule
        if rounds is None:
            self._rounds = self._keySchedule._rounds
        elif rounds <= self._keySchedule._rounds:
            self._rounds = rounds
        else:
            raise ValueError('number of rounds must be less or equal to the '
                             'number of rounds of the key schedule')
        self._blocksize = 64

    def __call__(self, B, K, algorithm='encrypt'):
        r"""
        Apply DES encryption or decryption on ``B`` using the key ``K``.
        The flag ``algorithm`` controls what action is to be performed on
        ``B``.

        INPUT:

        - ``B`` -- integer or bit list-like; the plaintext or ciphertext

        - ``K`` -- integer or bit list-like; the key

        - ``algorithm`` -- string (default: ``'encrypt'``); a flag to signify
          whether encryption or decryption is to be applied to ``B``. The
          encryption flag is ``'encrypt'`` and the decryption flag is
          ``'decrypt'``

        OUTPUT:

        - The plaintext or ciphertext corresponding to ``B``, obtained using
          the key ``K``. If ``B`` is an integer the output will be too. If
          ``B`` is list-like the output will be a bit vector.
        """
        if algorithm == 'encrypt':
            return self.encrypt(B, K)
        elif algorithm == 'decrypt':
            return self.decrypt(B, K)
        else:
            raise ValueError('Algorithm must be \'encrypt\' or \'decrypt\' and'
                             ' not \'%s\'' % algorithm)

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        DES objects are the same if all attributes are the same.
        """
        if not isinstance(other, DES):
            return False
        else:
            return self.__dict__ == other.__dict__

    def __repr__(self):
        r"""
        A string representation of this DES.
        """
        raise NotImplementedError

    def encrypt(self, P, K):
        r"""
        Return the ciphertext corresponding to the plaintext ``P``,
        using DES encryption with key ``K``.
        """
        raise NotImplementedError

    def decrypt(self, C, K):
        r"""
        Return the plaintext corresponding to the ciphertext ``C``,
        using DES decryption with key ``K``.
        """
        raise NotImplementedError


class DES_KS(SageObject):
    r"""
    This class implements the DES key schedules described in [BKLPPRSV2007]_.
    """

    def __init__(self, rounds=16, master_key=None):
        r"""
        Construct an instance of DES_KS.
        """
        self._rounds = rounds
        self._master_key = master_key
        self._keysize = 56

    def __call__(self, K):
        r"""
        Return all round keys in a list.

        INPUT:

        - ``K`` -- integer or bit list-like; the key

        OUTPUT:

        - A list containing the round keys

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: K = vector(GF(2),[0,0,0,1,0,0,1,1,0,0,1,1,0,1,0,0,0,1,0,1,0,1,1,1,0,1,1,1,1,0,0,1,1,0,0,1,1,0,1,1,1,0,1,1,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,1])
            sage: ks = DES_KS(16, K)
            sage: [k for k in ks]
            [(0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0),
             (0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1),
             ...
             (1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1)]
            sage: K = vector(GF(2),[0,0,0,1,0,0,1,0,0,1,1,0,1,0,0,1,0,1,0,1,1,0,1,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,0,0,0])
            sage: ks = DES_KS()
            sage: ks(K)
            [(0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0),
             (0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1),
             ...
             (1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1)]
            sage: K = 0x1feded93da9648
            sage: ks = DES_KS(master_key=K)
            sage: [k.hex() for k in ks]
            ['4e0e3ff740d8',
             'a793db9b759e',
             ...
             'afe870d1bcd3']

        .. NOTE::

            If you want to use a DES_KS object as an iterable you have to
            pass a ``master_key`` value on initialisation. Otherwise you can
            omit ``master_key`` and pass a key when you call the object.
        """
        K, inputType = _convert_to_vector(K, 64)
        roundKeys = []
        # ensure that K is a 64 bit vector
        if not any(K[56:]):
            # delete msbs and insert 'parity' bits
            K = list(K)[:56]
            for i in range(7, 64, 8):
                K.insert(i, 0)
            K = vector(GF(2), 64, K)
        C, D = self._pc1(K)
        for i in range(16):
            C, D = self._left_shift(C, i), self._left_shift(D, i)
            roundKeys.append(self._pc2(list(C)+list(D)))
        return roundKeys if inputType == 'vector' else [ZZ(list(k), 2) for k in
                                                        roundKeys]

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        DES_KS objects are the same if all attributes are the same.
        """
        if not isinstance(other, DES_KS):
            return False
        else:
            return self.__dict__ == other.__dict__

    def __repr__(self):
        r"""
        A string representation of this DES_KS.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: DES_KS() # indirect doctest
            Original DES key schedule with 56-bit keys and 16 rounds
        """
        return ('Original DES key schedule with %s-bit keys and %s rounds'
                % (self._keysize, self._rounds))

    def __getitem__(self, r):
        r"""
        Computes the sub key for round ``r`` derived from initial master key.

        The key schedule object has to have been initialised with the
        `master_key` argument.

        INPUT:

        - ``r`` integer; the round for which the sub key is computed
        """
        if self._master_key is None:
            raise ValueError('Key not set during initialisation')
        return self(self._master_key)[r]

    def __iter__(self):
        r"""
        Iterate over the DES round keys, derived from
        `master_key`
       """
        if self._master_key is None:
            raise ValueError('Key not set during initialisation')
        return iter(self(self._master_key))

    def _pc1(self, K):
        r"""
        Compute Permuted Choice 1 of ``key``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: ks = DES_KS()
            sage: K = vector(GF(2),[0,0,0,1,0,0,1,1,0,0,1,1,0,1,0,0,0,1,0,1,0,1,1,1,0,1,1,1,1,0,0,1,1,0,0,1,1,0,1,1,1,0,1,1,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,1])
            sage: C, D = ks._pc1(K)
            sage: C
            (1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1)
            sage: D
            (0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1)
        """
        orderC = [57, 49, 41, 33, 25, 17,  9,
                   1, 58, 50, 42, 34, 26, 18,
                  10,  2, 59, 51, 43, 35, 27,
                  19, 11,  3, 60, 52, 44, 46]
        orderD = [63, 55, 47, 39, 31, 23, 15,
                   7, 62, 54, 46, 38, 30, 22,
                  14,  6, 61, 53, 45, 37, 29,
                  21, 13,  5, 28, 20, 12,  4]
        C = vector(GF(2), 28, [K[i-1] for i in orderC])
        D = vector(GF(2), 28, [K[i-1] for i in orderD])
        return C, D

    def _pc2(self, K):
        r"""
        Compute Permuted Choice 2 of ``key``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: ks = DES_KS()
            sage: K = vector(GF(2),[1,1,1,0,0,0,0,1,1,0,0,1,1,0,0,1,0,1,0,1,0,1,0,1,1,1,1,1,1,0,1,0,1,0,1,0,1,1,0,0,1,1,0,0,1,1,1,1,0,0,0,1,1,1,1,0])
            sage: ks._pc2(K)
            (0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0)
        """
        order = [14, 17, 11, 24,  1,  5,
                  3, 28, 15,  6, 21, 10,
                 23, 19, 12,  4, 26,  8,
                 16,  7, 27, 20, 13,  2,
                 41, 52, 31, 37, 47, 55,
                 30, 40, 51, 45, 33, 48,
                 44, 49, 39, 56, 34, 53,
                 46, 42, 50, 36, 29, 32]
        return vector(GF(2), 48, [K[i-1] for i in order])

    def _left_shift(self, half, i):
        r"""
        Shift ``half`` one or two positions to the left depending on the
        iteration number ``i``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: ks = DES_KS()
            sage: bits = vector(GF(2), 6, [1,0,1,0,1,0])
            sage: ks._left_shift(bits, 0)
            (0, 1, 0, 1, 0, 1)
            sage: bits
            (1, 0, 1, 0, 1, 0)
            sage: ks._left_shift(bits, 2)
            (1, 0, 1, 0, 1, 0)
        """
        amount = [1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1]
        return vector(GF(2),
                      list(half[amount[i]:]) + list(half[0:amount[i]]))


def _convert_to_vector(I, L):
    r"""
    Convert ``I`` to a bit vector of length ``L``.

    INPUT:

    - ``I`` -- integer or bit list-like

    - ``L`` -- integer; the desired bit length of the ouput

    OUTPUT: a tuple of

    - the ``L``-bit vector representation of ``I``

    - a flag indicating the input type. Either ``'integer'`` or ``'vector'``.

    EXAMPLES::

        sage: from sage.crypto.block_cipher.des import _convert_to_vector
        sage: _convert_to_vector(0x1F, 8)
        ((1, 1, 1, 1, 1, 0, 0, 0), 'integer')
        sage: v = vector(GF(2), 4, [1,0,1,0])
        sage: _convert_to_vector(v, 4)
        ((1, 0, 1, 0), 'vector')
    """
    try:
        state = vector(GF(2), L, ZZ(I).digits(2, padto=L))
        return state, 'integer'
    except TypeError:
        # ignore the error and try list-like types
        pass
    state = vector(GF(2), L, ZZ(list(I), 2).digits(2, padto=L))
    return state, 'vector'
