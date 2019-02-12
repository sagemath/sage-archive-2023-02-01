r"""
PRESENT

An ultra-lightweight block cipher.

AUTHORS:

- Lukas Stennes (2019-02-01): initial version
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
from sage.crypto.sboxes import PRESENT as PRESENTSBOX


def smallscale_present_linearlayer(nsboxes=16):
    """
    TODO: switch to sage.crypto.linearlayer as soon as it is included in sage
    """
    from sage.modules.free_module import VectorSpace
    from sage.modules.free_module_element import vector
    from sage.matrix.constructor import Matrix
    from sage.rings.finite_rings.finite_field_constructor import GF

    def present_llayer(n, x):
        dim = 4*n
        y = [0]*dim
        for i in range(dim-1):
            y[i] = x[(n * i) % (dim - 1)]
        y[dim-1] = x[dim-1]
        return vector(GF(2), y)

    m = Matrix(GF(2), [present_llayer(nsboxes, ei)
                       for ei in VectorSpace(GF(2), 4*nsboxes).basis()])
    return m


class PRESENT(SageObject):
    r"""
    This class implements PRESENT described in [BKLPPRSV2007]_.

    TESTS:

    Check test vectors given in [BKLPPRSV2007]_.

        sage: from sage.crypto.block_cipher.present import PRESENT
        sage: present = PRESENT(keySchedule=80)
        sage: p1 = 0x0
        sage: k1 = 0x0
        sage: c1 = 0x5579C1387B228445
        sage: present.encrypt(p1, k1) == c1
        True
        sage: present.decrypt(c1, k1) == p1
        True
        sage: p2 = 0x0
        sage: k2 = 0xFFFFFFFFFFFFFFFFFFFF
        sage: c2 = 0xE72C46C0F5945049
        sage: present.encrypt(p2, k2) == c2
        True
        sage: present.decrypt(c2, k2) == p2
        True
        sage: p3 = 0xFFFFFFFFFFFFFFFF
        sage: k3 = 0x0
        sage: c3 = 0xA112FFC72F68417B
        sage: present.encrypt(p3, k3) == c3
        True
        sage: present.decrypt(c3, k3) == p3
        True
        sage: p4 = 0xFFFFFFFFFFFFFFFF
        sage: k4 = 0xFFFFFFFFFFFFFFFFFFFF
        sage: c4 = 0x3333DCD3213210D2
        sage: present.encrypt(p4, k4) == c4
        True
        sage: present.decrypt(c4, k4) == p4
        True
    """

    def __init__(self, keySchedule=80):
        r"""
        Construct an instance of PRESENT.

        INPUT:

        - ``keySchedule`` --
        """
        if keySchedule == 80:
            self._keySchedule = PRESENT_KS()
        elif keySchedule == 128:
            self._keySchedule = PRESENT_KS(128)
        else:
            self._keySchedule = keySchedule
        self._blocksize = 64
        self._sbox = PRESENTSBOX
        self._inverseSbox = self._sbox.inverse()
        self._permutationMatrix = smallscale_present_linearlayer()
        self._inversePermutationMatrix = self._permutationMatrix.inverse()

    def __call__(self, B, K, algorithm="encrypt"):
        r"""
        Apply PRESENT encryption or decryption on ``B`` using the key ``K``.
        The flag ``algorithm`` controls what action is to be performed on
        ``B``.

        INPUT:

        - ``B`` -- integer or list-like; the plaintext or ciphertext

        - ``K`` -- integer or list-like; the key

        - ``algorithm`` -- (default: ``"encrypt"``) a string; a flag to signify
          whether encryption or decryption is to be applied to ``B``. The
          encryption flag is ``"encrypt"`` and the decryption flag is
          ``"decrypt"``

        OUTPUT:

        - The plaintext or ciphertext corresponding to ``B``, obtained using
          the key ``K``

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: present = PRESENT(keySchedule=80)
            sage: P = 0xFFFFFFFFFFFFFFFF
            sage: K = 0x0
            sage: present(present(P, K, "encrypt"), K, "decrypt") == P
            True
            sage: P = ZZ(0).digits(2,padto=64)
            sage: K = ZZ(0).digits(2,padto=80)
            sage: present(present(P, K, "encrypt"), K, "decrypt") == P
            True
            sage: P = ZZ(0).digits(16,padto=16)
            sage: K = ZZ(0).digits(16,padto=20)
            sage: present(present(P, K, "encrypt"), K, "decrypt") == P
            True
        """
        if algorithm == "encrypt":
            return self.encrypt(B, K)
        elif algorithm == "decrypt":
            return self.decrypt(B, K)
        else:
            raise ValueError("Algorithm must be \"encrypt\" or \"decrypt\" and"
                             "not \"%s\"" % algorithm)

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        PRESENT objects are the same if the keysize is the same.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: PRESENT(80) == PRESENT(80)
            True
            sage: PRESENT(80) == PRESENT(128)
            False
        """
        return self._keysize == other._keysize

    def __repr__(self):
        r"""
        A string representation of this PRESENT.
        """
        return "PRESENT block cipher with %s-bit keys" % self._keysize

    def _to_state(self, I, L):
        r"""
        Convert ``I`` to a bit vector of length ``L``.

        If ``I`` is list-like its elements will be treated as binary if
        ``len(I) == L`` or as hexadecimal if ``4*len(I) == L``.

        INPUT:

        - ``I`` -- integer or list-like

        - ``L`` -- integer; the desired length of the ouput

        OUTPUT:

        - A ``L``-bit vector representation of ``I``

        - A flag describing the type of ``I``. Either "int" or "bin" or "hex"

        EXAMPLES::

        sage: from sage.crypto.block_cipher.present import PRESENT
        sage: present = PRESENT()
        sage: present._to_state("0x1F", 8)
        ((1, 1, 1, 1, 1, 0, 0, 0), 'int')
        sage: present._to_state(["0x1","0xF"], 8)
        ((1, 0, 0, 0, 1, 1, 1, 1), 'hex')
        sage: present._to_state(0x1F, 8)
        ((1, 1, 1, 1, 1, 0, 0, 0), 'int')
        sage: v = vector(GF(2), 4, [1,0,1,0])
        sage: present._to_state(v, 4)
        ((1, 0, 1, 0), 'bin')
        sage: present._to_state(0x1F, 9)
        ((1, 1, 1, 1, 1, 0, 0, 0, 0), 'int')
        sage: present._to_state(["1","0xF"], 9)
        Traceback (most recent call last):
        ...
        ValueError: ['1', '0xF'] can not be converted to bit vector of length 9
        """
        try:
            state = vector(GF(2), L, ZZ(I).digits(2, padto=L))
            return state, "int"
        except TypeError:
            pass
        if len(I) == L:
            state = vector(GF(2), L, ZZ(list(I), 2).digits(2, padto=L))
            flag = "bin"
        elif 4*len(I) == L:
            state = vector(GF(2), L, ZZ(list(I), 16).digits(2, padto=L))
            flag = "hex"
        else:
            raise ValueError("%s can not be converted to bit vector of "
                             "length %s" % (I, L))
        return state, flag

    def decrypt(self, C, K):
        r"""
        Return an plaintext corresponding to the ciphertext ``C``,
        using PRESENT decryption with key ``K``.

        INPUT:

        - ``C`` -- integer or list-like; the plaintext that will be decrypted.

        - ``K`` -- integer or list-like; the key

        OUTPUT:

        - The plaintext corresponding to ``C``, obtained using the key ``K``.
        """
        state, inputType = self._to_state(C, 64)
        roundKeys = self.generate_round_keys(K)
        state = state + roundKeys[-1]
        for K in roundKeys[:-1][::-1]:
            state[0:] = self._inversePermutationMatrix * state
            for nibble in [slice(4*j, 4*j+4) for j in range(16)]:
                state[nibble] = self._inverseSbox(state[nibble][::-1])[::-1]
            state = state + K
        if inputType == "bin":
            P = list(state)
        elif inputType == "hex":
            P = ZZ(list(state), 2).digits(16, padto=16)
        elif inputType == "int":
            P = ZZ(list(state), 2)
        return P

    def encrypt(self, P, K):
        r"""
        Return an ciphertext corresponding to the plaintext ``P``,
        using PRESENT encryption with key ``K``.

        INPUT:

        - ``P`` -- integer or list-like; the plaintext that will be encrypted.

        - ``K`` -- integer or list-like; the key

        OUTPUT:

        - The ciphertext corresponding to ``P``, obtained using the key ``K``.
        """
        state, inputType = self._to_state(P, 64)
        roundKeys = self.generate_round_keys(K)
        for K in roundKeys[:-1]:
            state = state + K
            for nibble in [slice(4*j, 4*j+4) for j in range(16)]:
                state[nibble] = self._sbox(state[nibble][::-1])[::-1]
            state[0:] = self._permutationMatrix * state
        state = state + roundKeys[-1]
        if inputType == "bin":
            C = list(state)
        elif inputType == "hex":
            C = ZZ(list(state), 2).digits(16, padto=16)
        elif inputType == "int":
            C = ZZ(list(state), 2)
        return C


class PRESENT_KS(SageObject):
    r"""
    """
    def __init__(self, keysize=80, rounds=31):
        self._keysize = keysize
        self._rounds = rounds
        self._sbox = PRESENTSBOX

    def __call__(self, K):
        if self._keysize == 80:
            roundKeys = []
            for i in range(1, self._rounds+1):
                roundKeys.append(K[16:])
                K[0:] = list(K[19:]) + list(K[:19])
                K[76:] = self._sbox(K[76:][::-1])[::-1]
                rc = vector(GF(2), ZZ(i).digits(2, padto=5))
                K[15:20] = K[15:20] + rc
            roundKeys.append(K[16:])
            return roundKeys
        elif self._keysize == 128:
            roundKeys = []
            for i in range(1, self._rounds+1):
                roundKeys.append(K[64:])
                K[0:] = list(K[67:]) + list(K[:67])
                K[124:] = self._sbox(K[124:][::-1])[::-1]
                K[120:124] = self._sbox(K[120:124][::-1])[::-1]
                rc = vector(GF(2), ZZ(i).digits(2, padto=5))
                K[62:67] = K[62:67] + rc
            roundKeys.append(K[64:])
            return roundKeys

    def __eq__(self, other):
        return self.__dict__ == other.__dict__
