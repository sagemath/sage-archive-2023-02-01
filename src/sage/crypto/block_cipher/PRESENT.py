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
from itertools import chain


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

        sage: from sage.crypto.block_cipher.PRESENT import PRESENT
        sage: present = PRESENT(keysize=80)
        sage: p1 = "0000000000000000"
        sage: k1 = "00000000000000000000"
        sage: c1 = "5579C1387B228445"
        sage: present.encrypt(p1, k1) == c1
        True
        sage: present.decrypt(c1, k1) == p1
        True
        sage: p2 = "0000000000000000"
        sage: k2 = "FFFFFFFFFFFFFFFFFFFF"
        sage: c2 = "E72C46C0F5945049"
        sage: present.encrypt(p2, k2) == c2
        True
        sage: present.decrypt(c2, k2) == p2
        True
        sage: p3 = "FFFFFFFFFFFFFFFF"
        sage: k3 = "00000000000000000000"
        sage: c3 = "A112FFC72F68417B"
        sage: present.encrypt(p3, k3) == c3
        True
        sage: present.decrypt(c3, k3) == p3
        True
        sage: p4 = "FFFFFFFFFFFFFFFF"
        sage: k4 = "FFFFFFFFFFFFFFFFFFFF"
        sage: c4 = "3333DCD3213210D2"
        sage: present.encrypt(p4, k4) == c4
        True
        sage: present.decrypt(c4, k4) == p4
        True
    """

    def __init__(self, keysize=80):
        r"""
        Construct an instance of PRESENT.

        INPUT:

        - ``keysize`` -- (default: ``80``) the size of the keys that will be
          used in bits. It must be either 80 or 128.
        """
        if keysize != 80 and keysize != 128:
            raise ValueError("keysize must bei either 80 or 128 and not %s"
                             % keysize)
        self.keysize = keysize
        self.blocksize = 64
        from sage.crypto.sboxes import PRESENT as PRESENTSBOX
        self.sbox = PRESENTSBOX
        self.inverseSbox = self.sbox.inverse()
        self.permutationMatrix = smallscale_present_linearlayer()
        self.inversePermutationMatrix = self.permutationMatrix.inverse()

    def __call__(self, B, K, algorithm="encrypt"):
        r"""
        Apply PRESENT encryption or decryption on ``B`` using the key ``K``.
        The flag ``algorithm`` controls what action is to be performed on
        ``B``.

        INPUT:

        - ``B`` -- The plaintext or the ciphertext

        - ``K`` -- a string of 16 or 32 hex digits; The key that will be used.
          The rightmost hex digit represents `k_3k_2k_1k_0`.

        - ``algorithm`` -- (default: ``"encrypt"``) a string; a flag to signify
          whether encryption or decryption is to be applied to ``B``. The
          encryption flag is ``"encrypt"`` and the decryption flag is
          ``"decrypt"``.

        OUTPUT:

        - The plaintext or ciphertext corresponding to ``B``, obtained using
          the key ``K``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.PRESENT import PRESENT
            sage: present = PRESENT(keysize=80)
            sage: P = "FFFFFFFFFFFFFFFF"
            sage: K = "00000000000000000000"
            sage: present(present(P, K, "encrypt"), K, "decrypt") == P
            True
        """
        if algorithm == "encrypt":
            return self.encrypt(B, K)
        elif algorithm == "decrypt":
            return self.decrypt(B, K)
        else:
            raise ValueError("Algorithm mus be \"encrypt\" or \"decrypt\" and"
                             "not \"%s\"" % algorithm)

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        PRESENT objects are the same if the keysize is the same.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.PRESENT import PRESENT
            sage: PRESENT(80) == PRESENT(80)
            True
            sage: PRESENT(80) == PRESENT(128)
            False
        """
        return self.keysize == other.keysize

    def __repr__(self):
        r"""
        A string representation of this PRESENT.
        """
        return "PRESENT block cipher with %s-bit keys" % self.keysize

    def generateRoundKeys(self, K):
        r"""
        Return all round keys `K_1 \dots K_{32}` in a list.

        INPUT:

        - ``K`` --  a string of 16 or 32 hex digits; The rightmost hex digit
          represents `k_3k_2k_1k_0`.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.PRESENT import PRESENT
            sage: present = PRESENT(80)
            sage: K = present.generateRoundKeys("00000000000000000000")
            sage: ZZ(K[0], 2) == ZZ("0000000000000000", 16)
            True
            sage: ZZ(K[1], 2) == ZZ("c000000000000000", 16)
            True
            sage: ZZ(K[2], 2) == ZZ("5000180000000001", 16)
            True
            sage: ZZ(K[31], 2) == ZZ("6dab31744f41d700", 16)
            True


        Description of the key schedule for 64-bit keys from [BKLPPRSV2007]_:

        At round `i` the 64-bit round key `K_i = \kappa_{63}\kappa_{62} \dots
        \kappa_0` consists of the 64 leftmost bits fo the current contents of
        register `K`. After extracting the round key `K_i`, the key register
        `K = k_{79}k_{78} \dots k_0` is updated as follows.

        .. MATH::
           :nowrap:

            \begin{aligned}
             \ [k_{79}k_{78} \dots k_{1}k_{0}] &=
             [k_{18}k_{17} \dots k_{20}k_{19}] \\
             [k_{79}k_{78}k_{77}k_{76}] &= S[k_{79}k_{78}k_{77}k_{76}] \\
             [k_{19}k_{18}k_{17}k_{16}k_{15}] &=
             [k_{19}k_{18}k_{17}k_{16}k_{15}] \oplus round\_counter
            \end{aligned}

        Thus, the key register is rotated by 61 bit positions to the left, the
        left-most four bits are passed through the PRESENT S-box, and the
        round_counter value `i` is exclusive-ored with bits `k_{19}k_{18}k_{17}
        k_{16}k_{15}` of `K` with the least significant bit of round_counter on
        the right.

        The key schedule for 128-bit keys is not implemented yet.
        """
        if len(K) == 20 and self.keysize == 80:
            # convert K to list of key bits [k_0,...,k_79]
            roundKeys = []
            K = (ZZ(K, 16).bits() + [0] * (80 - ZZ(K, 16).nbits()))
            for i in range(1, 33):
                roundKeys.append(K[16:])
                K = K[19:] + K[:19]
                K[76:] = self.sbox((K[76:])[::-1])[::-1]
                rc = (ZZ(i).bits() + [0] * (5 - ZZ(i).nbits()))
                K[15:20] = [int(K[15+j]) ^ int(rc[j]) for j in range(5)]
            return roundKeys
        elif len(K) == 32 and self.keysize == 128:
            raise NotImplementedError("The key schedule for 128-bit keys is"
                                      " not implemented yet.")
        else:
            raise ValueError("Key must be %s-bits long" % self.keysize)

    def decrypt(self, C, K):
        r"""
        Return an plaintext corresponding to the ciphertext ``C``,
        using PRESENT decryption with key ``K``.

        INPUT:

        - ``C`` -- The ciphertext that will be encrypted.

        - ``K`` -- a string of 16 or 32 hex digits; The key that will be used
          to encrypt ``P``. The rightmost hex digit represents `k_3k_2k_1k_0`.

        OUTPUT:

        - The plaintext corresponding to ``C``, obtained using the key ``K``.

        """
        state = (ZZ(C, 16).bits() + [0] * (64 - ZZ(C, 16).nbits()))
        K = self.generateRoundKeys(K)
        state = [int(state[j]) ^ int(K[31][j]) for j in range(64)]
        for i in range(30, -1, -1):
            state = self.inversePermutationMatrix * vector(GF(2), state)
            state = list(chain.from_iterable([self.inverseSbox(
                state[4*j:4*j+4][::-1])[::-1] for j in range(16)]))
            state = [int(state[j]) ^ int(K[i][j]) for j in range(64)]
        return ZZ(state, 2).hex().upper()

    def encrypt(self, P, K):
        r"""
        Return an ciphertext corresponding to the plaintext ``P``,
        using PRESENT encryption with key ``K``.

        INPUT:

        - ``P`` -- The plaintext that will be encrypted.

        - ``K`` -- a string of 16 or 32 hex digits; The key that will be used
          to encrypt ``P``. The rightmost hex digit represents `k_3k_2k_1k_0`.

        OUTPUT:

        - The ciphertext corresponding to ``P``, obtained using the key ``K``.
        """
        state = (ZZ(P, 16).bits() + [0] * (64 - ZZ(P, 16).nbits()))
        K = self.generateRoundKeys(K)
        for i in range(0, 31):
            state = [int(state[j]) ^ int(K[i][j]) for j in range(64)]
            state = list(chain.from_iterable(
                [self.sbox(state[4*j:4*j+4][::-1])[::-1] for j in range(16)]))
            state = self.permutationMatrix * vector(GF(2), state)
        state = [int(state[j]) ^ int(K[31][j]) for j in range(64)]
        return ZZ(state, 2).hex().upper()
