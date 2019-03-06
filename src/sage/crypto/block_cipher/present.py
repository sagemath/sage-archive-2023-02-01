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

    EXAMPLES::

        sage: from sage.crypto.block_cipher.present import PRESENT
        sage: present = PRESENT()
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
        sage: present = PRESENT(128)
        sage: P = 0x0123456789abcdef
        sage: K = 0x00112233445566778899aabbccddeeff
        sage: present(present(P, K, "encrypt"), K, "decrypt") == P
        True
    """

    def __init__(self, keySchedule=80):
        r"""
        Construct an instance of PRESENT.

        INPUT:

        - ``keySchedule`` -- the key schedule that will be used for encryption
            and decryption. Use ``80`` or ``128`` to use the original key
            schedules from [BKLPPRSV2007]_

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: from sage.crypto.block_cipher.present import PRESENT_KS
            sage: PRESENT() # indirect doctest
            PRESENT block cipher with the following key schedule:
            Original PRESENT key schedule with 80-bit keys and 31 rounds
            sage: PRESENT(128) # indirect doctest
            PRESENT block cipher with the following key schedule:
            Original PRESENT key schedule with 128-bit keys and 31 rounds
            sage: PRESENT(PRESENT_KS(80, 15)) # indirect doctest
            PRESENT block cipher with the following key schedule:
            Original PRESENT key schedule with 80-bit keys and 15 rounds

        .. SEEALSO::

            :class: `PRESENT_KS`
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

        PRESENT objects are the same if the keySchedule is the same.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: PRESENT(80) == PRESENT(80) # indirect doctest
            True
            sage: PRESENT(80) == PRESENT(128) # indirect doctest
            False
            sage: PRESENT(80) == 80 # indirect doctest
            False
            sage: present = PRESENT()
            sage: present._inverseSbox = present._sbox
            sage: present == PRESENT() # indirect doctest
            False
        """
        try:
            return self.__dict__ == other.__dict__
        except AttributeError:
            # if other has not attribute __dict__
            return False

    def __repr__(self):
        r"""
        A string representation of this PRESENT.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: PRESENT() # indirect doctest
            PRESENT block cipher with the following key schedule:
            Original PRESENT key schedule with 80-bit keys and 31 rounds
        """
        return ("PRESENT block cipher with the following key schedule:\n%s"
                % self._keySchedule.__repr__())

    def decrypt(self, C, K):
        r"""
        Return an plaintext corresponding to the ciphertext ``C``,
        using PRESENT decryption with key ``K``.

        INPUT:

        - ``C`` -- integer or list-like; the plaintext that will be decrypted.

        - ``K`` -- integer or list-like; the key

        OUTPUT:

        - The plaintext corresponding to ``C``, obtained using the key ``K``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: present = PRESENT()
            sage: p1 = 0x0
            sage: k1 = 0x0
            sage: c1 = 0x5579C1387B228445
            sage: present.decrypt(c1, k1) == p1
            True
            sage: p2 = 0x0
            sage: k2 = 0xFFFFFFFFFFFFFFFFFFFF
            sage: c2 = 0xE72C46C0F5945049
            sage: present.decrypt(c2, k2) == p2
            True
            sage: p3 = 0xFFFFFFFFFFFFFFFF
            sage: k3 = 0x0
            sage: c3 = 0xA112FFC72F68417B
            sage: present.decrypt(c3, k3) == p3
            True
            sage: p4 = 0xFFFFFFFFFFFFFFFF
            sage: k4 = 0xFFFFFFFFFFFFFFFFFFFF
            sage: c4 = 0x3333DCD3213210D2
            sage: present.decrypt(c4, k4) == p4
            True
        """
        state, inputType = convert_to_vector(C, 64)
        K, _ = convert_to_vector(K, self._keySchedule._keysize)
        roundKeys = self._keySchedule(K)
        state = state + roundKeys[-1]
        for K in roundKeys[:-1][::-1]:
            state[0:] = self._inversePermutationMatrix * state
            for nibble in [slice(4*j, 4*j+4) for j in range(16)]:
                state[nibble] = self._inverseSbox(state[nibble][::-1])[::-1]
            state = state + K
        P = convert_vector_to_type(state, inputType)
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

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: present = PRESENT()
            sage: present = PRESENT(keySchedule=80)
            sage: p1 = 0x0
            sage: k1 = 0x0
            sage: c1 = 0x5579C1387B228445
            sage: present.encrypt(p1, k1) == c1
            True
            sage: p2 = 0x0
            sage: k2 = 0xFFFFFFFFFFFFFFFFFFFF
            sage: c2 = 0xE72C46C0F5945049
            sage: present.encrypt(p2, k2) == c2
            True
            sage: p3 = 0xFFFFFFFFFFFFFFFF
            sage: k3 = 0x0
            sage: c3 = 0xA112FFC72F68417B
            sage: present.encrypt(p3, k3) == c3
            True
            sage: p4 = 0xFFFFFFFFFFFFFFFF
            sage: k4 = 0xFFFFFFFFFFFFFFFFFFFF
            sage: c4 = 0x3333DCD3213210D2
            sage: present.encrypt(p4, k4) == c4
            True
        """
        state, inputType = convert_to_vector(P, 64)
        K, _ = convert_to_vector(K, self._keySchedule._keysize)
        roundKeys = self._keySchedule(K)
        for K in roundKeys[:-1]:
            state = state + K
            for nibble in [slice(4*j, 4*j+4) for j in range(16)]:
                state[nibble] = self._sbox(state[nibble][::-1])[::-1]
            state[0:] = self._permutationMatrix * state
        state = state + roundKeys[-1]
        C = convert_vector_to_type(state, inputType)
        return C


class PRESENT_KS(SageObject):
    r"""
    This class implements the PRESENT key schedules for both 80-bit and 128-bit
    keys as described in [BKLPPRSV2007]_.

    EXAMPLES::

        sage: from sage.crypto.block_cipher.present import PRESENT_KS
        sage: ks = PRESENT_KS()
        sage: K = ks(0x0)
        sage: K[0] == 0x0
        True
        sage: K[31] == 0x6dab31744f41d700
        True
        sage: ks = PRESENT_KS(128)
        sage: K = ks(0x00112233445566778899aabbccddeeff)
        sage: K[0] == 0x0011223344556677
        True
        sage: K[31] == 0x091989a5ae8eab21
        True

    Description of the key schedule for 64-bit and 128-bit keys from
    [BKLPPRSV2007]_:

    The key schedule for 64-bit keys works as follows:

    At round `i` the 64-bit round key `K_i = \kappa_{63}\kappa_{62} \dots
    \kappa_0` consists of the 64 leftmost bits fo the current contents of
    register `K`. Thus at round `i` we have that:

    .. MATH::

        K_i = \kappa_{63}\kappa_{62}\dots \kappa_0 = k_{79}k_{78}\dots k_{16}

    After extracting the round key `K_i`, the key register `K = k_{79}k_{78}
    \dots k_0` is updated as follows.

    .. MATH::
        :nowrap:

        \begin{aligned}
        \ [k_{79}k_{78}\dots k_{1}k_{0}] &= [k_{18}k_{17}\dots k_{20}k_{19}] \\
        [k_{79}k_{78}k_{77}k_{76}] &= S[k_{79}k_{78}k_{77}k_{76}] \\
        [k_{19}k_{18}k_{17}k_{16}k_{15}] &= [k_{19}k_{18}k_{17}k_{16}k_{15}]
                                            \oplus round\_counter
        \end{aligned}

    Thus, the key register is rotated by 61 bit positions to the left, the
    left-most four bits are passed through the PRESENT S-box, and the
    round_counter value `i` is exclusive-ored with bits `k_{19}k_{18}k_{17}
    k_{16}k_{15}` of `K` with the least significant bit of round_counter on
    the right.

    The key schedule for 128-bit keys works as follows:

    At round `i` the 64-bit round key `K_i = \kappa_{63}\kappa_{62} \dots
    \kappa_0` consists of the 64 leftmost bits fo the current contents of
    register `K`. Thus at round `i` we have that:

    .. MATH::

        K_i = \kappa_{63}\kappa_{62}\dots \kappa_0 = k_{127}k_{126}\dots k_{64}

    After extracting the round key `K_i`, the key register `K = k_{127}k_{126}
    \dots k_0` is updated as follows:

    .. MATH::
        :nowrap:

        \begin{aligned}
        \ [k_{127}k_{126}\dots k_{1}k_{0}] &= [k_{66}k_{65}\dots k_{68}k_{67}]
        \\
        [k_{127}k_{126}k_{125}k_{124}] &= S[k_{127}k_{126}k_{125}k_{124}]\\
        [k_{123}k_{122}k_{121}k_{120}] &= S[k_{123}k_{122}k_{121}k_{120}]\\
        [k_{66}k_{65}k_{64}k_{63}k_{62}] &= [k_{66}k_{65}k_{64}k_{63}k_{62}]
                                            \oplus round\_counter
        \end{aligned}

    Thus, the key register is rotated by 61 bit positions to the left, the
    left-most eight bits are passed through two PRESENT S-boxes, and the
    round_counter value `i` is exclusive-ored with bits `k_{66}k_{65}k_{64}
    k_{63}k_{62}` of `K` with the least significant bit of round_counter on
    the right.
    """

    def __init__(self, keysize=80, rounds=31, master_key=None):
        r"""
        Construct an instance of PRESENT_KS.

        INPUT:

         - ``keysize`` -- integer; the size of the keys that will be
           used in bits. It must be either 80 or 128 (default: ``80``).

        - ``rounds`` -- integer; the number of rounds ``self`` can create keys
          for.

        - ``master_key`` -- integer of list-like; the key that will be used
        """
        if keysize != 80 and keysize != 128:
            raise ValueError("keysize must bei either 80 or 128 and not %s"
                             % keysize)
        self.keysize = keysize
        self._keysize = keysize
        self._rounds = rounds
        self._sbox = PRESENTSBOX
        self._master_key = master_key

    def __call__(self, K):
        r"""
        Return all round keys in a list.

        INPUT:

        - ``K`` -- integer or list-like; the key

        OUTPUT:

        - A list containing ``rounds + 1`` round keys.
        """
        K, inputType = convert_to_vector(K, self._keysize)
        roundKeys = []
        if self._keysize == 80:
            for i in range(1, self._rounds+1):
                roundKeys.append(convert_vector_to_type(K[16:], inputType))
                K[0:] = list(K[19:]) + list(K[:19])
                K[76:] = self._sbox(K[76:][::-1])[::-1]
                rc = vector(GF(2), ZZ(i).digits(2, padto=5))
                K[15:20] = K[15:20] + rc
            roundKeys.append(convert_vector_to_type(K[16:], inputType))
        elif self._keysize == 128:
            for i in range(1, self._rounds+1):
                roundKeys.append(convert_vector_to_type(K[64:], inputType))
                K[0:] = list(K[67:]) + list(K[:67])
                K[124:] = self._sbox(K[124:][::-1])[::-1]
                K[120:124] = self._sbox(K[120:124][::-1])[::-1]
                rc = vector(GF(2), ZZ(i).digits(2, padto=5))
                K[62:67] = K[62:67] + rc
            roundKeys.append(convert_vector_to_type(K[64:], inputType))
        return roundKeys

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        PRESENT_KS objects are the same if all attributes are the same.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT_KS
            sage: PRESENT_KS(80) == PRESENT_KS(80) # indirect doctest
            True
            sage: PRESENT_KS(80) == PRESENT_KS(128) # indirect doctest
            False
            sage: PRESENT_KS(80) == 80 # indirect doctest
            False
        """
        try:
            return self.__dict__ == other.__dict__
        except AttributeError:
            # if other has no attribute __dict__
            return False

    def __repr__(self):
        r"""
        A string representation of this PRESENT_KS.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT_KS
            sage: PRESENT_KS() # indirect doctest
            Original PRESENT key schedule with 80-bit keys and 31 rounds
        """
        return ("Original PRESENT key schedule with %s-bit keys and %s rounds"
                % (self._keysize, self._rounds))

    def __getitem__(self, r):
        r"""
        Computes the sub key for round ``r`` derived from initial master key.

        The key schedule object has to have been initialised with the
        `master_key` argument.

        INPUT:

        - ``r`` integer; the round for which the sub key is computed

        EXAMPLES:

        Check against test vectors.::

            sage: from sage.crypto.block_cipher.present import PRESENT_KS
            sage: ks = PRESENT_KS(master_key=0x0)
            sage: ks[0] ==  0x0 # indirect doctest
            True
            sage: ks[31] ==  0x6dab31744f41d700 # indirect doctest
            True
        """
        if self._master_key is None:
            raise ValueError("Key not set during initialisation")
        return self(self._master_key)[r]

    def __iter__(self):
        """
        iterate over the ``self._rounds + 1`` PRESENT round keys, derived from
        `master_key`

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT_KS
            sage: K = [k for k in PRESENT_KS(master_key=0x0)]
            sage: K[0] == 0x0 # indirect doctest
            True
            sage: K[31] == 0x6dab31744f41d700 # indirect doctest
            True
        """
        if self._master_key is None:
            raise ValueError("Key not set during initialisation")
        return iter(self(self._master_key))


def convert_to_vector(I, L):
    r"""
    Convert ``I`` to a bit vector of length ``L``.

    If ``I`` is list-like its elements will be treated as binary if
    ``len(I) == L`` or as hexadecimal if ``4*len(I) == L``.

    INPUT:

    - ``I`` -- integer or list-like

    - ``L`` -- integer; the desired length of the ouput

    OUTPUT:

    - A ``L``-bit vector representation of ``I``

    - A tuple containing the type of ``I`` and if necessary its representation
      e.g. "bin" or "hex"

    EXAMPLES::

    sage: from sage.crypto.block_cipher.present import convert_to_vector
    sage: convert_to_vector("0x1F", 8)
    ((1, 1, 1, 1, 1, 0, 0, 0), (<type 'str'>, None))
    sage: convert_to_vector(["0x1","0xF"], 8)
    ((1, 0, 0, 0, 1, 1, 1, 1), (<type 'list'>, 'hex'))
    sage: convert_to_vector(0x1F, 8)
    ((1, 1, 1, 1, 1, 0, 0, 0), (<type 'sage.rings.integer.Integer'>, None))
    sage: v = vector(GF(2), 4, [1,0,1,0])
    sage: convert_to_vector(v, 4)
    ((1, 0, 1, 0),
     (<type 'sage.modules.vector_mod2_dense.Vector_mod2_dense'>, 'bin'))
    sage: convert_to_vector(0x1F, 9)
    ((1, 1, 1, 1, 1, 0, 0, 0, 0), (<type 'sage.rings.integer.Integer'>, None))
    sage: convert_to_vector(["1","0xF"], 9)
    Traceback (most recent call last):
    ...
    ValueError: ['1', '0xF'] can not be converted to bit vector of length 9
    """
    try:
        state = vector(GF(2), L, ZZ(I).digits(2, padto=L))
        return state, (type(I), None)
    except TypeError:
        # ignore the error and try list-like types
        pass
    if len(I) == L:
        state = vector(GF(2), L, ZZ(list(I), 2).digits(2, padto=L))
        rep = "bin"
    elif 4*len(I) == L:
        state = vector(GF(2), L, ZZ(list(I), 16).digits(2, padto=L))
        rep = "hex"
    else:
        raise ValueError("%s can not be converted to bit vector of "
                         "length %s" % (I, L))
    return state, (type(I), rep)


def convert_vector_to_type(V, T):
    r"""
    Convert the bit vector ``B`` to something of type ``T``.

    If ``I`` is list-like its elements will be treated as binary if
    ``len(I) == L`` or as hexadecimal if ``4*len(I) == L``.

    INPUT:

    - ``V`` -- bit vector

    - ``T`` -- tuple; A tuple containing the type of ``I`` and if necessary
      its representation e.g. "bin" or "hex"


    OUTPUT:

    - A representation of V of type ``T``


    EXAMPLES::

    sage: from sage.crypto.block_cipher.present import convert_vector_to_type
    sage: v = vector(GF(2), [0,1,0,0])
    sage: t = (type(0xF), None)
    sage: convert_vector_to_type(v, t) == 0x2
    True
    sage: t = (type(v), None)
    sage: convert_vector_to_type(v, t) == v
    True
    sage: t = (type([]), "bin")
    sage: convert_vector_to_type(v, t) == list(v)
    True
    sage: t = (type([]), "oct")
    sage: convert_vector_to_type(v, t)
    Traceback (most recent call last):
    ...
    ValueError: can not convert `V` to <type 'list'>
    """
    # TODO there must be a better way than to compare str(T) to the types
    if str(T[0]) == "<type 'sage.modules.vector_mod2_dense.Vector_mod2_dense'>":
        return V
    elif str(T[0]) == "<type 'sage.rings.integer.Integer'>":
        return ZZ(list(V), 2)
    elif str(T[0]) == "<type 'list'>":
        if T[1] == "bin":
            return list(V)
        elif T[1] == "hex":
            return ZZ(list(V), 2).digits(16, padto=len(V)/4)
        else:
            raise ValueError("can not convert `V` to %s" % T[0])
    else:
        raise ValueError("can not convert `V` to %s" % T[0])
