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

    def __call__(self, B, K, algorithm="encrypt"):
        r"""
        This class implements PRESENT described in [BKLPPRSV2007]_.
        """
        pass

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        PRESENT objects are the same if ???
        """
        pass

    def __repr__(self):
        r"""
        A string representation of this PRESENT.
        """
        return "PRESENT block cipher with %s-bit keys" % self.keysize

    def decrypt(self, C, K):
        r"""
        Return an plaintext corresponding to the ciphertext ``C``,
        using PRESENT decryption with key ``K``.
        """
        raise NotImplementedError("Decryption is not implemented yet!")

    def encrypt(self, P, K):
        r"""
        Return an ciphertext corresponding to the plaintext ``P``,
        using PRESENT encryption with key ``K``.

        INPUT:

        - ``P`` -- The plaintext that will be encrypted.

        - ``K`` -- a string of 16 or 32 hex digits; The key that will be used
          to encrypt ``P``.

        OUTPUT:

        - The ciphertext corresponding to ``P``, obtained using the key ``K``.
        """
        raise NotImplementedError("Encryption is not implemented yet!")
