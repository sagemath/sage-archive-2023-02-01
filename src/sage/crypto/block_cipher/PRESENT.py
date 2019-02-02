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


class PRESENT(SageObject):
    r"""
    This class implements PRESENT described in [BKLPPRSV2007]_.
    """

    def __init__(self, keysize=80):
        r"""
        Construct an instance of PRESENT.

        INPUT:

        - ``keysize`` -- (default: ``80``) the size of the keys that will be
          used in bits. It must be either 80 or 128.
        """
        self.keysize = 80

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
        pass

    def decrypt(self, C, K):
        r"""
        Return an plaintext corresponding to the ciphertext ``C``,
        using PRESENT decryption with key ``K``.
        """
        pass

    def encrypt(self, P, K):
        r"""
        Return an ciphertext corresponding to the plaintext ``P``,
        using PRESENT encryption with key ``K``.
        """
        pass
