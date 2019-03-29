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
    raise NotImplementedError
