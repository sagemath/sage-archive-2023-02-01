"""
Cryptosystems.
"""

#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.structure.parent_old as parent_old
from sage.sets.set import Set_generic

class Cryptosystem(parent_old.Parent, Set_generic):
      r"""
      A cryptosystem is a pair of maps

      .. math::

                 E : {\mathcal K} \rightarrow {\rm Hom}({\mathcal M},{\mathcal C})



      .. math::

                 D : {\mathcal K} \rightarrow {\rm Hom}({\mathcal C},{\mathcal M})


      where `{\mathcal K}` is the keyspace,
      `{\mathcal M}` is the plaintext or message space, and
      `{\mathcal C}` is the ciphertext space. In many instances
      `{\mathcal M} = {\mathcal C}` and the images will lie in
      `{\rm Aut}({\mathcal M})`. An element of the image of
      `E` is called a cipher.

      We may assume that `E` and `D` are injective, hence
      identify a key `K` in `{\mathcal K}` with its image
      `E_K := E(K)` in
      `\mathrm{Hom}({\mathcal M},{\mathcal C})`.

      The cryptosystem has the property that for every encryption key
      `K_1` there is a decryption key `K_2` such that
      `D_{K_2} \circ E_{K_1}`. A cryptosystem with the
      property that `K := K_2 = K_1`, is called a symmetric
      cryptosystem. Otherwise, if the key `K_2 \ne K_1`, nor is
      `K_2` easily derived from `K_1`, we call the
      cryptosystem asymmetric of public key. In that case, `K_1`
      is called the public key and `K_2` is called the private
      key.
      """
      def __init__(self, plaintext_space, ciphertext_space, key_space, block_length = 1, period = None):
	  self._cipher_domain = plaintext_space
	  self._cipher_codomain = ciphertext_space
	  self._key_space = key_space
	  self._block_length = block_length
	  self._period = period

      def __eq__(self,right):
            return type(self) == type(right) and  \
                   self._cipher_domain == right._cipher_domain and \
                   self._cipher_codomain == right._cipher_codomain and \
                   self._key_space ==  right._key_space and \
                   self._block_length == right._block_length and \
                   self._period == right._period


      def plaintext_space(self):
	  return self._cipher_domain

      def cipher_domain(self):
	  return self._cipher_domain

      def ciphertext_space(self):
	  return self._cipher_codomain

      def cipher_codomain(self):
	  return self._cipher_codomain

      def key_space(self):
	  return self._key_space

      def block_length(self):
          return self._block_length

      def period(self):
	  if self._period is None:
	      raise TypeError, "Argument has no associated period."
	  return self._period

class SymmetricKeyCryptosystem(Cryptosystem):
      """

      """

class PublicKeyCryptosystem(Cryptosystem):
      """

      """


