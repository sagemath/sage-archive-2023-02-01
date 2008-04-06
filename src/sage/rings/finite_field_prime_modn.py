"""
Finite Prime Fields

AUTHORS:
     -- William Stein: initial version
     -- Martin Albrecht (2008-01): refactoring

TESTS:
    sage: k = GF(3)
    sage: loads(dumps(k)) == k
    True
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Martin Albrecht <malb@informatik.uni-bremen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sys

from ring import FiniteField as FiniteField_generic
from sage.structure.parent_gens import normalize_names, ParentWithGens

import integer_mod_ring
import integer
import rational
import integer_mod
import arith


class FiniteField_prime_modn(FiniteField_generic, integer_mod_ring.IntegerModRing_generic):
    def __init__(self, p, name=None):
        """
        Return a new finite field of order $p$ where $p$ is prime.

        INPUT:
            p -- an integer >= 2
            name -- ignored

        EXAMPLES:
            sage: FiniteField(3)
            Finite Field of size 3

            sage: FiniteField(next_prime(1000))
            Finite Field of size 1009
        """
        p = integer.Integer(p)
        if not arith.is_prime(p):
            raise ArithmeticError, "p must be prime"
        integer_mod_ring.IntegerModRing_generic.__init__(self, p)
        import sage.structure.factorization as factorization
        self._IntegerModRing_generic__factored_order = factorization.Factorization([(p,1)], integer.Integer(1))
        self._kwargs = {}
        self.__char = p
        self.__gen = self(1)  # self(int(pari.pari(p).znprimroot().lift()))
        ParentWithGens.__init__(self, self, ('x',), normalize=False)

    def __cmp__(self, other):
        r"""
        Compare \code{self} with \code{other}. Two finite prime fields
        are considered equal if their characteristic is equal.

        EXAMPLE:
            sage: K = FiniteField(3)
            sage: copy(K) == K
            True
            sage: copy(K) is K
            False
        """
        if not isinstance(other, FiniteField_prime_modn):
            return cmp(type(self), type(other))
        return cmp(self.__char, other.__char)

    def _is_valid_homomorphism_(self, codomain, im_gens):
        """
        This is called implicitly by the hom constructor.

        EXAMPLES:
            sage: k = GF(73^2,'a')
            sage: f = k.modulus()
            sage: r = f.change_ring(k).roots()
            sage: k.hom([r[0][0]]) # indirect doctest
            Ring endomorphism of Finite Field in a of size 73^2
              Defn: a |--> 72*a + 3
        """
        try:
            return im_gens[0] == codomain._coerce_(self.gen(0))
        except TypeError:
            return False

    def _coerce_impl(self, x):
        """
        This is called implicitly by arithmetic methods.

        EXAMPLES:
            sage: k = GF(7)
            sage: e = k(6)
            sage: e * 2 # indirect doctest
            5
            sage: 12 % 7
            5
        """
        if isinstance(x, (int, long, integer.Integer)):
            return self(x)
        if isinstance(x, integer_mod.IntegerMod_abstract) and \
               x.parent().characteristic() == self.characteristic():
            return self(x)
        raise TypeError, "no canonical coercion of x"

    def characteristic(self):
        r"""
        Return the characteristic of \code{self}.

        EXAMPLE:
            sage: k = GF(7)
            sage: k.characteristic()
            7
        """
        return self.__char

    def modulus(self):
        """
        Return the minimal polynomial of self, which is allways $x - 1$.

        EXAMPLE:
            sage: k = GF(199)
            sage: k.modulus()
            x + 198
        """
        try:
            return self.__modulus
        except AttributeError:
            x = self['x'].gen()
            self.__modulus = x - 1
        return self.__modulus

    def is_prime_field(self):
        """
        Return True

        EXAMPLE:
            sage: k.<a> = GF(3)
            sage: k.is_prime_field()
            True

            sage: k.<a> = GF(3^2)
            sage: k.is_prime_field()
            False
        """
        return True

    def polynomial(self, name=None):
        """
        Returns the polynomial \var{name}.

        EXAMPLE:
            sage: k.<a> = GF(3)
            sage: k.polynomial()
            x
        """
        if name is None:
            name = self.variable_name()
        try:
            return self.__polynomial[name]
        except  AttributeError:
            from sage.rings.finite_field import FiniteField
            R = FiniteField(self.characteristic())[name]
            f = self[name]([0,1])
            try:
                self.__polynomial[name] = f
            except (KeyError, AttributeError):
                self.__polynomial = {}
                self.__polynomial[name] = f
            return f

    def order(self):
        """
        Return the order of this finite field.

        EXAMPLE:
            sage: k = GF(5)
            sage: k.order()
            5
        """
        return self.__char

    def gen(self, n=0):
        """
        Return generator of this finite field.

        EXAMPLES:
            sage: k = GF(13)
            sage: k.gen()
            1
            sage: k.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: only one generator
        """
        if n != 0:
            raise IndexError, "only one generator"
        return self.__gen

    def __iter__(self):
        """
        EXAMPLES:
            sage: list(GF(7))
            [0, 1, 2, 3, 4, 5, 6]

        We can even start iterating over something that would be too big
        to actually enumerate:
            sage: K = GF(next_prime(2^256))
            sage: all = iter(K)
            sage: all.next()
            0
            sage: all.next()
            1
            sage: all.next()
            2
        """
        yield self(0)
        i = one = self(1)
        while i:
            yield i
            i += one

    def degree(self):
        """
        Returns the degree of the finite field, which is a positive
        integer.

        EXAMPLES:
            sage: FiniteField(3).degree()
            1
        """
        return 1
