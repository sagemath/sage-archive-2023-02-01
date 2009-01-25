"""
The set of prime numbers.

TESTS:
    sage: loads(dumps(Primes())) == Primes()
    True
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

from sage.rings.all import Integer, ZZ, infinity

from set import Set_generic

class Primes_class(Set_generic):
    """
    The set of prime numbers.

    EXAMPLES::

        sage: P = Primes(); P
        Set of all prime numbers: 2, 3, 5, 7, ...
        sage: loads(P.dumps()) == P
        True
    """
    def __init__(self, proof=True):
        """
        EXAMPLES::

            sage: P = Primes(); P
            Set of all prime numbers: 2, 3, 5, 7, ...
        """
        self.__proof = proof

    def cardinality(self):
        """
        There is no largest prime number, so we say the set
        has infinite cardinality.

        EXAMPLES::

            sage: P = Primes()
            sage: P.cardinality()
            +Infinity
        """
        return infinity

    def __len__(self):
        """
        TESTS::

            sage: len(Primes())
            Traceback (most recent call last):
            ...
            TypeError
        """
        raise TypeError


    def __cmp__(self, right):
        """
        The set of primes can be compared to various things,
        but is only equal to itself.

        EXAMPLES::

            sage: P = Primes()
            sage: R = Primes()
            sage: P.__cmp__(R)
            0
            sage: P == R
            True
            sage: P != R
            False
            sage: Q=[1,2,3]
            sage: Q != P        # indirect doctest
            True
            sage: R.<x>=ZZ[]
            sage: P!=x^2+x
            True

        Make sure changing order changes the comparison with something
        of a different type::

            sage: cmp('foo', Primes()) != cmp(Primes(), 'foo')
            True
        """
        if isinstance(right, Primes_class):
            return 0
        return -1

    def _repr_(self):
        """
        Representation of the set of primes.

        EXAMPLES::

            sage: P = Primes(); P
            Set of all prime numbers: 2, 3, 5, 7, ...
        """
        return "Set of all prime numbers: 2, 3, 5, 7, ..."

    def __iter__(self):
        """
        Iterator for the set of primes.  This is an infinite set,
        so USE WITH CAUTION!  That is, do not do things like
        ``[p for p in Primes()]``.

        EXAMPLES::

            sage: P = Primes()
            sage: iter(P).next()
            2
        """
        p = Integer(2)
        while True:
            yield p
            p = p.next_prime(self.__proof)

    def __contains__(self, x):
        """
        Checks whether an object is a prime number.  If
        it is not an integer, returns False.

        EXAMPLES::

            sage: P = Primes()
            sage: 5 in P
            True
            sage: 100 in P
            False
            sage: 1.5 in P
            False
            sage: e in P
            False
        """
        try:
            if not x in ZZ:
                return False
            return ZZ(x).is_prime()
        except TypeError:
            return False

the_set_of_primes = {True: Primes_class(proof=True), False: Primes_class(proof=False)}

def Primes(proof=True):
    """
    Return the set of prime numbers.

    EXAMPLES::

        sage: P = Primes(); P
        Set of all prime numbers: 2, 3, 5, 7, ...

    We show various methods about the primes::

        sage: P.cardinality()
        +Infinity
        sage: R = Primes()
        sage: P == R
        True
        sage: 5 in P
        True
        sage: 100 in P
        False
    """
    return the_set_of_primes[proof]
