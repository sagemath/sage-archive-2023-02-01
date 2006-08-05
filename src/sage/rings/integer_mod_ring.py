r"""
Ring $\Z/n\Z$ of integers modulo $n$

EXAMPLES:
    sage: R = Integers(97)
    sage: a = R(5)
    sage: a**100000000000000000000000000000000000000000000000000000000000000
    61

AUTHORS
    -- William Stein (initial code)
    -- David Joyner (2005-12-22): most examples
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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

import random
import weakref

import arith
import commutative_ring
import field
import integer_mod
import integer
import integer_ring
import rational
import quotient_ring
import ideal

import sage.interfaces.all

_objsIntegerModRing = {}
def IntegerModRing(order=0, check_prime=True):
    r"""
    INPUT:
        order -- integer (default: 0)
        check_prime -- bool (default: True); if False do not test for
                       primality of the order in constructing the
                       residue class ring (thus always constructing
                       the generic integer_mod_ring that doesn't
                       have special field functionality and implementation).
                       Do this if the modulus is huge (thousands of digits).

    EXAMPLES:
        sage: IntegerModRing(15)
        Ring of integers modulo 15

    The following example illustrates the \code{check_prime} option.
    Without it, just defining R would take a very long time.

        sage: n = 5*2^23473+1
        sage: len(str(n))
        7067
        sage: R = IntegerModRing(n, check_prime=False)
        sage: type(R)
        <class 'sage.rings.integer_mod_ring.IntegerModRing_generic'>

    Note that you can also user \code{Integers}, which is a synonym
    for \code{IntegerModRing}.
        sage: Integers(18)
        Ring of integers modulo 18

    """
    if order == 0:
        return integer_ring.IntegerRing()
    global _objsIntegerModRing
    if _objsIntegerModRing.has_key(order):
        x = _objsIntegerModRing[order]()
        if x != None: return x
    if check_prime and arith.is_prime(order):
        R = IntegerModRing_field(order)
    else:
        R = IntegerModRing_generic(order)
    _objsIntegerModRing[order] = weakref.ref(R)
    return R

def is_IntegerModRing(x):
    """
    Return True if x is an integer modulo ring.

    EXAMPLES:
        sage: R = IntegerModRing(17)
        sage: is_IntegerModRing(R)
        True
        sage: is_IntegerModRing(GF(13))
        True
        sage: is_IntegerModRing(GF(4))
        False
        sage: is_IntegerModRing(10)
        False
        sage: is_IntegerModRing(ZZ)
        False
    """
    return isinstance(x, IntegerModRing_generic)

class IntegerModRing_generic(quotient_ring.QuotientRing_generic):
    """
    The ring of integers modulo N, e.g., when N is prime
    this is a prime finite field.

    EXAMPLES:
        sage: R = IntegerModRing(97)
        sage: a = R(5)
        sage: a**(10^62)
        61
    """
    def __init__(self, order):
        """
        Create with the command
              IntegerModRing(order)

        INPUT:
            order -- an integer > 1

        OUTPUT:
            IntegerModRing -- the ring of integers modulo N.

        EXAMPLES:

        First we compute with integers modulo $17$.
            sage: FF = IntegerModRing(17)
            sage: FF
            Ring of integers modulo 17
            sage: FF.is_field()
            True
            sage: FF.characteristic()
            17
            sage: FF.order()
            17
            sage: gens = FF.unit_gens()
            sage: a = gens[0]
            sage: a
            3
            sage: a.is_square()
            False
            sage: def pow(i): return a**i
            sage: [pow(i) for i in range(16)]
            [1, 3, 9, 10, 13, 5, 15, 11, 16, 14, 8, 7, 4, 12, 2, 6]

        Next we compute with the integers modulo $16$.
            sage: Z16 = IntegerModRing(16)
            sage: Z16.is_field()
            False
            sage: Z16.order()
            16
            sage: Z16.characteristic()
            16
            sage: gens = Z16.unit_gens()
            sage: gens
            [15, 5]
            sage: a = gens[0]
            sage: b = gens[1]
            sage: def powa(i): return a**i
            sage: def powb(i): return b**i
            sage: gp_exp = FF.unit_group_exponent()
            sage: gp_exp
            16
            sage: [powa(i) for i in range(15)]
            [1, 15, 1, 15, 1, 15, 1, 15, 1, 15, 1, 15, 1, 15, 1]
            sage: [powb(i) for i in range(15)]
            [1, 5, 9, 13, 1, 5, 9, 13, 1, 5, 9, 13, 1, 5, 9]
            sage: a.multiplicative_order()
            2
            sage: b.multiplicative_order()
            4

        Saving and loading:
            sage: R = Integers(100000)
            sage: loads(R.dumps()) == R
            True
        """
        ZZ = integer_ring.IntegerRing()
        order = ZZ(order)
        if order <= 0:
            raise ZeroDivisionError, "order must be positive"
        self.__order = order
        self.__unit_group_exponent = None
        self.__factored_order = None
        quotient_ring.QuotientRing_generic.__init__(self, ZZ, ZZ.ideal(order))

    def is_finite(self):
        """
        EXAMPLES:
            sage: R = IntegerModRing(18)
            sage: R.is_finite()
            True
        """
        return True

    def is_field(self):
        """
        EXAMPLES:
            sage: R = IntegerModRing(18)
            sage: R.is_field()
            False
            sage: FF = IntegerModRing(17)
            sage: FF.is_field()
            True
        """
        return False

    def factored_order(self):
        """
        EXAMPLES:
            sage: R = IntegerModRing(18)
            sage: FF = IntegerModRing(17)
            sage: R.factored_order()
            2 * 3^2
            sage: FF.factored_order()
            17
        """
        if self.__factored_order != None:
            return self.__factored_order
        self.__factored_order = arith.factor(self.__order, int_=True)
        return self.__factored_order

    def characteristic(self):
        """
        EXAMPLES:
            sage: R = IntegerModRing(18)
            sage: FF = IntegerModRing(17)
            sage: FF.characteristic()
            17
            sage: R.characteristic()
            18
        """
        return self.__order

    def _repr_(self):
        return "Ring of integers modulo %s"%self.__order

    def modulus(self):
        r"""
        Return the polynomial $x - 1$ over this ring.

        \note{This function exists for consistency with the
        finite-field modulus function.}

        EXAMPLES:
            sage: R = IntegerModRing(18)
            sage: R.modulus()
            x + 17
            sage: R = IntegerModRing(17)
            sage: R.modulus()
            x + 16
        """
        try:
            return self.__modulus
        except AttributeError:
            from polynomial_ring import PolynomialRing
            x = PolynomialRing(self).gen()
            self.__modulus = x - 1
            return self.__modulus

    def order(self):
        return self.__order

    def __call__(self, x, construct=False):
        global TypeError
        if sage.interfaces.all.is_GapElement(x):
            import finite_field
            try:
                return finite_field.gap_to_sage(x, self)
            except (ValueError, IndexError, TypeError):
                raise TypeError, "error coercing %s to finite field"%x
        return integer_mod.IntegerMod(self, x, construct=construct)
        try:
            return integer_mod.IntegerMod(self, x, construct=construct)
        except (RuntimeError, TypeError), msg:
            raise TypeError, "no way to coerce %s of type %s into %s."%(
                x, type(x), self)

    def _coerce_(self, x):
        r"""
        Canonical coercion.

        EXAMPLES:
            sage: R = IntegerModRing(17)
            sage: a = R(3)
            sage: b = R._coerce_(3)
            sage: b
            3
            sage: a==b
            True

        This is allowed:
            sage: R(2/3)
            12

        But this is not, since there is no (canonical or not!) ring homomorphism
        from $\Q$ to $\GF(17)$.

            sage: R._coerce_(2/3)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion of 2/3 into Ring of integers modulo 17        """
        if isinstance(x, (int, long, integer.Integer)):
            return integer_mod.IntegerMod(self, x)
        if isinstance(x, integer_mod.IntegerMod) and x.parent().order() % self.modulus() == 0:
            return integer_mod.IntegerMod(self, x)
        raise TypeError, "no canonical coercion of %s into %s"%(x, self)

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: F = GF(11)
            sage: F
            Finite Field of size 11
            sage: R = IntegerModRing(11)
            sage: R == F
            True
        """
        if not isinstance(other, IntegerModRing_generic):
            return -1
        if self.__order < other.__order:
            return -1
        elif self.__order > other.__order:
            return 1
        return 0

    ## This conflicts with polynomial ring constructor
##     def __getitem__(self, n):
##         if n < 0 or n >= self.order():
##             raise IndexError
##         return self(n)

    # The following __unit_gens functions are here since I just factored
    # them out from the unit_gens function.  They are only called by
    # the unit_gens function.
    def __unit_gens_primecase(self, p):
        if p==2:
            return integer_mod.Mod(1,p)
        P = arith.prime_divisors(p-1)
        ord = integer.Integer(p-1)
        one = integer_mod.Mod(1,p)
        x = 2
        while x < p:
            generator = True
            z = integer_mod.Mod(x,p)
            for q in P:
                if z**(ord/q) == one:
                    generator = False
                    break
            if generator:
                return z
            x += 1
        #end for
        assert False, "didn't find primitive root for p=%s"%p

    def __unit_gens_primepowercase(self, p, r):
        r"""
        Find smallest generator for $(\Z/p^r\Z)^*$.
        """
        if r==1:
            return [self.__unit_gens_primecase(p)]
        elif p==2:
            if r==1:
                return []
            elif r==2:
                return [integer_mod.Mod(-1,2**r)]
            elif r>=3:
                pr=2**r
                a = integer_mod.Mod(5, pr)
                return [integer_mod.Mod(-1,pr), a]
            assert False, "p=2, r=%s should be >=1"%r
        else:  # odd prime
            pr = p**r
            R = IntegerModRing(pr)
            x = R(self.__unit_gens_primecase(p).lift())
            n = p**(r-2)*(p-1)
            one = integer_mod.Mod(1,pr)
            for b in range(0,p):
                z = x+R(b*p)
                if z**n != one:
                    a = integer_mod.Mod(z,pr)
                    return [a]
            assert False, "p=%s, r=%s, couldn't find generator"%(p,r)

    def unit_gens(self):
        r"""
        Returns generators for the unit group $(\Z/N\Z)^*$.

        We compute the list of generators using a deterministic
        algorithm, so the generators list will always be the same.
        Each generator corresponds to a prime divisor of $N$ (or
        possibly two prime divisors for p=2).

        INPUT: (none)
        OUTPUT:
            list -- a list of elements of self

        EXAMPLES:
            sage: R = IntegerModRing(18)
            sage: R.unit_gens()
            [1, 11]
            sage: R = IntegerModRing(17)
            sage: R.unit_gens()
            [3]
        """
        try:
            return self.__unit_gens
        except AttributeError:
            self.__unit_gens = []
        n = self.__order
        if n == 1:
            self.__unit_gens = [self(1)]
            return self.__unit_gens
        for p,r in self.factored_order():
            m = n/(p**r)
            for g in self.__unit_gens_primepowercase(p, r):
                x = g.crt(integer_mod.Mod(1,m))
                self.__unit_gens.append(x)
        return self.__unit_gens

    def unit_group_exponent(self):
        """
        EXAMPLES:
            sage: R = IntegerModRing(17)
            sage: R.unit_group_exponent()
            16
            sage: R = IntegerModRing(18)
            sage: R.unit_group_exponent()
            6
        """
        if self.__unit_group_exponent != None:
            return self.__unit_group_exponent
        a = []
        for p, r in self.factored_order():
            if p != 2:
                a.append((p-1)*(p**(r-1)))   # phi(p**r)
            else:  # p=2
                if r==2:
                    a.append(2)
                elif r>2:
                    a.append(2**(r-2))
            #endif
        #endfor
        self.__unit_group_exponent = int(arith.LCM(a))
        return self.__unit_group_exponent

    def unit_group_order(self):
        """
        Return the order of the unit group of this residue class ring.

        EXAMPLES;
            sage: R = Integers(500)
            sage: R.unit_group_order()
            200
        """
        return arith.euler_phi(self.order())

    def random_element(self, bound=None):
        """
        Return a random element of this ring.

        If bound is not None, return the coercion of an integer in the
        interval [-bound, bound] into this ring.


        EXAMPLES:
            sage: R = IntegerModRing(18)
            sage: R.random_element()
            15
        """
        if not (bound is None):
            return commutative_ring.CommutativeRing.random_element(self, bound)
        a = random.randrange(0,self.order()-1)
        return self(a)

    #######################################################
    # Suppose for interfaces
    #######################################################
    def _gap_init_(self):
        """
        EXAMPLES:
           sage: R = Integers(12345678900)
           sage: R
           Ring of integers modulo 12345678900
           sage: gap(R)
           (Integers mod 12345678900)
        """
        return 'ZmodnZ(%s)'%self.order()

    def _magma_init_(self):
        """
        EXAMPLES:
            sage: R = Integers(12345678900)
            sage: R
            Ring of integers modulo 12345678900
            sage: magma(R)                                          # optional
            Residue class ring of integers modulo 12345678900
        """
        return 'Integers(%s)'%self.order()

class IntegerModRing_field(field.Field, IntegerModRing_generic):
    def __init__(self, order):
        IntegerModRing_generic.__init__(self, order)


Zmod = IntegerModRing
Integers = IntegerModRing

## def GF(p):
##     """
##     EXAMPLES:
##         sage: F = GF(11)
##         sage: F
##         Finite field of size 11
##     """
##     if not arith.is_prime(p):
##         raise NotImplementedError, "Only prime fields currently implemented."
##     return IntegerModRing(p)

def crt(v):
    """
    INPUT: v -- (list) a lift of elements of rings.IntegerMod(n), for
                 various coprime moduli n.
    """
    if len(v) == 0:
        return IntegerModRing(1)(1)
    x = v[0]
    for i in range(1,len(v)):
        x = x.crt(v[i])
    return x

