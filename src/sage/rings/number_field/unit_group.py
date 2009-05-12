"""
The unit group of a number field.

EXAMPLES:
    sage: x = polygen(QQ)
    sage: K.<a> = NumberField(x^4-8*x^2+36)
    sage: UK = UnitGroup(K); UK
    Unit group with structure C4 x Z of Number Field in a with defining polynomial x^4 - 8*x^2 + 36

    # The first generator is a primitive root of unity in the field:
    sage: UK.gens() # random
    [1/12*a^3 - 1/6*a, 1/24*a^3 + 1/4*a^2 - 1/12*a - 1]
    sage: [u.multiplicative_order() for u in UK.gens()]
    [4, +Infinity]

    sage: UK.rank()
    1
    sage: UK.ngens()
    2

    # Units in the field can be converted into elements of the unit
    # group represented as elements of an abstract multiplicative
    # group:
    sage: UK(1)
    1
    sage: UK(-1)
    u0^2
    sage: [UK(u) for u in (x^4-1).roots(K,multiplicities=False)]
    [1, u0^2, u0, u0^3]

    sage: UK.fundamental_units() # random
    [1/24*a^3 + 1/4*a^2 - 1/12*a - 1]
    sage: UK.torsion_generator()
    1/12*a^3 - 1/6*a
    sage: UK.zeta_order()
    4
    sage: UK.roots_of_unity()
    [1/12*a^3 - 1/6*a, -1, -1/12*a^3 + 1/6*a, 1]

    # exp and log functions provide maps between units as field
    # elements and exponent vectors with respect to the generators:

    sage: u = UK.exp([13,10]); u # random
    -41/8*a^3 - 55/4*a^2 + 41/4*a + 55
    sage: UK.log(u)
    [1, 10]
    sage: u = UK.fundamental_units()[0]
    sage: [UK.log(u^k) == [0,k] for k in range(10)]
    [True, True, True, True, True, True, True, True, True, True]
    sage: all([UK.log(u^k) == [0,k] for k in range(10)])
    True

    sage: K.<a> = NumberField(x^5-2,'a')
    sage: UK = UnitGroup(K)
    sage: UK.rank()
    2
    sage: UK.fundamental_units()
    [a^3 + a^2 - 1, a - 1]

    # A relative number field example:
    sage: L.<a, b> = NumberField([x^2 + x + 1, x^4 + 1])
    sage: UL = L.unit_group(); UL
    Unit group with structure C24 x Z x Z x Z of Number Field in a with defining polynomial x^2 + x + 1 over its base field
    sage: UL.gens() # random
    [-b^3*a - b^3, -b^3*a + b, (-b^3 - b^2 - b)*a - b - 1, (-b^3 - 1)*a - b^2 + b - 1]
    sage: UL.zeta_order()
    24
    sage: UL.roots_of_unity()
    [-b^3*a - b^3,
    -b^2*a,
    b,
    a + 1,
    -b^3*a,
    b^2,
    b*a + b,
    a,
    b^3,
    b^2*a + b^2,
    b*a,
    -1,
    b^3*a + b^3,
    b^2*a,
    -b,
    -a - 1,
    b^3*a,
    -b^2,
    -b*a - b,
    -a,
    -b^3,
    -b^2*a - b^2,
    -b*a,
    1]


    # A relative extension example, which worked thanks to the code review by F.W.Clarke:
    sage: PQ.<X> = QQ[]
    sage: F.<a, b> = NumberField([X^2 - 2, X^2 - 3])
    sage: PF.<Y> = F[]
    sage: K.<c> = F.extension(Y^2 - (1 + a)*(a + b)*a*b)
    sage: K.unit_group()
    Unit group with structure C2 x Z x Z x Z x Z x Z x Z x Z of Number Field in c with defining polynomial Y^2 + (-2*b - 3)*a - 2*b - 6 over its base field

AUTHOR:
    -- John Cremona
"""
#*****************************************************************************
#       Copyright (C) 2009 William Stein, John Cremona
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

from sage.groups.abelian_gps.abelian_group import AbelianGroup_class
from sage.groups.abelian_gps.abelian_group_element import AbelianGroupElement

from sage.structure.sequence import Sequence
from sage.structure.element import MultiplicativeGroupElement
from sage.structure.proof.proof import get_flag
from sage.libs.all import pari, pari_gen
from sage.misc.misc import prod
from sage.rings.integer_ring import ZZ

class UnitGroup(AbelianGroup_class):
    """
    The unit group of a number field.
    """
    def __init__(self, number_field, proof=True):
        """
        Create a unit group of a number field.

        INPUT:
            number_field - a number field
            proof - boolean (default True): proof flag

        The proof flag is passed to pari via the pari_bnf() function
        which computes the unit group.  See the documentation for the
        number_field module.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2-38)
            sage: UK = K.unit_group(); UK
            Unit group with structure C2 x Z of Number Field in a with defining polynomial x^2 - 38
            sage: UK.gens()
            [-1, 6*a - 37]

            sage: K.<a> = QuadraticField(-3)
            sage: UK = K.unit_group(); UK
            Unit group with structure C6 of Number Field in a with defining polynomial x^2 + 3
            sage: UK.gens()
            [-1/2*a + 1/2]

            sage: K.<z> = CyclotomicField(13)
            sage: UK = K.unit_group(); UK
            Unit group with structure C26 x Z x Z x Z x Z x Z of Cyclotomic Field of order 13 and degree 12
            sage: UK.gens() # random
            [-z^11, z^5 + z^3, z^6 + z^5, z^9 + z^7 + z^5, z^9 + z^5 + z^4 + 1, z^5 + z]
            """
        proof = get_flag(proof, "number_field")
        K = number_field
        pK = K.pari_bnf(proof)
        self.__number_field = K

        # compute the units via pari:
        fu = [K(u) for u in pK.bnfunit()]

        # compute a torsion generator and pick the 'simplest' one:
        n, z = pK.nfrootsof1()
        n = ZZ(n)
        self.__ntu = n

        # For an absolute field we can now set z = K(z), but this
        # does not work for relative fields, so we work harder:
        if K.is_absolute():
            z = K(z)
        else:
            zk = pK.getattr('zk')
            z = z.mattranspose()
            cc = [z[0,i] for i in range(z.ncols())]
            z = sum([K(c*d) for d,c in zip(zk,cc)])

        # If we replaced z by another torsion generator we would need
        # to allow for this in the dlog function!  So we do not.

        # Store the actual generators (torsion first):
        gens = [z] + fu
        self.__nfu = len(fu)
        self.__gens = Sequence(gens, immutable=True, universe=self, check=False)
        # Construct the abtract group:
        AbelianGroup_class.__init__(self, 1+len(fu), [n]+[0]*len(fu), 'u')


    def __call__(self, u):
        """
        Returns the abstract group element corresponding to the unit u.

        INPUT:
            u -- Any object from which an element of the unit group's
            number field K may be constructed; an error is raised if
            an element of K cannot be constructed from u, or if the
            element constructed is not a unit.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2-38)
            sage: UK = UnitGroup(K)
            sage: UK(1)
            1
            sage: UK(-1)
            u0
            sage: UK.gens()
            [-1, 6*a - 37]
            sage: UK.ngens()
            2
            sage: [UK(u) for u in UK.gens()]
            [u0, u1]
            sage: [UK(u).list() for u in UK.gens()]
            [[1, 0], [0, 1]]
            sage: UK(a)
            Traceback (most recent call last):
            ...
            ValueError: a is not a unit
        """
        K = self.__number_field

        try:
            u = K(u)
        except TypeError:
            raise ValueError, "%s is not an element of %s"%(u,K)
        if not u.is_integral() or u.norm().abs() != 1:
            raise ValueError, "%s is not a unit"%u
        m = K.pari_bnf().bnfisunit(pari(u)).mattranspose()
        # convert column matrix to a list:
        m = [ZZ(m[0,i].python()) for i in range(m.ncols())]
        # NB pari puts the torsion at the end!
        m.insert(0,m.pop())
        return AbelianGroupElement(self, m)

    def _coerce_impl(self, x):
        """
        Canonical coercion of x into this unit group.

        EXAMPLES:

        """
        return self(x)

    def gens(self):
        """
        Return generators for the unit group, as a list.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^4 + 23)
            sage: K.unit_group().gens() # random
            [-1, 1/4*a^3 - 7/4*a^2 + 17/4*a - 19/4]
        """
        return self.__gens

    def ngens(self):
        """
        Return the number of generators of the unit group.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: U = NumberField(x^2 + x + 23899, 'a').unit_group(); U
            Unit group with structure C2 of Number Field in a with defining polynomial x^2 + x + 23899
            sage: U.ngens()
            1
        """
        return len(self.__gens)

    def rank(self):
        """
        Return the rank of the unit group.

        EXAMPLES:
        sage: K.<z> = CyclotomicField(13)
        sage: UnitGroup(K).rank()
        5
        """
        return len(self.__gens)-1

    def gen(self, i=0):
        """
        Return the i-th generator for this unit group.

        NOTE: i=0 gives the torsion generator, i.e. a primitive root of unity.

        EXAMPLES:
            sage: K.<z> = CyclotomicField(13)
            sage: UK = K.unit_group()
            sage: UK.ngens()
            6
            sage: UK.gen(0) # random
            -z^11
            sage: UK.gen(1) # random
            z^5 + z^3
            sage: UK.gen(2) # random
            z^6 + z^5
            sage: UK.gen(3) # random
            z^9 + z^7 + z^5
            sage: UK.gen(4) # random
            z^9 + z^5 + z^4 + 1
            sage: UK.gen(5) # random
            z^5 + z
        """
        if i < 0 or i >= len(self.__gens):
            raise IndexError
        return self.__gens[i]

    def _repr_(self):
        """
        Return string representation of this unit group.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: U = UnitGroup(NumberField(x^3 - 2, 'a'))
            sage: U
            Unit group with structure C2 x Z of Number Field in a with defining polynomial x^3 - 2
            sage: U._repr_()
            'Unit group with structure C2 x Z of Number Field in a with defining polynomial x^3 - 2'
        """
        return 'Unit group with structure %s of %s'%(
            self._group_notation(self.invariants()),
            self.number_field())

    def fundamental_units(self):
        """
        Return generators for the free part of the unit group, as a list.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^4 + 23)
            sage: U = UnitGroup(K)
            sage: U.fundamental_units()  # random
            [1/4*a^3 - 7/4*a^2 + 17/4*a - 19/4]
        """
        return self.__gens[1:]

    def roots_of_unity(self):
        """
        Return all the roots of unity in this unit group, primitive or not.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: K.<b> = NumberField(x^2+1)
            sage: U = UnitGroup(K)
            sage: zs = U.roots_of_unity(); zs
            [b, -1, -b, 1]
            sage: [ z**U.zeta_order() for z in zs ]
            [1, 1, 1, 1]
        """
        z = self.__gens[0]
        n = self.__ntu
        return [ z**k for k in range(1, n+1) ]

    def torsion_generator(self):
        """
        Return a generator for the torsion part of the unit group.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^4 - x^2 + 4)
            sage: U = UnitGroup(K)
            sage: U.torsion_generator() # random
            -1/4*a^3 - 1/4*a + 1/2
        """
        return self.__gens[0]

    def zeta_order(self):
        """
        Returns the order of the torsion part of the unit group.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^4 - x^2 + 4)
            sage: U = UnitGroup(K)
            sage: U.zeta_order()
            6
        """
        return self.__ntu

    def zeta(self, n=2, all=False):
        """
        Return one, or a list of all, primitive n-th root of unity in this unit group.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: K.<z> = NumberField(x^2 + 3)
            sage: U = UnitGroup(K)
            sage: U.zeta(1)
            1
            sage: U.zeta(2)
            -1
            sage: U.zeta(2, all=True)
            [-1]
            sage: U.zeta(3)
            -1/2*z - 1/2
            sage: U.zeta(3, all=True)
            [-1/2*z - 1/2, 1/2*z - 1/2]
            sage: U.zeta(4)
            Traceback (most recent call last):
            ...
            ValueError: n (=4) does not divide order of generator

            sage: r.<x> = QQ[]
            sage: K.<b> = NumberField(x^2+1)
            sage: U = UnitGroup(K)
            sage: U.zeta(4)
            b
            sage: U.zeta(4,all=True)
            [b, -b]
            sage: U.zeta(3)
            Traceback (most recent call last):
            ...
            ValueError: n (=3) does not divide order of generator
            sage: U.zeta(3,all=True)
            []

        """
        N = self.__ntu
        K = self.number_field()
        n = ZZ(n)
        if n <= 0:
            raise ValueError, "n (=%s) must be positive"%n
        if n == 1:
            if all:
                return [K(1)]
            else:
                return K(1)
        elif n == 2:
            if all:
                return [K(-1)]
            else:
                return K(-1)
        if n.divides(N):
            z = self.torsion_generator() ** (N//n)
            if all:
                return [z**i for i in n.coprime_integers(n)]
            else:
                return z
        else:
            if all:
                return []
            else:
                raise ValueError, "n (=%s) does not divide order of generator"%n

    def number_field(self):
        """
        Return the number field associated with this unit group.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: U = UnitGroup(NumberField(x^2 + 23, 'w')); U
            Unit group with structure C2 of Number Field in w with defining polynomial x^2 + 23
            sage: U.number_field()
            Number Field in w with defining polynomial x^2 + 23
        """
        return self.__number_field


    def log(self, u):
        """
        Return the exponents of the unit u with respect to group generators.

        INPUT:
            u -- Any object from which an element of the unit group's
            number field K may be constructed; an error is raised if
            an element of K cannot be constructed from u, or if the
            element constructed is not a unit.

        OUTPUT: a list of integers giving the exponents of u with
        respect to the unit group's basis.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: K.<z> = CyclotomicField(13)
            sage: UK = UnitGroup(K)
            sage: [UK.log(u) for u in UK.gens()]
            [[1, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1]]
            sage: vec = [65,6,7,8,9,10]
            sage: unit = UK.exp(vec); unit
            -253576*z^11 + 7003*z^10 - 395532*z^9 - 35275*z^8 - 500326*z^7 - 35275*z^6 - 395532*z^5 + 7003*z^4 - 253576*z^3 - 59925*z - 59925
            sage: UK.log(unit)
            [13, 6, 7, 8, 9, 10]
        """
        return self(u).list()

    def exp(self, exponents):
        """
        Return unit with given exponents with respect to group generators.

        INPUT:
            u -- Any object from which an element of the unit group's
            number field K may be constructed; an error is raised if
            an element of K cannot be constructed from u, or if the
            element constructed is not a unit.

        OUTPUT: a list of integers giving the exponents of u with
        respect to the unit group's basis.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: K.<z> = CyclotomicField(13)
            sage: UK = UnitGroup(K)
            sage: [UK.log(u) for u in UK.gens()]
            [[1, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1]]
            sage: vec = [65,6,7,8,9,10]
            sage: unit = UK.exp(vec)
            sage: UK.log(unit)
            [13, 6, 7, 8, 9, 10]
            sage: UK.exp(UK.log(u)) == u
            True
        """
        return prod([u**e for u,e in zip(self.gens(),exponents)], self.number_field().one_element())


