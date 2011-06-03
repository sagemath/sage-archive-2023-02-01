# -*- coding: utf-8 -*-
r"""
Class Groups of Number Fields

An element of a class group is stored as a pair consisting of both an explicit
ideal in that ideal class, and a list of exponents giving that ideal class in
terms of the generators of the parent class group. These can be accessed with
the ``ideal()`` and ``list()`` methods respectively.

EXAMPLES::

    sage: K.<a> = NumberField(x^2 + 23)
    sage: I = K.class_group().gen(); I
    Fractional ideal class (2, 1/2*a - 1/2)
    sage: J = I * I; J
    Fractional ideal class (2, 1/2*a + 1/2)
    sage: J.list()
    [2]
    sage: O = K.OK(); O
    Maximal Order in Number Field in a with defining polynomial x^2 + 23
    sage: O*(2, 1/2*a + 1/2)
    Fractional ideal (2, 1/2*a + 1/2)
    sage: (O*(2, 1/2*a + 1/2)).is_principal()
    False
    sage: (O*(2, 1/2*a + 1/2))^3
    Fractional ideal (1/2*a - 3/2)
"""

from sage.groups.abelian_gps.abelian_group import AbelianGroup_class
from sage.structure.sequence import Sequence
from sage.structure.element import MonoidElement
from sage.groups.abelian_gps.abelian_group_element import AbelianGroupElement
from sage.groups.group import Group
from sage.rings.arith import LCM

class ClassGroup(AbelianGroup_class):
    r"""
    The class group of a number field.
    """
    def __init__(self, invariants, names, number_field, gens, proof=True):
        r"""
        Create a class group.

        Note that the error in the test suite below is caused by the fact that
        there is no category of additive abelian groups.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 23)
            sage: G = K.class_group(); G
            Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^2 + 23

            sage: G.category()
            Category of groups
            sage: TestSuite(G).run() # see #7945
              Failure in _test_category:
            ...
            The following tests failed: _test_elements
        """
        AbelianGroup_class.__init__(self, len(invariants), invariants, names)
        self._proof_flag = proof
        self.__number_field = number_field
        self.__gens = Sequence([FractionalIdealClass(self, x) for x in gens], immutable=True,
                               universe=self, check=False)

    def __call__(self, *args, **kwds):
        r"""
        Call method. This exists *purely* to override the old-fashioned
        behaviour of the parent AbelianGroup class and ensure that
        :meth:`element_constructor` gets called.

        EXAMPLE::

            sage: K.<b> = NumberField(x^2 + 389)
            sage: C = K.class_group()
            sage: C(K.ideal(b))
            Trivial principal fractional ideal class
        """
        return Group.__call__(self, *args, **kwds)

    def _element_constructor_(self, *args, **kwds):
        r"""
        Create an element of this class group from the given data. This may be:
        an ideal class in this number field; an ideal class in a subfield; or
        anything from which an ideal in this number field can be constructed.

        EXAMPLES::

            sage: K.<b> = NumberField(x^2 + 389)
            sage: C = K.class_group()
            sage: C(K.ideal(b)) # indirect doctest
            Trivial principal fractional ideal class
            sage: C(K.ideal(59049, b + 35312)) # indirect doctest
            Fractional ideal class (59049, b + 35312)
            sage: C((59049, b + 35312)) # indirect doctest
            Fractional ideal class (59049, b + 35312)
            sage: C(59049, b + 35312) # indirect doctest
            Fractional ideal class (59049, b + 35312)

            sage: K.<a> = QuadraticField(-23)
            sage: L.<b> = K.extension(x^2 - 2)
            sage: CK = K.class_group()
            sage: CL = L.class_group()
            sage: [CL(I).list() for I in CK]
            [[0], [2], [4]]
        """
        if isinstance(args[0], FractionalIdealClass):
            return FractionalIdealClass(self, self.__number_field.ideal(args[0].ideal()))
        else:
            I = self.__number_field.ideal(*args, **kwds)
            if I.is_zero(): raise TypeError, "The zero ideal is not a fractional ideal"
            return FractionalIdealClass(self, I)

    def gens(self):
        r"""
        Return generators for the class group.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 + 23)
            sage: K.class_group().gens()   # random gens (platform dependent)
            [Fractional ideal class (2, 1/2*a^2 - a + 3/2)]
        """
        return self.__gens

    def ngens(self):
        r"""
        Return the number of generators of the class group.

        EXAMPLES::

            sage: C = NumberField(x^2 + x + 23899, 'a').class_group(); C
            Class group of order 68 with structure C34 x C2 of Number Field in a with defining polynomial x^2 + x + 23899
            sage: C.ngens()
            2
        """
        return len(self.invariants())

    def gen(self, i=0):
        r"""
        Return the i-th generator for this class group.

        EXAMPLES::

            sage: C = NumberField(x^2 + 120071, 'a').class_group(); C
            Class group of order 500 with structure C250 x C2 of Number Field in a with defining polynomial x^2 + 120071
            sage: C.gen(0) # random
            Fractional ideal class (130, 1/2*a + 137/2)
            sage: C.gen(1) # random
            Fractional ideal class (7, a)
        """
        if i < 0 or i >= len(self.__gens):
            raise IndexError
        return self.__gens[i]

    def __iter__(self):
        r"""
        Return an iterator of all ideal classes in this class group.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 + 23)
            sage: G = K.class_group()
            sage: G
            Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^4 + 23
            sage: list(G)
            [Trivial principal fractional ideal class, Fractional ideal class (2, 1/4*a^3 - 1/4*a^2 + 1/4*a - 1/4), Fractional ideal class (2, 1/2*a^2 + 1/2)]
            sage: G.list()
            [Trivial principal fractional ideal class, Fractional ideal class (2, 1/4*a^3 - 1/4*a^2 + 1/4*a - 1/4), Fractional ideal class (2, 1/2*a^2 + 1/2)]

        TESTS::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: G = K.class_group()
            sage: G
            Class group of order 1 of Number Field in a with defining polynomial x^2 + 1
            sage: list(G)
            [Trivial principal fractional ideal class]
            sage: G.list()
            [Trivial principal fractional ideal class]
        """
        from sage.misc.mrange import mrange
        invs = self.invariants()
        T = mrange(invs)
        g = self.gens()
        for t in T:
            I = self(1)
            for i, j in enumerate(t):
                I *= g[i]**j
            yield I
        if not T:
            yield self(1)

    def _repr_(self):
        r"""
        Return string representation of self.

        EXAMPLES::

            sage: C = NumberField(x^2 + 23, 'a').class_group()
            sage: C._repr_()
            'Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^2 + 23'
        """
        s = 'Class group of order %s '%self.order()
        if self.order() > 1:
            s += 'with structure %s '%self._group_notation(self.invariants())
        s += 'of %s'%self.number_field()
        return s

    def number_field(self):
        r"""
        Return the number field that this class group is attached to.

        EXAMPLES::

            sage: C = NumberField(x^2 + 23, 'w').class_group(); C
            Class group of order 3 with structure C3 of Number Field in w with defining polynomial x^2 + 23
            sage: C.number_field()
            Number Field in w with defining polynomial x^2 + 23
        """
        return self.__number_field


class FractionalIdealClass(AbelianGroupElement):
    r"""
    A fractional ideal class in a number field.

    EXAMPLES::

        sage: G = NumberField(x^2 + 23,'a').class_group(); G
        Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^2 + 23
        sage: I = G.0; I
        Fractional ideal class (2, 1/2*a - 1/2)
    """
    def __init__(self, parent, ideal, element=None):
        """
        Returns the ideal class of this fractional ideal.

        EXAMPLE::

            sage: K.<a> = NumberField(x^2 + 23,'a'); G = K.class_group()
            sage: G(K.ideal(13, a + 4))
            Fractional ideal class (13, 1/2*a + 17/2)
        """
        self.__ideal = ideal
        if element is None:
            element = map(int, ideal._ideal_class_log(proof=parent._proof_flag))
        AbelianGroupElement.__init__(self, parent, element)

    def _repr_(self):
        r"""
        Return string representation of this fractional ideal class.

         EXAMPLE::

            sage: K.<a> = NumberField(x^2 + 23,'a'); G = K.class_group()
            sage: G(K.ideal(13, a + 4))._repr_()
            'Fractional ideal class (13, 1/2*a + 17/2)'
            sage: G(K.ideal(59, a+6))._repr_()
            'Trivial principal fractional ideal class'
        """
        if self.is_principal():
            return 'Trivial principal fractional ideal class'
        return 'Fractional ideal class %s'%self.__ideal._repr_short()

    def _mul_(self, other):
        r"""
        Multiplication of two ideal classes.

        EXAMPLE::

            sage: G = NumberField(x^2 + 23,'a').class_group(); G
            Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^2 + 23
            sage: I = G.0; I
            Fractional ideal class (2, 1/2*a - 1/2)
            sage: I*I # indirect doctest
            Fractional ideal class (2, 1/2*a + 1/2)
            sage: I*I*I # indirect doctest
            Trivial principal fractional ideal class
        """
        m = AbelianGroupElement._mul_(self, other)
        return FractionalIdealClass(self.parent(), (self.__ideal * other.__ideal).reduce_equiv(), m.list())

    def __pow__(self, n):
        r"""
        Raise this element to the power n.

        EXAMPLE::

            sage: K.<a> = NumberField(x^3 - 3*x + 8)
            sage: C=K.class_group()
            sage: c = C(2, a)
            sage: c^2
            Fractional ideal class (2, a^2 + 2*a - 1)
            sage: c^3
            Trivial principal fractional ideal class
            sage: c^1000
            Fractional ideal class (2, a)
            sage: (c^2)^2
            Fractional ideal class (2, a)
        """
        # We use MonoidElement's __pow__ routine, since that does
        # repeated squaring, and hence the ideal gets reduced as
        # we go along; actually computing self.__ideal ** n would
        # be disastrous.
        n = n % self.order()
        return MonoidElement.__pow__(self, n)

    def inverse(self):
        r"""
        Return the multiplicative inverse of this ideal class.

        EXAMPLE::

            sage: K.<a> = NumberField(x^3 - 3*x + 8); G = K.class_group()
            sage: G(2, a).inverse()
            Fractional ideal class (2, a^2 + 2*a - 1)
        """
        m = AbelianGroupElement.inverse(self)
        return FractionalIdealClass(self.parent(), (~self.__ideal).reduce_equiv(), m.list())

    def __invert__(self):
        r"""
        Return the multiplicative inverse of this ideal class.

        EXAMPLE::

            sage: K.<a> = NumberField(x^3 - 3*x + 8); G = K.class_group()
            sage: ~G(2, a)
            Fractional ideal class (2, a^2 + 2*a - 1)
        """
        return self.inverse()

    def is_principal(self):
        r"""
        Returns True iff this ideal class is the trivial (principal) class

        EXAMPLES::

            sage: K.<w>=QuadraticField(-23)
            sage: OK=K.ring_of_integers()
            sage: C=OK.class_group()
            sage: P2a,P2b=[P for P,e in (2*OK).factor()]
            sage: c=C(P2a)
            sage: c.is_principal()
            False
            sage: (c^2).is_principal()
            False
            sage: (c^3).is_principal()
            True
        """
        return self.is_one()

    def reduce(self):
        r"""
        Return representative for this ideal class that has been
        reduced using PARI's idealred.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 20072); G = k.class_group(); G
            Class group of order 76 with structure C38 x C2 of Number Field in a with defining polynomial x^2 + 20072
            sage: I = G.0; I
            Fractional ideal class (41, 1/2*a + 5)
            sage: J = G(I.ideal()^5); J
            Fractional ideal class (115856201, 1/2*a + 40407883)
            sage: J.reduce()
            Fractional ideal class (57, 1/2*a + 44)
        """
        return self.parent()(self.__ideal.reduce_equiv())

    def order(self):
        r"""
        Return the order of this ideal class in the class group.

        EXAMPLE::

            sage: K.<w>=QuadraticField(-23)
            sage: OK=K.ring_of_integers()
            sage: C=OK.class_group()
            sage: [c.order() for c in C]
            [1, 3, 3]

            sage: k.<a> = NumberField(x^2 + 20072); G = k.class_group(); G
            Class group of order 76 with structure C38 x C2 of Number Field in a with defining polynomial x^2 + 20072
            sage: [c.order() for c in G.gens()]
            [38, 2]

        """
        # an old method with a new docstring
        return AbelianGroupElement.order(self)

    def multiplicative_order(self):
        r"""
        Alias for :meth:`order`.

        EXAMPLE::

            sage: K.<w>=QuadraticField(-23)
            sage: K.class_group()(K.primes_above(2)[0]).multiplicative_order()
            3
        """
        return self.order()

    def ideal(self):
        r"""
        Return a representative ideal in this ideal class.

        EXAMPLE::

            sage: K.<w>=QuadraticField(-23)
            sage: OK=K.ring_of_integers()
            sage: C=OK.class_group()
            sage: P2a,P2b=[P for P,e in (2*OK).factor()]
            sage: c=C(P2a); c
            Fractional ideal class (2, 1/2*w - 1/2)
            sage: c.ideal()
            Fractional ideal (2, 1/2*w - 1/2)
        """
        return self.__ideal

    def gens(self):
        r"""
        Return generators for a representative ideal in this
        ideal class.

        EXAMPLES::

            sage: K.<w>=QuadraticField(-23)
            sage: OK=K.ring_of_integers()
            sage: C=OK.class_group()
            sage: P2a,P2b=[P for P,e in (2*OK).factor()]
            sage: c=C(P2a); c
            Fractional ideal class (2, 1/2*w - 1/2)
            sage: c.gens()
            (2, 1/2*w - 1/2)
        """
        return self.ideal().gens()


class SClassGroup(ClassGroup):
    r"""
    The S-class group of a number field.
    """

    def __init__(self, invariants, names, number_field, gens, S, proof=True):
        r"""
        Create an S-class group.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: K.S_class_group(S)
            S-class group of order 2 with structure C2 of Number Field in a with defining polynomial x^2 + 14
            sage: K.<a> = QuadraticField(-105)
            sage: K.S_class_group([K.ideal(13, a + 8)])
            S-class group of order 4 with structure C2 x C2 of Number Field in a with defining polynomial x^2 + 105
        """
        AbelianGroup_class.__init__(self, len(invariants), invariants, names)
        self._proof_flag = proof
        self.__number_field = number_field
        self.__S = S
        self.__gens = Sequence([SFractionalIdealClass(self, x) for x in gens], immutable=True,
                               universe=self, check=False)

    def S(self):
        r"""
        Return the set (or rather tuple) of primes used to define this class group.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S);CS
            S-class group of order 2 with structure C2 of Number Field in a with defining polynomial x^2 + 14
            sage: T = tuple([])
            sage: CT = K.S_class_group(T);CT
            S-class group of order 4 with structure C4 of Number Field in a with defining polynomial x^2 + 14
            sage: CS.S()
            (Fractional ideal (2, a),)
            sage: CT.S()
            ()
        """
        return self.__S

    def __call__(self, *args, **kwds):
        r"""
        Call method.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: J = K.ideal(7,a)
            sage: G = K.ideal(3,a+1)
            sage: CS(I)
            Trivial S-ideal class
            sage: CS(J)
            Trivial S-ideal class
            sage: CS(G)
            Fractional S-ideal class (3, a + 1)
        """
        return Group.__call__(self, *args, **kwds)

    def _element_constructor_(self, *args, **kwds):
        r"""
        Create an element of this class group from the given data.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: J = K.ideal(7,a)
            sage: G = K.ideal(3,a+1)
            sage: CS(I)
            Trivial S-ideal class
            sage: CS(J)
            Trivial S-ideal class
            sage: CS(G)
            Fractional S-ideal class (3, a + 1)
        """
        if isinstance(args[0], FractionalIdealClass):
            return SFractionalIdealClass(self, args[0].ideal())
        else:
            I = self.number_field().ideal(*args, **kwds)
            if I.is_zero(): raise TypeError, "The zero ideal is not a fractional ideal"
            return SFractionalIdealClass(self, I)

    def gens(self):
        r"""
        Return generators for the S-class group.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: S = K.primes_above(2)
            sage: K.S_class_group(S).gens()   # random gens (platform dependent)
            [Fractional S-ideal class (3, a + 2)]
        """
        return self.__gens

    def gen(self, i=0):
        r"""
        Return the i-th generator for this class group.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-974)
            sage: CS = K.S_class_group(K.primes_above(2)); CS
            S-class group of order 18 with structure C6 x C3 of Number Field in a with defining polynomial x^2 + 974
            sage: CS.gen(0) # random
            Fractional S-ideal class (3, a + 2)
            sage: CS.gen(1) # random
            Fractional S-ideal class (31, a + 24)
        """
        if i < 0 or i >= len(self.__gens):
            raise IndexError
        return self.__gens[i]

    def _repr_(self):
        r"""
        Return string representation of this S-class group.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: CS = K.S_class_group(K.primes_above(2))
            sage: CS._repr_()
            'S-class group of order 2 with structure C2 of Number Field in a with defining polynomial x^2 + 14'
        """
        s = 'S-class group of order %s ' % self.order()
        if self.order() > 1:
            s += 'with structure %s ' % self._group_notation(self.invariants())
        s += 'of %s' % self.number_field()
        return s

    def number_field(self):
        r"""
        Return the number field that this S-class group is attached to.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: CS = K.S_class_group(K.primes_above(2))
            sage: CS.number_field()
            Number Field in a with defining polynomial x^2 + 14
        """
        return self.__number_field


class SFractionalIdealClass(FractionalIdealClass):
    r"""
    An S-fractional ideal class in a number field for a tuple of primes S.
    """

    def __init__(self, parent, ideal, element=None):
        r"""
        Returns the S-ideal class of this fractional ideal.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: J = K.ideal(7,a)
            sage: G = K.ideal(3,a+1)
            sage: CS(I)
            Trivial S-ideal class
            sage: CS(J)
            Trivial S-ideal class
            sage: CS(G)
            Fractional S-ideal class (3, a + 1)
        """
        self.__ideal = ideal
        if element is None:
            element = ideal._S_ideal_class_log(parent.S())
        AbelianGroupElement.__init__(self, parent, element)

    def _repr_(self):
        r"""
        Returns a string representation of the S-ideal class of this fractional ideal.

        EXAMPLE::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: J = K.ideal(3, a + 2)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: CS(J)
            Fractional S-ideal class (3, a + 2)
            sage: CS(J^2)
            Trivial S-ideal class
        """
        if not any(self.list()):
            return 'Trivial S-ideal class'
        return 'Fractional S-ideal class %s' % self.__ideal._repr_short()

    def ideal(self):
        r"""
        Returns a representative ideal for this S-ideal class.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: J = K.ideal(7,a)
            sage: G = K.ideal(3,a+1)
            sage: CS(I).ideal()
            Fractional ideal (2, a)
            sage: CS(J).ideal()
            Fractional ideal (7, a)
            sage: CS(G).ideal()
            Fractional ideal (3, a + 1)
        """
        return self.__ideal

    def order(self):
        r"""
        Finds the order of the given S-ideal class.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: J = K.ideal(7,a)
            sage: G = K.ideal(3,a+1)
            sage: CS(I).order()
            1
            sage: CS(J).order()
            1
            sage: CS(G).order()
            2
        """
        return LCM([e.additive_order() for e in self.list()])

    def _mul_(self, other):
        r"""
        Multiplies together two S-ideal classes.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: G = K.ideal(3,a+1)
            sage: CS(G)*CS(G)
            Trivial S-ideal class
        """

        m = AbelianGroupElement._mul_(self, other)
        return SFractionalIdealClass(self.parent(), (self.ideal() * other.ideal()).reduce_equiv(), m.list())

    def inverse(self):
        r"""
        Finds the inverse of the given S-ideal class.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: G = K.ideal(3,a+1)
            sage: CS(G).inverse()
            Fractional S-ideal class (3, a + 2)
        """
        m = AbelianGroupElement.inverse(self)
        inv_ideal = self.parent().number_field().class_group()(self.ideal()).inverse().ideal()
        return SFractionalIdealClass(self.parent(), inv_ideal , m.list())
