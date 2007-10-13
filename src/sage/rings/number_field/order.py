"""
Orders in number fields.

AUTHORS:
    -- William Stein and Robert Bradshaw (2007-09): initial version
"""

from sage.rings.ring import IntegralDomain, DedekindDomain
from sage.structure.sequence import Sequence
from sage.rings.integer_ring import ZZ
from sage.structure.element import is_Element

from number_field_element import OrderElement_absolute, OrderElement_relative
from sage.rings.monomials import monomials


def is_NumberFieldOrder(R):
    """
    Return True if R an order in a number field or R is the ring ZZ of integers.

    EXAMPLES:
        sage: is_NumberFieldOrder(NumberField(x^2+1,'a').maximal_order())
        True
        sage: is_NumberFieldOrder(ZZ)
        True
        sage: is_NumberFieldOrder(QQ)
        False
        sage: is_NumberFieldOrder(45)
        False
    """
    return isinstance(R, Order) or R == ZZ

class Order(IntegralDomain):
    r"""
    An order in a number field.

    An order is a subring of the number field that has $\ZZ$-rank equal
    to the degree of the number field over $\QQ$.

    EXAMPLES:
        sage: K.<theta> = NumberField(x^4 + x + 17)
        sage: K.maximal_order()
        Order with module basis 1, theta, theta^2, theta^3 in Number Field in theta with defining polynomial x^4 + x + 17
        sage: K.order(17*theta)
        Order with module basis 1, 17*theta, 289*theta^2, 4913*theta^3 in Number Field in theta with defining polynomial x^4 + x + 17
        sage: K.order(17*theta, 13*theta)
        Order with module basis 1, theta, theta^2, theta^3 in Number Field in theta with defining polynomial x^4 + x + 17
        sage: K.order([34*theta, 17*theta + 17])
        Order with module basis 1, 17*theta, 289*theta^2, 4913*theta^3 in Number Field in theta with defining polynomial x^4 + x + 17
        sage: K.<b> = NumberField(x^4 + x^2 + 2)
        sage: (b^2).charpoly().factor()
        (x^2 + x + 2)^2
        sage: K.order(b^2)
        Traceback (most recent call last):
        ...
        ValueError: the rank of the span of gens is wrong
    """
    def __init__(self, K, is_maximal):
        """
        This is called when creating an order to set the ambient field.

        EXAMPLES:
            sage: k = CyclotomicField(5)
            sage: k.maximal_order()
            Order with module basis 1, zeta5, zeta5^2, zeta5^3 in Cyclotomic Field of order 5 and degree 4
        """
        self._K = K
        self._is_maximal = is_maximal
        DedekindDomain.__init__(self, base = K.base(), names = K.variable_names(), normalize = False) # base should probably change

    def __mul__(self, right):
        """
        Create an ideal in this order using the notation Ok*gens

        EXAMPLES:
            sage: k.<a> = NumberField(x^2 + 5077); G = k.class_group(); G
            Class group of order 22 with structure C22 of Number Field in a with defining polynomial x^2 + 5077
            sage: G.0   # random output
            Fractional ideal class (11, a - 4) of Number Field in a with defining polynomial x^2 + 5077
            sage: Ok = k.maximal_order(); Ok
            Order with module basis 1, a in Number Field in a with defining polynomial x^2 + 5077
            sage: Ok*(11, a - 4)
            Fractional ideal (11, a - 4) of Number Field in a with defining polynomial x^2 + 5077
            sage: (11, a - 4) * Ok
            Fractional ideal (11, a - 4) of Number Field in a with defining polynomial x^2 + 5077
        """
        if self.is_maximal():
            return self._K.ideal(right)
        raise TypeError

    def __rmul__(self, left):
        """
        Create an ideal in this order using the notation gens*Ok.

        EXAMPLES:
            sage: k.<a> = NumberField(x^2 + 431); G = k.class_group(); G
            Class group of order 21 with structure C21 of Number Field in a with defining polynomial x^2 + 431
            sage: G.0   # random output
            Fractional ideal class (6, 1/2*a + 11/2) of Number Field in a with defining polynomial x^2 + 431
            sage: Ok = k.maximal_order(); Ok
            Order with module basis 1/2*a + 1/2, a in Number Field in a with defining polynomial x^2 + 431
            sage: (6, 1/2*a + 11/2)*Ok    # random output
            Fractional ideal (6, 1/2*a + 11/2) of Number Field in a with defining polynomial x^2 + 431
            sage: 17*Ok
            Principal ideal (17) of Order with module basis 1/2*a + 1/2, a in Number Field in a with defining polynomial x^2 + 431
        """
        return self.__mul__(left)

    def is_maximal(self):
        """
        Returns True if this is the maximal order.

            sage: k.<i> = NumberField(x^2 + 1)
            sage: O3 = k.order(3*i); O5 = k.order(5*i); Ok = k.maximal_order(); Osum = O3 + O5
            sage: Osum.is_maximal()
            True
            sage: O3.is_maximal()
            False
            sage: O5.is_maximal()
            False
            sage: Ok.is_maximal()
            True
        """
        if self._is_maximal is None:
            self._is_maximal = (self.discriminant() == self._K.discriminant())
        return self._is_maximal

    def is_integrally_closed(self):
        """
        Return True if this ring is integrally closed, i.e., is equal
        to the maximal order.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + 189*x + 394)
            sage: R = K.order(2*a)
            sage: R.is_integrally_closed()
            False
            sage: R
            Order with module basis 1, 2*a in Number Field in a with defining polynomial x^2 + 189*x + 394
            sage: S = K.maximal_order(); S
            Order with module basis 1, a in Number Field in a with defining polynomial x^2 + 189*x + 394
            sage: S.is_integrally_closed()
            True
        """
        return self.is_maximal()


    def integral_closure(self):
        """
        Return the integral closure of this order.

        EXAMPLES:
            sage: K.<a> = QuadraticField(5)
            sage: O2 = K.order(2*a); O2
            Order with module basis 1, 2*a in Number Field in a with defining polynomial x^2 - 5
            sage: O2.integral_closure()
            Order with module basis 1/2*a + 1/2, a in Number Field in a with defining polynomial x^2 - 5
            sage: OK = K.maximal_order()
            sage: OK is OK.integral_closure()
            True
        """
        if self.is_maximal():
            return self
        else:
            return self.number_field().maximal_order()

    def gen(self, i):
        """
        Return i-th module generator of this order.

        EXAMPLES:
            sage: K.<c> = NumberField(x^3 + 2*x + 17)
            sage: O = K.maximal_order(); O
            Order with module basis 1, c, c^2 in Number Field in c with defining polynomial x^3 + 2*x + 17
            sage: O.gen(1)
            c
            sage: O.gen(2)
            c^2
            sage: O.gen(5)
            Traceback (most recent call last):
            ...
            IndexError: no 5th generator
            sage: O.gen(-1)
            Traceback (most recent call last):
            ...
            IndexError: no -1th generator
        """
        b = self.basis()
        if i < 0 or i >= len(b):
            raise IndexError, "no %sth generator"%i
        return self.basis()[i]

    def gens(self):
        """
        Return a list of the module generators of this order.

        NOTE: For a (much smaller) list of ring generators use
        \code{self.ring_generators()}.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
            sage: O = K.maximal_order()
            sage: O.gens()
            [1, 1/2*a^2 + 1/2*a, a^2]
        """
        return self.basis()

    def ngens(self):
        """
        Return the number of module generators of this order.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
            sage: O = K.maximal_order()
            sage: O.ngens()
            3
        """
        return self.absolute_degree()

    def basis(self):  # this must be defined in derived class
        """
        Return a basis over ZZ of this order.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + x^2 - 16*x + 16)
            sage: O = K.maximal_order(); O
            Order with module basis 1, 1/4*a^2 + 1/4*a, a^2 in Number Field in a with defining polynomial x^3 + x^2 - 16*x + 16
            sage: O.basis()
            [1, 1/4*a^2 + 1/4*a, a^2]
        """
        raise NotImplementedError

    def ring_generators(self):
        """
        Return generators for self as a ring.

        EXAMPLES:
            sage: K.<i> = NumberField(x^2 + 1)
            sage: O = K.maximal_order(); O
            Order with module basis 1, i in Number Field in i with defining polynomial x^2 + 1
            sage: O.ring_generators()
            [i]

        This is an example where 2 generators are required (because 2 is an essential
        discriminant divisor).
            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
            sage: O = K.maximal_order(); O.basis()
            [1, 1/2*a^2 + 1/2*a, a^2]
            sage: O.ring_generators()
            [1/2*a^2 + 1/2*a, a^2]
        """
        try:
            return self.__ring_generators
        except AttributeError:
            K = self._K
            n = K.degree()
            V, from_V, to_V = self._K.vector_space()
            A = ZZ**self.rank()
            remaining = [x for x in self.basis() if x != 1]
            gens = []
            while len(remaining) > 0:
                gens.append(remaining[0])
                del remaining[0]
                W = A.span([to_V(x) for x in monomials(gens, n)])
                remaining = [x for x in remaining if not to_V(x) in W]
            self.__ring_generators = Sequence(gens,immutable=True)
            return self.__ring_generators



    def number_field(self):
        """
        Return the number field of this order, which is the ambient
        number field that this order is embedded in.

        EXAMPLES:
            sage: K.<b> = NumberField(x^4 + x^2 + 2)
            sage: O = K.order(2*b); O
            Order with module basis 1, 2*b, 4*b^2, 8*b^3 in Number Field in b with defining polynomial x^4 + x^2 + 2
            sage: O.number_field()
            Number Field in b with defining polynomial x^4 + x^2 + 2
            sage: O.number_field() is K
            True
        """
        return self._K

    def ambient(self):
        r"""
        Return the ambient number field that contains self.

        This is the same as \code{self.number_field()} and
        \code{self.fraction_field()}

        EXAMPLES:
            sage: k.<z> = NumberField(x^2 - 389)
            sage: o = k.order(389*z + 1)
            sage: o
            Order with module basis 1, 389*z in Number Field in z with defining polynomial x^2 - 389
            sage: o.ambient()
            Number Field in z with defining polynomial x^2 - 389
        """
        return self._K

    def residue_field(self, prime, name = None, check = False):
        """
        Return the residue field of this number field at a given prime, ie $O_K / p O_K$.

        INPUT:
            prime -- a prime ideal of the maximal order in this number field.
            name -- the name of the variable in the residue field
            check -- whether or not to check the primality of prime.
        OUTPUT:
            The residue field at this prime.

        EXAMPLES:
        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(x^4+3*x^2-17)
        sage: P = K.ideal(61).factor()[0][0]
        sage: OK = K.maximal_order()
        sage: OK.residue_field(P)
        Residue field of Fractional ideal (-2*a^2 + 1) of Number Field in a with defining polynomial x^4 + 3*x^2 - 17
        """
        import sage.rings.residue_field
        return sage.rings.residue_field.ResidueField(prime)


    def fraction_field(self):
        """
        Return the fraction field of this order, which is the
        ambient number field.

        EXAMPLES:
        sage: K.<b> = NumberField(x^4 + 17*x^2 + 17)
        sage: O = K.order(17*b); O
        Order with module basis 1, 17*b, 289*b^2, 4913*b^3 in Number Field in b with defining polynomial x^4 + 17*x^2 + 17
        sage: O.fraction_field()
        Number Field in b with defining polynomial x^4 + 17*x^2 + 17
        """
        return self._K

    def degree(self):
        r"""
        Return the degree of this order, which is the rank
        of this order as a $\ZZ$-module.

        EXAMPLES:
        sage: k.<c> = NumberField(x^3 + x^2 - 2*x+8)
        sage: o = k.maximal_order()
        sage: o.degree()
        3
        sage: o.rank()
        3
        """
        return self._K.degree()

    def rank(self):
        r"""
        Return the rank of this order, which is the rank of
        the underlying $\ZZ$-module, or the degree of the ambient
        number field that contains this order.

        This is a synonym for \code{self.degree()}.

        EXAMPLES:
            sage: k.<c> = NumberField(x^5 + x^2 + 1)
            sage: o = k.maximal_order(); o
            Order with module basis 1, c, c^2, c^3, c^4 in Number Field in c with defining polynomial x^5 + x^2 + 1
            sage: o.rank()
            5
        """
        return self.degree()

    def class_group(self, proof=None, names='c'):
        r"""
        Return the class group of this order.

        (Currently only implemented for the maximal order.)

        EXAMPLES:
            sage: k.<a> = NumberField(x^2 + 5077)
            sage: O = k.maximal_order(); O
            Order with module basis 1, a in Number Field in a with defining polynomial x^2 + 5077
            sage: O.class_group()
            Class group of order 22 with structure C22 of Number Field in a with defining polynomial x^2 + 5077
        """
        if self.is_maximal():
            return self.number_field().class_group(proof=proof, names=names)
        else:
            raise NotImplementedError

    def is_suborder(self, other):
        """
        Return True if self and other are both orders in the
        same ambient number field and self is a subset of other.

        EXAMPLES:
            sage: W.<i> = NumberField(x^2 + 1)
            sage: O5 = W.order(5*i)
            sage: O10 = W.order(10*i)
            sage: O15 = W.order(15*i)
            sage: O15.is_suborder(O5)
            True
            sage: O5.is_suborder(O15)
            False
            sage: O10.is_suborder(O15)
            False

        We create another isomorphic but different field:
            sage: W2.<j> = NumberField(x^2 + 1)
            sage: P5 = W2.order(5*j)

        This is False because the ambient number fields are not equal.
            sage: O5.is_suborder(P5)
            False

        We create a field that contains (in no natural way!) W,
        and of course again is_suborder returns False:
            sage: K.<z> = NumberField(x^4 + 1)
            sage: M = K.order(5*z)
            sage: O5.is_suborder(M)
            False
        """
        if not isinstance(other, Order):
            return False
        if other.number_field() != self.number_field():
            return False
        return self.module().is_submodule(other.module())

    def __cmp__(self, other):
        r"""
        Compare the order self to other.

        NOTE: This is a well defined way to compare any two objects,
        but it is not the partial inclusion ordering!.  Thus self <
        other being True does not necessarily mean that self is
        contained in other.  Use \code{self.is_suborder(other)} to
        determine inclusion.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + 2)
            sage: O1 = K.order(a); O1
            Order with module basis 1, a, a^2 in Number Field in a with defining polynomial x^3 + 2
            sage: O2 = K.order(a^2); O2
            Order with module basis 1, 2*a, a^2 in Number Field in a with defining polynomial x^3 + 2
            sage: O1 == O2
            False
            sage: O1 < O2
            True
            sage: O2 < O1
            False

        Note that "less than" does not mean "is a subset":
            sage: O2.is_suborder(O1)
            True
            sage: O1 == K
            False
            sage: K == O1
            False
        """
        if not isinstance(other, Order):
            return cmp(type(self), type(other))
        if self._K != other._K:
            return cmp(self._K, other._K)
        return cmp(self._module_rep, other._module_rep)

    def absolute_degree(self):
        """
        Returns the absolute degree of this order, ie the degree of this order over ZZ.

        EXAMPLES:
        sage: K.<a> = NumberField(x^3 + 2)
        sage: O = K.maximal_order()
        sage: O.absolute_degree()
        3
        sage: K.<a> = NumberField([x^3 + 2, x^2 - 3])
        sage: O = K.maximal_order()
        Traceback (most recent call last):
        ...
        NotImplementedError
        """
        return self.number_field().absolute_degree()

##     def absolute_polynomial(self):
##         """
##         Returns the absolute polynomial of this order, which is just the absolute polynomial of the number field.

##         EXAMPLES:
##         sage: K.<a, b> = NumberField([x^2 + 1, x^3 + x + 1]); OK = K.maximal_order()
##         Traceback (most recent call last):
##         ...
##         NotImplementedError

##         #sage: OK.absolute_polynomial()
##         #x^6 + 5*x^4 - 2*x^3 + 4*x^2 + 4*x + 1
##         """
##         return self.number_field().absolute_polynomial()

##     def polynomial(self):
##         """
##         Returns the polynomial defining the number field that contains self.
##         """
##         return self.number_field().polynomial()

##     def polynomial_ntl(self):
##         """
##         Return defining polynomial of the parent number field as a
##         pair, an ntl polynomial and a denominator.

##         This is used mainly to implement some internal arithmetic.

##         EXAMPLES:
##             sage: NumberField(x^2 + 1,'a').maximal_order().polynomial_ntl()
##             ([1 0 1], 1)
##         """
##         return self.number_field().polynomial_ntl()

class AbsoluteOrder(Order):

    def __init__(self, K, module_rep, is_maximal=None, check=True):
        """
        EXAMPLES:
            sage: from sage.rings.number_field.order import *
            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^3+2)
            sage: V, from_v, to_v = K.vector_space()
            sage: M = span(ZZ, [to_v(a^2), to_v(a), to_v(1)])
            sage: O = AbsoluteOrder(K, M); O
            Order with module basis 1, a, a^2 in Number Field in a with defining polynomial x^3 + 2

            sage: M = span(ZZ, [to_v(a^2), to_v(a), to_v(2)])
            sage: O = AbsoluteOrder(K, M); O
            Traceback (most recent call last):
            ...
            ValueError: 1 is not in the span of the module, hence not an order.

            sage: loads(dumps(O)) == O
            True
        """
        Order.__init__(self, K, is_maximal=is_maximal)
        self._module_rep = module_rep
        V, from_v, to_v = self._K.vector_space()
        if check:
            if not K.is_absolute():
                raise ValueError, "AbsoluteOrder must be called with an absolute number field."
            if to_v(1) not in module_rep:
                raise ValueError, "1 is not in the span of the module, hence not an order."
            if module_rep.rank() != self._K.degree():
                raise ValueError, "the module must have full rank."

    def __call__(self, x):
        """
        Coerce x into this order.

        EXAMPLES:
            sage: k.<z> = NumberField(x^2 - 389)
            sage: m = k.order(3*z); m
            Order with module basis 1, 3*z in Number Field in z with defining polynomial x^2 - 389
            sage: m(6*z)
            6*z
            sage: k(m(6*z))
            6*z
        """
        if is_Element(x) and x.parent() is self:
            return x
        if not is_Element(x) or x.parent() is not self._K:
            x = self._K(x)
        V, _, embedding = self._K.vector_space()
        if not embedding(x) in self._module_rep:
            raise TypeError, "Not an element of the order."
        return OrderElement_absolute(self, x)

    def __add__(left, right):
        """
        Add two orders.

        EXAMPLES:
            sage: K.<a> = NumberField(polygen(QQ,'z')^3 - 2)
            sage: O6 = K.order(6*a); O6
            Order with module basis 1, 6*a, 36*a^2 in Number Field in a with defining polynomial z^3 - 2
            sage: O15 = K.order(15*a^2); O15
            Order with module basis 1, 450*a, 15*a^2 in Number Field in a with defining polynomial z^3 - 2
            sage: O6 + O15
            Order with module basis 1, 6*a, 3*a^2 in Number Field in a with defining polynomial z^3 - 2
        """
        if not isinstance(left, AbsoluteOrder) or not isinstance(right, AbsoluteOrder):
            raise NotImplementedError
        if left.number_field() != right.number_field():
            raise TypeError, "Number fields don't match."
        if left._is_maximal:
            return left
        elif right._is_maximal:
            return right
        return AbsoluteOrder(left._K, left._module_rep + right._module_rep, None)

    def __and__(left, right):
        """
        Intersect orders.

        EXAMPLES:
            sage: K.<i> = QuadraticField(-1)
            sage: O3 = K.order(3*i); O5 = K.order(5*i)
            sage: O3 & O5
            Order with module basis 1, 15*i in Number Field in i with defining polynomial x^2 + 1
            sage: O3.intersection(O5)
            Order with module basis 1, 15*i in Number Field in i with defining polynomial x^2 + 1
        """
        if not isinstance(left, AbsoluteOrder) or not isinstance(right, AbsoluteOrder):
            raise NotImplementedError
        if left.number_field() != right.number_field():
            raise TypeError, "Number fields don't match."
        return AbsoluteOrder(left._K, left._module_rep.intersection(right._module_rep), False)

    def discriminant(self):
        """
        Return the discriminant of this order.

        EXAMPLES:
            sage: K.<a> = NumberField(x^8 + x^3 - 13*x + 26)
            sage: O = K.maximal_order()
            sage: factor(O.discriminant())
            3 * 11 * 13^2 * 613 * 1575917857
            sage: L = K.order(13*a^2)
            sage: factor(L.discriminant())
            3^3 * 5^2 * 11 * 13^60 * 613 * 733^2 * 1575917857
            sage: factor(L.index_in(O))
            3 * 5 * 13^29 * 733
            sage: L.discriminant() / O.discriminant() == L.index_in(O)^2
            True
        """
        try:
            return self.__discriminant
        except AttributeError:
            if self._is_maximal:
                D = self._K.discriminant()
            else:
                D = self._K.discriminant(self.basis())
            self.__discriminant = D
            return D

    def index_in(self, other):
        """
        Return the index of self in other.  This is a lattice index,
        so it is a rational number if self isn't contained in other.

        INPUT:
            other -- another absolute order with the same ambient
            number field.

        OUTPUT:
            a rational number

        EXAMPLES:
            sage: k.<i> = NumberField(x^2 + 1)
            sage: O1 = k.order(i)
            sage: O5 = k.order(5*i)
            sage: O5.index_in(O1)
            5

            sage: k.<a> = NumberField(x^3 + x^2 - 2*x+8)
            sage: o = k.maximal_order()
            sage: o
            Order with module basis 1, 1/2*a^2 + 1/2*a, a^2 in Number Field in a with defining polynomial x^3 + x^2 - 2*x + 8
            sage: O1 = k.order(a); O1
            Order with module basis 1, a, a^2 in Number Field in a with defining polynomial x^3 + x^2 - 2*x + 8
            sage: O1.index_in(o)
            2
            sage: O2 = k.order(1+2*a); O2
            Order with module basis 1, 2*a, 4*a^2 in Number Field in a with defining polynomial x^3 + x^2 - 2*x + 8
            sage: o.index_in(O2)
            1/16
        """
        if not isinstance(other, AbsoluteOrder):
            raise TypeError, "other must be an absolute order."
        if other.ambient() != self.ambient():
            raise ValueError, "other must have the same ambient number field as self."
        return self._module_rep.index_in(other._module_rep)

    def module(self):
        """
        Returns the underlying free module corresponding to this
        order, embedded in the vector space corresponding to the
        ambient number field.

        EXAMPLES:
            sage: k.<a> = NumberField(x^3 + x + 3)
            sage: m = k.order(3*a); m
            Order with module basis 1, 3*a, 9*a^2 in Number Field in a with defining polynomial x^3 + x + 3
            sage: m.module()
            Free module of degree 3 and rank 3 over Integer Ring
            Echelon basis matrix:
            [1 0 0]
            [0 3 0]
            [0 0 9]
        """
        return self._module_rep

    def intersection(self, other):
        """
        Return the intersection of this order with another order.

        EXAMPLES:
            sage: k.<i> = NumberField(x^2 + 1)
            sage: O6 = k.order(6*i)
            sage: O9 = k.order(9*i)
            sage: O6.intersection(O9)
            Order with module basis 1, 18*i in Number Field in i with defining polynomial x^2 + 1
            sage: O6 & O9
            Order with module basis 1, 18*i in Number Field in i with defining polynomial x^2 + 1
            sage: O6 + O9
            Order with module basis 1, 3*i in Number Field in i with defining polynomial x^2 + 1
        """
        return self & other

    def _repr_(self):
        """
        Return print representation of this absolute order.

        EXAMPLES:
            sage: K.<a> = NumberField(x^4 - 5)
            sage: O = K.maximal_order()
            sage: O._repr_()
            'Order with module basis 1/2*a^2 + 1/2, 1/2*a^3 + 1/2*a, a^2, a^3 in Number Field in a with defining polynomial x^4 - 5'
        """
        return "Order with module basis %s in %r" % (", ".join([str(b) for b in self.basis()]), self._K)

    def basis(self):
        """
        Return the basis over ZZ for this order.

        EXAMPLES:
            sage: k.<c> = NumberField(x^3 + x^2 + 1)
            sage: O = k.maximal_order(); O
            Order with module basis 1, c, c^2 in Number Field in c with defining polynomial x^3 + x^2 + 1
            sage: O.basis()
            [1, c, c^2]

        The basis is an immutable sequence:
            sage: type(O.basis())
            <class 'sage.structure.sequence.Sequence'>

        The generator functionality uses the basis method:
            sage: O.0
            1
            sage: O.1
            c
            sage: O.gens()
            [1, c, c^2]
            sage: O.ngens()
            3
        """
        try:
            return self.__basis
        except AttributeError:
            V, from_V, to_V = self._K.vector_space()
            B = Sequence([from_V(b) for b in self._module_rep.basis()], immutable=True)
            self.__basis = B
        return B

    def absolute_order(self):
        """
        Return the absolute order associated to this order, which is
        just this order again since this is an absolute order.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + 2)
            sage: O1 = K.order(a); O1
            Order with module basis 1, a, a^2 in Number Field in a with defining polynomial x^3 + 2
            sage: O1.absolute_order() is O1
            True
        """
        return self


class RelativeOrder(Order):
    """
    A relative order in a number field.

    A relative order is an order in some relative number field, and a
    specific choice of order in the base of the relative number field.

    Invariants of this order may be computed with respect to the
    contained order.
    """
    def __init__(self, K, absolute_order, base, is_maximal=None, check=True):
        """
        Create the relative order.
        """
        Order.__init__(self, K, is_maximal=is_maximal)
        self._absolute_order = absolute_order
        self._base = base

    def __call__(self, x):
        """
        Coerce an element into this relative order.
        """
        if x.parent() is not self._K:
            x = self._K(x)
        x = self._absolute_order(x) # will test membership
        return OrderElement_relative(self, x)

    def _repr_(self):
        """
        Return print representation of this relative order.
        """
        return "Relative Order with ZZ-module basis %s in %s" % (", ".join([str(b) for b in self.basis()]), self._K)

    def absolute_order(self):
        """
        Return underlying absolute order associated to this relative
        order.
        """
        return self._absolute_order

    def basis(self):
        """
        Return module basis for this relative order.  This is a list
        of elements that generate this order over the base order.

        WARNING: For now this basis is actually just a basis over ZZ.

        EXAMPLES:
            sage: K.<a,b> = NumberField([x^2+1, x^2+3])
            sage: O = K.order([a,b])
            sage: O.basis()
            [(-1/2*b - 5/2)*a + 3/2*b - 1/2,
             (-1/2*b - 7/2)*a + 2*b - 1,
             (-b)*a + -2,
             (-5)*a + 3*b]
            sage: z = O.0; z
            (-1/2*b - 5/2)*a + 3/2*b - 1/2
            sage: z.absolute_minpoly()
            x^4 + 2*x^3 + 26*x^2 - 20*x + 4
        """
        try:
            return self.__basis
        except AttributeError:
            pass
        O = self._absolute_order
        K = O.number_field()
        from_K, _ = K.structure()
        self.__basis = [from_K(a) for a in O.basis()]
        return self.__basis

    def __add__(left, right):
        """
        Add two relative orders or a relative order to an absolute
        order (which always results in an absolute order).

        EXAMPLES:
            sage: K.<a,b> = NumberField([x^2+1, x^2+3])
            sage: O2 = K.order([2*a]); O2.absolute_discriminant()
            2304
            sage: O3 = K.order([3*a, 2*b]); O3.absolute_discriminant()
            11664
            sage: O = (O2 + O3); O
            Relative Order with ZZ-module basis (-1/2*b - 5/2)*a + 3/2*b - 1/2, (-1/2*b - 7/2)*a + 2*b - 1, (-b)*a + -2, (-5)*a + 3*b in Number Field in a with defining polynomial x^2 + 1 over its base field
            sage: O.absolute_discriminant()
            144
            sage: O.is_suborder(O2)
            False
            sage: O2.is_suborder(O)
            True
            sage: O3.is_suborder(O)
            True
        """
        if isinstance(left, AbsoluteOrder):
            return left + right._absolute_order
        elif isinstance(right, AbsoluteOrder):
            return left._absolute_order + right
        elif isinstance(left, RelativeOrder) and isinstance(right, RelativeOrder):
            if left._K != right._K:
                raise TypeError, "Number fields don't match."
            if left._base != right._base:
                raise TypeError, "Bases don't match."
            return RelativeOrder(left._K, left._absolute_order + right._absolute_order,
                                 left._base, check=False)
        else:
            raise NotImplementedError

    def __and__(left, right):
        """
        Intersect two relative orders or a relative and absolute order
        (which always results in an absolute order).
        """
        if isinstance(left, AbsoluteOrder):
            return left & right._absolute_order
        elif isinstance(right, AbsoluteOrder):
            return left._absolute_order & right
        elif isinstance(left, RelativeOrder) and isinstance(right, RelativeOrder):
            if left._K != right._K:
                raise TypeError, "Number fields don't match."
            if left._base != right._base:
                raise TypeError, "Bases don't match."
            return RelativeOrder(left._K, left._absolute_order & right._absolute_order,
                                 left._base, check=False)
        else:
            raise NotImplementedError

    def absolute_discriminant(self):
        return self.absolute_order().discriminant()

    def is_suborder(self, other):
        return self.absolute_order().is_suborder(other.absolute_order())


def each_is_integral(v):
    """
    Return True if each element of the list v of elements of a number
    field is integral.

    EXAMPLES:
        sage: W.<sqrt5> = NumberField(x^2 - 5)
        sage: from sage.rings.number_field.order import each_is_integral
        sage: each_is_integral([sqrt5, 2, (1+sqrt5)/2])
        True
        sage: each_is_integral([sqrt5, (1+sqrt5)/3])
        False
    """
    for x in v:
        if not x.is_integral():
            return False
    return True

def absolute_order_from_ring_generators(gens, check_is_integral=True,
                                        check_rank=True, is_maximal=None,
                                        allow_subfield=False):
    """
    INPUT:
        gens -- list of integral elements of an absolute order.
        check_is_integral -- bool (default: True), whether to check
                             that each generator is integral.
        check_rank -- bool (default: True), whether to check that
                      the ring generated by gens is of full rank.
        is_maximal -- bool (or None); set if maximality of the generated order is known
        allow_subfield -- bool (default: False), if True and the generators
              do not generate an order, i.e., they generate a subring
              of smaller rank, instead of raising an error, return
              an order in a smaller number field.

    EXAMPLES:
        sage: K.<a> = NumberField(x^4 - 5)
        sage: K.order(a)
        Order with module basis 1, a, a^2, a^3 in Number Field in a with defining polynomial x^4 - 5

    We have to explicitly import this function, since typically
    it is called with \code{K.order} as above.
        sage: from sage.rings.number_field.order import absolute_order_from_ring_generators
        sage: absolute_order_from_ring_generators([a])
        Order with module basis 1, a, a^2, a^3 in Number Field in a with defining polynomial x^4 - 5
        sage: absolute_order_from_ring_generators([3*a, 2, 6*a+1])
        Order with module basis 1, 3*a, 9*a^2, 27*a^3 in Number Field in a with defining polynomial x^4 - 5

    If one of the inputs is non-integral, it is an error.
        sage: absolute_order_from_ring_generators([a/2])
        Traceback (most recent call last):
        ...
        ValueError: each generator must be integral

    If the gens do not generate an order, i.e., generate a ring of full
    rank, then it is an error.
        sage: absolute_order_from_ring_generators([a^2])
        Traceback (most recent call last):
        ...
        ValueError: the rank of the span of gens is wrong

    Both checking for integrality and checking for full rank can be
    turned off in order to save time, though one can get nonsense as
    illustrated below.
        sage: absolute_order_from_ring_generators([a/2], check_is_integral=False)
        Order with module basis 1, 1/2*a, 1/4*a^2, 1/8*a^3 in Number Field in a with defining polynomial x^4 - 5
        sage: absolute_order_from_ring_generators([a^2], check_rank=False)
        Order with module basis 1, a^2 in Number Field in a with defining polynomial x^4 - 5
    """
    if check_is_integral and not each_is_integral(gens):
        raise ValueError, "each generator must be integral"
    gens = Sequence(gens)
    K = gens.universe()
    n = K.degree()
    module_gens = monomials(gens, n)
    return absolute_order_from_module_generators(module_gens,
               check_integral=False, check_is_ring=False,
               check_rank=check_rank, is_maximal=is_maximal,
               allow_subfield = allow_subfield)


def absolute_order_from_module_generators(gens,
              check_integral=True, check_rank=True,
              check_is_ring=True, is_maximal=None,
              allow_subfield = False):
    """
    INPUT:
        gens -- list of elements of an absolute number field
                that generates an order in that number field as a ZZ
                *module*.
        check_integral -- check that each gen is integral
        check_rank -- check that the gens span a module of the correct rank
        check_is_ring -- check that the module is closed under multiplication
                         (this is very expensive)
        is_maximal -- bool (or None); set if maximality of the generated order is known

    OUTPUT:
        an absolute order

    EXAMPLES:
        sage: K.<a> = NumberField(x^4 - 5)

    We have to explicitly import the function, since it isn't meant for regular usage:
        sage: from sage.rings.number_field.order import absolute_order_from_module_generators
        sage: O = K.maximal_order(); O
        Order with module basis 1/2*a^2 + 1/2, 1/2*a^3 + 1/2*a, a^2, a^3 in Number Field in a with defining polynomial x^4 - 5
        sage: O.module()
        Free module of degree 4 and rank 4 over Integer Ring
        Echelon basis matrix:
        [1/2   0 1/2   0]
        [  0 1/2   0 1/2]
        [  0   0   1   0]
        [  0   0   0   1]
        sage: g = O.gens(); g
        [1/2*a^2 + 1/2, 1/2*a^3 + 1/2*a, a^2, a^3]
        sage: absolute_order_from_module_generators(g)
        Order with module basis 1/2*a^2 + 1/2, 1/2*a^3 + 1/2*a, a^2, a^3 in Number Field in a with defining polynomial x^4 - 5

    We illustrate each check flag -- the output is the same but in case the function would
    run ever so slightly faster:
        sage: absolute_order_from_module_generators(g,  check_is_ring=False)
        Order with module basis 1/2*a^2 + 1/2, 1/2*a^3 + 1/2*a, a^2, a^3 in Number Field in a with defining polynomial x^4 - 5
        sage: absolute_order_from_module_generators(g,  check_rank=False)
        Order with module basis 1/2*a^2 + 1/2, 1/2*a^3 + 1/2*a, a^2, a^3 in Number Field in a with defining polynomial x^4 - 5
        sage: absolute_order_from_module_generators(g,  check_integral=False)
        Order with module basis 1/2*a^2 + 1/2, 1/2*a^3 + 1/2*a, a^2, a^3 in Number Field in a with defining polynomial x^4 - 5

    Next we illustrate constructing "fake" to illustrate turning off various check flags:
        sage: k.<i> = NumberField(x^2 + 1)
        sage: absolute_order_from_module_generators([2, 2*i],  check_is_ring=False)
        Order with module basis 2, 2*i in Number Field in i with defining polynomial x^2 + 1
        sage: absolute_order_from_module_generators([k(1)],  check_rank=False)
        Order with module basis 1 in Number Field in i with defining polynomial x^2 + 1

    If the order contains a non-integral element, even if we don't check that, we'll
    find that the rank is wrong or that the order isn't closed under multiplication:
        sage: absolute_order_from_module_generators([1/2, i],  check_integral=False)
        Traceback (most recent call last):
        ...
        ValueError: the module span of the gens is not closed under multiplication.
        sage: absolute_order_from_module_generators([1/2, i],  check_is_ring=False, check_integral=False)
        Order with module basis 1/2, i in Number Field in i with defining polynomial x^2 + 1

    We turn off all check flags and make a really messed up order.
        sage: absolute_order_from_module_generators([1/2, i],  check_is_ring=False, check_integral=False, check_rank=False)
        Order with module basis 1/2, i in Number Field in i with defining polynomial x^2 + 1
    """
    if len(gens) == 0:
        raise ValueError, "gens must span an order over ZZ"
    gens = Sequence(gens)
    if check_integral and not each_is_integral(gens):
        raise ValueError, "each generator must be integral"

    K = gens.universe()
    V, from_V, to_V = K.vector_space()
    mod_gens = [to_V(x) for x in gens]
    ambient = ZZ**V.dimension()
    W = ambient.span(mod_gens)

    if allow_subfield:
        if W.rank() < K.degree():
            # We have to make the order in a smaller field.
            # We do this by choosing a random element of W,
            # moving it back to K, and checking that it defines
            # a field of degree equal to the degree of W.
            # Then we move everything into that field, where
            # W does define an order.
            while True:
                z = W.random_element()
                alpha = from_V(z)
                if alpha.minpoly() == W.rank():
                    break
            # Now alpha generates a subfield there W is an order
            # (with the right rank).
            # We move each element of W to this subfield.
            c = alpha.coordinates_in_terms_of_powers()

    elif check_rank:
        if W.rank() != K.degree():
            raise ValueError, "the rank of the span of gens is wrong"

    if check_is_ring:
        # Is there a faster way?
        alg = [to_V(x) for x in monomials(gens, K.degree())]
        if ambient.span(alg) != W:
            raise ValueError, "the module span of the gens is not closed under multiplication."

    return AbsoluteOrder(K, W, check=False, is_maximal=is_maximal)  # we have already checked everything





def relative_order_from_ring_generators(gens, base,
                                        check_is_integral=True,
                                        check_rank=True,
                                        ensure_contain_base=True,
                                        is_maximal = None,
                                        allow_subfield=False):
    """
    INPUT:
        gens -- list of integral elements of an absolute order.
        base -- order in the base field.
        check_is_integral -- bool (default: True), whether to check
                             that each generator is integral.
        check_rank -- bool (default: True), whether to check that
                      the ring generated by gens is of full rank.
        ensure_contain_base -- bool (default: True), whether to
                      generate over the base, i.e., to throw in
                      generators for the base.
        is_maximal -- bool (or None); set if maximality of the generated order is known

    EXAMPLES:

    """
    if check_is_integral and not each_is_integral(gens):
        raise ValueError, "each generator must be integral"
    gens = Sequence(gens)

    # The top number field that contains the order.
    K = gens.universe()

    # The absolute version of that field.
    Kabs = K.absolute_field('z')
    from_Kabs, to_Kabs = Kabs.structure()

    n = K.degree()
    module_gens = [to_Kabs(a) for a in gens]
    if ensure_contain_base:
        module_gens = module_gens + [to_Kabs(a) for a in base.gens()]
    absolute_order_module_gens = monomials(module_gens, n)

    abs_order =  absolute_order_from_module_generators(absolute_order_module_gens,
                                                       check_integral=False, check_is_ring=False,
                                                       check_rank=check_rank)

    return RelativeOrder(K, abs_order, base, check=False, is_maximal=is_maximal)
