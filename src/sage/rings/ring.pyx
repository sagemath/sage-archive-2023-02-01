"""
Rings

This module provides the abstract base class :class:`Ring` from which all rings
in Sage are derived, as well as a selection of more specific base classes. The
class inheritance hierarchy is:

- :class:`Ring`

  - :class:`Algebra`
  - :class:`CommutativeRing`

    - :class:`NoetherianRing`
    - :class:`CommutativeAlgebra`
    - :class:`IntegralDomain`

      - :class:`DedekindDomain`
      - :class:`PrincipalIdealDomain`

        - :class:`EuclideanDomain`
        - :class:`Field`

          - :class:`FiniteField`

Some aspects of this structure may seem strange, but this is an unfortunate
consequence of the fact that Cython classes do not support multiple
inheritance. Hence, for instance, :class:`Field` cannot be a subclass of both
:class:`NoetherianRing` and :class:`PrincipalIdealDomain`, although all fields
are Noetherian PIDs.

(A distinct but equally awkward issue is that sometimes we may not know *in
advance* whether or not a ring belongs in one of these classes; e.g. some
orders in number fields are Dedekind domains, but others are not, and we still
want to offer a unified interface, so orders are never instances of the
DedekindDomain class.)

AUTHORS:

- David Harvey (2006-10-16): changed :class:`CommutativeAlgebra` to derive from
  :class:`CommutativeRing` instead of from :class:`Algebra`
- David Loeffler (2009-07-09): documentation fixes, added to reference manual
"""

#*****************************************************************************
#       Copyright (C) 2005,2007 William Stein <wstein@gmail.com>
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

include "../ext/stdsage.pxi"
include "../ext/python_bool.pxi"

import re

from sage.structure.parent_gens cimport ParentWithGens
from sage.structure.parent cimport Parent
from sage.misc.prandom import randint, randrange

cdef class Ring(ParentWithGens):
    """
    Generic ring class.
    """
    def __iter__(self):
        r"""
        Return an iterator through the elements of self. Not implemented in general.

        EXAMPLES::

            sage: sage.rings.ring.Ring.__iter__(ZZ)
            Traceback (most recent call last):
            ...
            NotImplementedError: object does not support iteration
        """
        raise NotImplementedError, "object does not support iteration"

    def __len__(self):
        r"""
        Return the cardinality of this ring if it is finite, else raise a TypeError.

        EXAMPLE::

            sage: len(Integers(24))
            24
            sage: len(RR)
            Traceback (most recent call last):
            ...
            TypeError: len() of unsized object
        """
        if self.is_finite():
            return self.cardinality()
        raise TypeError, 'len() of unsized object'

    def __getitem__(self, x):
        """
        Create a polynomial or power series ring over self and inject
        the variables into the global module scope.

        If x is an algebraic element, this will return an extension of self
        that contains x.

        EXAMPLES:

        We create several polynomial rings::

            sage: ZZ['x']
            Univariate Polynomial Ring in x over Integer Ring
            sage: QQ['x']
            Univariate Polynomial Ring in x over Rational Field
            sage: GF(17)['abc']
            Univariate Polynomial Ring in abc over Finite Field of size 17
            sage: GF(17)['a,b,c']
            Multivariate Polynomial Ring in a, b, c over Finite Field of size 17

        We can also create power series rings (in one variable) by
        using double brackets::

            sage: QQ[['t']]
            Power Series Ring in t over Rational Field
            sage: ZZ[['W']]
            Power Series Ring in W over Integer Ring

        Use ``Frac`` (for fraction field) to obtain a Laurent series ring::

            sage: Frac(QQ[['t']])
            Laurent Series Ring in t over Rational Field

        This can be used to create number fields too::

            sage: QQ[I]
            Number Field in I with defining polynomial x^2 + 1
            sage: QQ[sqrt(2)]
            Number Field in sqrt2 with defining polynomial x^2 - 2
            sage: QQ[sqrt(2),sqrt(3)]
            Number Field in sqrt2 with defining polynomial x^2 - 2 over its base field


        and orders in number fields::

            sage: ZZ[I]
            Order in Number Field in I with defining polynomial x^2 + 1
            sage: ZZ[sqrt(5)]
            Order in Number Field in sqrt5 with defining polynomial x^2 - 5
            sage: ZZ[sqrt(2)+sqrt(3)]
            Order in Number Field in a with defining polynomial x^4 - 10*x^2 + 1
        """

        from sage.rings.polynomial.polynomial_element import is_Polynomial
        if is_Polynomial(x):
            x = str(x)

        if not isinstance(x, str):
            if isinstance(x, tuple):
                v = x
            else:
                v = (x,)

            minpolys = None
            try:
                minpolys = [a.minpoly() for a in v]
            except (AttributeError, NotImplementedError, ValueError, TypeError), err:
                pass

            if minpolys:
                R = self
                # how to pass in names?
                # TODO: set up embeddings
                name_chr = 97 # a

                if len(minpolys) > 1:
                    w = []
                    names = []
                    for poly, var in zip(minpolys, v):
                        w.append(poly)
                        n, name_chr = gen_name(repr(var), name_chr)
                        names.append(n)
                else:
                    w = minpolys
                    name, name_chr = gen_name(repr(v[0]), name_chr)
                    names = [name]

                names = tuple(names)
                if len(w) > 1:
                    try:
                        # Doing the extension all at once is best, if possible.
                        return R.extension(w, names)
                    except (TypeError, ValueError):
                        pass
                for poly, var in zip(w, names):
                    R = R.extension(poly, var)
                return R

        if not isinstance(x, list):
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            P = PolynomialRing(self, x)
            return P

        P = None
        if isinstance(x, list):
            if len(x) != 1:
                raise NotImplementedError, "Power series rings only implemented in 1 variable"
            x = (str(x[0]), )
            from sage.rings.power_series_ring import PowerSeriesRing
            P = PowerSeriesRing

        # TODO: is this code ever used? Should it be?

        elif isinstance(x, (tuple, str)):
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            P = PolynomialRing
            if isinstance(x, tuple):
                y = []
                for w in x:
                    y.append(str(w))
                x = tuple(y)

        else:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            P = PolynomialRing
            x = (str(x),)

        if P is None:
            raise NotImplementedError

        if isinstance(x, tuple):
            v = x
        else:
            v = x.split(',')

        if len(v) > 1:
            R = P(self, len(v), names=v)
        else:
            R = P(self, x)

        return R

    def __xor__(self, n):
        r"""
        Trap the operation ^. It's next to impossible to test this since ^ is
        intercepted first by the preparser.

        EXAMPLE::

            sage: RR^3 # not tested
        """
        raise RuntimeError, "Use ** for exponentiation, not '^', which means xor\n"+\
              "in Python, and has the wrong precedence."

    def base_extend(self, R):
        """
        EXAMPLES::

            sage: QQ.base_extend(GF(7))
            Traceback (most recent call last):
            ...
            TypeError: no base extension defined
            sage: ZZ.base_extend(GF(7))
            Finite Field of size 7
        """
        if R.has_coerce_map_from(self):
            return R
        raise TypeError, 'no base extension defined'

    def category(self):
        """
        Return the category to which this ring belongs.

        EXAMPLES::

            sage: QQ['x,y'].category()
            Category of rings
        """
        from sage.categories.all import Rings
        return Rings()

    def ideal(self, *x, **kwds):
        """
        Return the ideal defined by x, i.e., generated by x.

        INPUT:

        - ``*x`` -- list or tuple of generators (or several input arguments)
        - ``coerce`` -- bool (default: True); this must be a keyword argument.
          Only set it to False if you are certain that each generator is
          already in the ring.

        TODO: For noncommutative rings, distinguish between ideals, right
        ideals, and left ideals.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: R.ideal(x,y)
            Ideal (x, y) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: R.ideal(x+y^2)
            Ideal (y^2 + x) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: R.ideal( [x^3,y^3+x^3] )
            Ideal (x^3, x^3 + y^3) of Multivariate Polynomial Ring in x, y over Rational Field
        """
        C = self._ideal_class_()
        if len(x) == 1 and isinstance(x[0], (list, tuple)):
            x = x[0]
        return C(self, x, **kwds)

    def __mul__(self, x):
        """
        Return the ideal x*R generated by x, where x is either an element
        or tuple or list of elements.

        TODO: For noncommutative rings, distinguish between ideals, right
        ideals, and left ideals.

        EXAMPLES::

            sage: R.<x,y,z> = GF(7)[]
            sage: (x+y)*R
            Ideal (x + y) of Multivariate Polynomial Ring in x, y, z over Finite Field of size 7
            sage: (x+y,z+y^3)*R
            Ideal (x + y, y^3 + z) of Multivariate Polynomial Ring in x, y, z over Finite Field of size 7
        """
        if isinstance(self, Ring):
            return self.ideal(x)
        else:
            return x.ideal(self)    # switched because this is Pyrex / extension class

# This will never get called -- _r_action is supposed to be an attribute of
# elements, not parents -- DL
#
#    def _r_action(self, x):
#        r"""
#        Return the ideal x * self.
#
#        TODO: For noncommutative rings, distinguish between ideals, right
#        ideals, and left ideals.
#
#        EXAMPLE::
#
#            sage: ZZ._r_action(3)
#
#        """
#        return self.ideal(x)

    def _ideal_class_(self):
        r"""
        Return a callable object that can be used to create ideals in this
        ring. For generic rings, this returns the factory function
        :func:`sage.rings.ideal.Ideal`, which does its best to be clever about
        what is required.

        TODO: For noncommutative rings, distinguish between ideals, right
        ideals, and left ideals.

        EXAMPLES::

            sage: RR._ideal_class_()
            <function Ideal at ...>
        """
        import sage.rings.ideal
        return sage.rings.ideal.Ideal

    def principal_ideal(self, gen, coerce=True):
        """
        Return the principal ideal generated by gen.

        EXAMPLES::

            sage: R.<x,y> = ZZ[]
            sage: R.principal_ideal(x+2*y)
            Ideal (x + 2*y) of Multivariate Polynomial Ring in x, y over Integer Ring
        """
        return self.ideal([gen], coerce=coerce)

    def unit_ideal(self):
        """
        Return the unit ideal of this ring.

        EXAMPLES::

            sage: Zp(7).unit_ideal()
            Principal ideal (1 + O(7^20)) of 7-adic Ring with capped relative precision 20
        """
        if self._unit_ideal is None:
            I = Ring.ideal(self, [self(1)], coerce=False)
            self._unit_ideal = I
            return I
        return self._unit_ideal

    def zero_ideal(self):
        """
        Return the zero ideal of this ring (cached).

        EXAMPLES::

            sage: ZZ.zero_ideal()
            Principal ideal (0) of Integer Ring
            sage: QQ.zero_ideal()
            Principal ideal (0) of Rational Field
            sage: QQ['x'].zero_ideal()
            Principal ideal (0) of Univariate Polynomial Ring in x over Rational Field

        The result is cached::

            sage: ZZ.zero_ideal() is ZZ.zero_ideal()
            True
        """
        if self._zero_ideal is None:
            I = Ring.ideal(self, [self(0)], coerce=False)
            self._zero_ideal = I
            return I
        return self._zero_ideal

    def zero_element(self):
        """
        Return the zero element of this ring (cached).

        EXAMPLES::

            sage: ZZ.zero_element()
            0
            sage: QQ.zero_element()
            0
            sage: QQ['x'].zero_element()
            0

        The result is cached::

            sage: ZZ.zero_element() is ZZ.zero_element()
            True
        """
        if self._zero_element is None:
            x = self(0)
            self._zero_element = x
            return x
        return self._zero_element

    def one_element(self):
        """
        Return the one element of this ring (cached), if it exists.

        EXAMPLES::

            sage: ZZ.one_element()
            1
            sage: QQ.one_element()
            1
            sage: QQ['x'].one_element()
            1

        The result is cached::

            sage: ZZ.one_element() is ZZ.one_element()
            True
        """
        if self._one_element is None:
            x = self(1)
            self._one_element = x
            return x
        return self._one_element

    def is_atomic_repr(self):
        """
        True if the elements have atomic string representations, in the sense
        that they print if they print at s, then -s means the negative of s.
        For example, integers are atomic but polynomials are not.

        EXAMPLES::

            sage: Zp(7).is_atomic_repr()
            False
            sage: QQ.is_atomic_repr()
            True
            sage: CDF.is_atomic_repr()
            False
        """
        return False

    def is_zero(self):
        """
        True if this is the zero ring.

        EXAMPLES::

            sage: Integers(1).is_zero()
            True
            sage: Integers(2).is_zero()
            False
            sage: QQ.is_zero()
            False
            sage: R.<x> = ZZ[]
            sage: R.quo(1).is_zero()
            True
            sage: R.<x> = GF(101)[]
            sage: R.quo(77).is_zero()
            True
            sage: R.quo(x^2+1).is_zero()
            False
        """
        return self.one_element() == self.zero_element()

    def is_commutative(self):
        """
        Return True if this ring is commutative.

        EXAMPLES::

            sage: QQ.is_commutative()
            True
            sage: QQ['x,y,z'].is_commutative()
            True
            sage: Q.<i,j,k> = QuaternionAlgebra(QQ, -1,-1)
            sage: Q.is_commutative()
            False
        """
        if self.is_zero():
            return True
        raise NotImplementedError

    def is_field(self):
        """
        Return True if this ring is a field.

        EXAMPLES::

            sage: QQ.is_field()
            True
            sage: GF(9,'a').is_field()
            True
            sage: ZZ.is_field()
            False
            sage: QQ['x'].is_field()
            False
            sage: Frac(QQ['x']).is_field()
            True
        """
        if self.is_zero():
            return False
        raise NotImplementedError

    cpdef bint is_exact(self) except -2:
        """
        Return True if elements of this ring are represented exactly, i.e.,
        there is no precision loss when doing arithmetic.

        .. note::

            This defaults to True, so even if it does return True you have
            no guarantee (unless the ring has properly overloaded this).

        EXAMPLES::

            sage: QQ.is_exact()
            True
            sage: ZZ.is_exact()
            True
            sage: Qp(7).is_exact()
            False
            sage: Zp(7, type='capped-abs').is_exact()
            False
        """
        return True

    def is_subring(self, other):
        """
        Return True if the canonical map from self to other is injective.

        Raises a NotImplementedError if not known.

        EXAMPLES::

            sage: ZZ.is_subring(QQ)
            True
            sage: ZZ.is_subring(GF(19))
            False
        """
        try:
            return self.Hom(other).natural_map().is_injective()
        except TypeError:
            return False

    def is_prime_field(self):
        r"""
        Return True if this ring is one of the prime fields `\QQ` or
        `\GF{p}`.

        EXAMPLES::

            sage: QQ.is_prime_field()
            True
            sage: GF(3).is_prime_field()
            True
            sage: GF(9,'a').is_prime_field()
            False
            sage: ZZ.is_prime_field()
            False
            sage: QQ['x'].is_prime_field()
            False
            sage: Qp(19).is_prime_field()
            False
        """
        return False

    def is_finite(self):
        """
        Return True if this ring is finite.

        EXAMPLES::

            sage: QQ.is_finite()
            False
            sage: GF(2^10,'a').is_finite()
            True
            sage: R.<x> = GF(7)[]
            sage: R.is_finite()
            False
            sage: S.<y> = R.quo(x^2+1)
            sage: S.is_finite()
            True
        """
        if self.is_zero():
            return True
        raise NotImplementedError

    def is_integral_domain(self):
        """
        Return True if this ring is an integral domain.

        EXAMPLES::

            sage: QQ.is_integral_domain()
            True
            sage: ZZ.is_integral_domain()
            True
            sage: ZZ['x,y,z'].is_integral_domain()
            True
            sage: Integers(8).is_integral_domain()
            False
            sage: Zp(7).is_integral_domain()
            True
            sage: Qp(7).is_integral_domain()
            True
        """
        if self.is_zero():
            return False
        return NotImplementedError

    def is_ring(self):
        """
        Return True since self is a ring.

        EXAMPLES::

            sage: QQ.is_ring()
            True
        """
        return True

    def is_noetherian(self):
        """
        Return True if this ring is Noetherian.

        EXAMPLES::

            sage: QQ.is_noetherian()
            True
            sage: ZZ.is_noetherian()
            True
        """
        raise NotImplementedError

    def characteristic(self):
        """
        Return the characteristic of this ring.

        EXAMPLES::

            sage: QQ.characteristic()
            0
            sage: GF(19).characteristic()
            19
            sage: Integers(8).characteristic()
            8
            sage: Zp(5).characteristic()
            0
        """
        from sage.rings.infinity import infinity
        from sage.rings.integer_ring import ZZ
        order_1 = self.one_element().additive_order()
        return ZZ.zero_element() if order_1 is infinity else order_1

    def order(self):
        """
        The number of elements of self.

        EXAMPLES::

            sage: GF(19).order()
            19
            sage: QQ.order()
            +Infinity
        """
        if self.is_zero():
            return 1
        raise NotImplementedError

    def __hash__(self):
        """
        EXAMPLES::

            sage: hash(QQ)
            -11115808
            sage: hash(ZZ)
            554590422
        """
        return hash(self.__repr__())

    def zeta(self, n=2, all=False):
        """
        Return an n-th root of unity in self if there is one,
        or raise an ArithmeticError otherwise.

        INPUT:

        - ``n`` -- positive integer
        - ``all`` -- bool, default: False.  If True, return a list of all n-th
          roots of 1.

        OUTPUT:

        element of self of finite order

        EXAMPLES::

            sage: QQ.zeta()
            -1
            sage: QQ.zeta(1)
            1
            sage: CyclotomicField(6).zeta()
            zeta6
            sage: CyclotomicField(3).zeta()
            zeta3
            sage: CyclotomicField(3).zeta().multiplicative_order()
            3
            sage: a = GF(7).zeta(); a
            3
            sage: a.multiplicative_order()
            6
            sage: a = GF(49,'z').zeta(); a
            z
            sage: a.multiplicative_order()
            48
            sage: a = GF(49,'z').zeta(2); a
            6
            sage: a.multiplicative_order()
            2
            sage: QQ.zeta(3)
            Traceback (most recent call last):
            ...
            ValueError: no n-th root of unity in rational field
            sage: Zp(7, prec=8).zeta()
            3 + 4*7 + 6*7^2 + 3*7^3 + 2*7^5 + 6*7^6 + 2*7^7 + O(7^8)
        """
        if n == 2:
            if all:
                return [self(-1)]
            else:
                return self(-1)
        elif n == 1:
            if all:
                return [self(1)]
            else:
                return self(1)
        else:
            f = self['x'].cyclotomic_polynomial(n)
            if all:
                return [-P[0] for P, e in f.factor() if P.degree() == 1]
            for P, e in f.factor():
                if P.degree() == 1:
                    return -P[0]
            raise ArithmeticError, "no %s-th root of unity in self"%n

    def zeta_order(self):
        """
        Return the order of the distinguished root of unity in self.

        EXAMPLES::

            sage: CyclotomicField(19).zeta_order()
            38
            sage: GF(19).zeta_order()
            18
            sage: GF(5^3,'a').zeta_order()
            124
            sage: Zp(7, prec=8).zeta_order()
            6
        """
        return self.zeta().multiplicative_order()

    def random_element(self, bound=2):
        """
        Return a random integer coerced into this ring, where the
        integer is chosen uniformly from the interval [-bound,bound].

        INPUT:

        - bound -- integer (default: 2)

        ALGORITHM:

        Uses Python's randint.

        TESTS:

        The following example returns a ``NotImplementedError`` since the
        generic ring class ``__call__`` function returns a
        ``NotImplementedError``. Note that
        ``sage.rings.ring.Ring.random_element`` performs a call in the generic
        ring class by a random integer::

            sage: R = sage.rings.ring.Ring(ZZ); R
            <type 'sage.rings.ring.Ring'>
            sage: R.random_element()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return self(randint(-bound,bound))


cdef class CommutativeRing(Ring):
    """
    Generic commutative ring.
    """
    def fraction_field(self):
        """
        Return the fraction field of self.

        EXAMPLES::

            sage: R = Integers(389)['x,y']
            sage: Frac(R)
            Fraction Field of Multivariate Polynomial Ring in x, y over Ring of integers modulo 389
            sage: R.fraction_field()
            Fraction Field of Multivariate Polynomial Ring in x, y over Ring of integers modulo 389
        """
        if self.is_field():
            return self
        elif not self.is_integral_domain():
            raise TypeError, "self must be an integral domain."
        if self.__fraction_field is not None:
            return self.__fraction_field
        else:
            import sage.rings.fraction_field
            K = sage.rings.fraction_field.FractionField_generic(self)
            self.__fraction_field = K
        return self.__fraction_field

    def _pseudo_fraction_field(self):
        r"""
        This method is used by the coercion model to determine if a / b should
        be treated as a * (1/b), for example when dividing an element of
        `\ZZ[x]` by an element of `\ZZ`.

        The default is to return the same value as ``self.fraction_field()``,
        but it may return some other domain in which division is usually
        defined (for example, ``\ZZ/n\ZZ`` for possibly composite `n`).

        EXAMPLES::

            sage: ZZ._pseudo_fraction_field()
            Rational Field
            sage: ZZ['x']._pseudo_fraction_field()
            Fraction Field of Univariate Polynomial Ring in x over Integer Ring
            sage: Integers(15)._pseudo_fraction_field()
            Ring of integers modulo 15
            sage: Integers(15).fraction_field()
            Traceback (most recent call last):
            ...
            TypeError: self must be an integral domain.
        """
        return self.fraction_field()

    def __pow__(self, n, _):
        """
        Return the free module of rank `n` over this ring.

        EXAMPLES::

            sage: QQ^5
            Vector space of dimension 5 over Rational Field
            sage: Integers(20)^1000
            Ambient free module of rank 1000 over Ring of integers modulo 20
        """
        import sage.modules.all
        return sage.modules.all.FreeModule(self, n)

    def is_commutative(self):
        """
        Return True, since this ring is commutative.

        EXAMPLES::

            sage: QQ.is_commutative()
            True
            sage: ZpCA(7).is_commutative()
            True
            sage: A = QuaternionAlgebra(QQ, -1, -3, names=('i','j','k')); A
            Quaternion Algebra (-1, -3) with base ring Rational Field
            sage: A.is_commutative()
            False
        """
        return True

    def krull_dimension(self):
        """
        Return the Krull dimension of this commutative ring.

        The Krull dimension is the length of the longest ascending chain
        of prime ideals.

        TESTS:

        ``krull_dimension`` is not implemented for generic commutative
        rings. Fields and PIDs, with Krull dimension equal to 0 and 1,
        respectively, have naive implementations of ``krull_dimension``.
        Orders in number fields also have Krull dimension 1::

            sage: R = CommutativeRing(ZZ)
            sage: R.krull_dimension()
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: QQ.krull_dimension()
            0
            sage: ZZ.krull_dimension()
            1
            sage: type(R); type(QQ); type(ZZ)
            <type 'sage.rings.ring.CommutativeRing'>
            <class 'sage.rings.rational_field.RationalField'>
            <type 'sage.rings.integer_ring.IntegerRing_class'>

        All orders in number fields have Krull dimension 1, including
        non-maximal orders::

            sage: K.<i> = QuadraticField(-1)
            sage: R = K.maximal_order(); R
            Maximal Order in Number Field in i with defining polynomial x^2 + 1
            sage: R.krull_dimension()
            1
            sage: R = K.order(2*i); R
            Order in Number Field in i with defining polynomial x^2 + 1
            sage: R.is_maximal()
            False
            sage: R.krull_dimension()
            1
        """
        raise NotImplementedError

    def ideal_monoid(self):
        """
        Return the monoid of ideals of this ring.

        EXAMPLES::

            sage: ZZ.ideal_monoid()
            Monoid of ideals of Integer Ring
            sage: R.<x>=QQ[]; R.ideal_monoid()
            Monoid of ideals of Univariate Polynomial Ring in x over Rational Field
        """
        if self.__ideal_monoid is not None:
            return self.__ideal_monoid
        else:
            from sage.rings.ideal_monoid import IdealMonoid
            M = IdealMonoid(self)
            #try:
            self.__ideal_monoid = M
            #except AttributeError:   # for pyrex classes
            #    pass
            return M

    def quotient(self, I, names=None):
        """
        Create the quotient of R by the ideal I.

        INPUT:

        - ``R`` -- a commutative ring
        - ``I`` -- an ideal of R
        - ``names`` -- (optional) names of the generators of the quotient (if
          there are multiple generators, you can specify a single character
          string and the generators are named in sequence starting with 0).

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
            sage: S = R.quotient(I, 'a')
            sage: S.gens()
            (a,)

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: S.<a,b> = R.quotient((x^2, y))
            sage: S
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2, y)
            sage: S.gens()
            (a, 0)
            sage: a == b
            False
        """
        import sage.rings.quotient_ring
        return sage.rings.quotient_ring.QuotientRing(self, I, names=names)

    def quo(self, I, names=None):
        """
        Create the quotient of R by the ideal I.  This is a synonym for
        :meth:`.quotient`

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: S.<a,b> = R.quo((x^2, y))
            sage: S
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2, y)
            sage: S.gens()
            (a, 0)
            sage: a == b
            False
        """
        return self.quotient(I, names=names)

    def extension(self, poly, name=None, names=None, embedding=None):
        """
        Algebraically extends self by taking the quotient self[x] / (f(x)).

        INPUT:

        - ``poly`` -- A polynomial whose coefficients are coercible into self
        - ``name`` -- (optional) name for the root of f

        EXAMPLES::

            sage: R = QQ['x']
            sage: y = polygen(R)
            sage: R.extension(y^2-5, 'a')
            Univariate Quotient Polynomial Ring in a over Univariate Polynomial Ring in x over Rational Field with modulus a^2 - 5
        """
        if name is None and names is not None:
            name = names
        elif name is None:
            name = str(poly.parent().gen(0))
        if embedding is not None:
            raise NotImplementedError
        R = self[str(name)]
        I = R.ideal(R(poly.list()))
        return R.quotient(I, name)

    def __div__(self, I):
        """
        Dividing one ring by another is not supported because there is no good
        way to specify generator names.

        EXAMPLES::

            sage: QQ / ZZ
            Traceback (most recent call last):
            ...
            TypeError: Use self.quo(I) or self.quotient(I) to construct the quotient ring.
        """
        raise TypeError, "Use self.quo(I) or self.quotient(I) to construct the quotient ring."
        #return self.quotient(I, names=None)

    def quotient_ring(self, I, names=None):
        """
        Return the quotient of self by the ideal I of self.
        (Synonym for self.quotient(I).)

        INPUT:

        - ``self`` -- a ring R
        - ``I`` -- an ideal of R
        - ``names`` -- (optional) names of the generators of the quotient. (If
          there are multiple generators, you can specify a single character
          string and the generators are named in sequence starting with 0.)

        OUTPUT:

        - R/I -- the quotient ring of R by the ideal I

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
            sage: S = R.quotient_ring(I, 'a')
            sage: S.gens()
            (a,)

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: S.<a,b> = R.quotient_ring((x^2, y))
            sage: S
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2, y)
            sage: S.gens()
            (a, 0)
            sage: a == b
            False
        """
        return self.quotient(I, names)


cdef class IntegralDomain(CommutativeRing):
    """
    Generic integral domain class.
    """
    def is_integral_domain(self):
        """
        Return True, since this ring is an integral domain.

        (This is a naive implementation for objects with type
        ``IntegralDomain``)

        EXAMPLES::

            sage: ZZ.is_integral_domain(); QQ.is_integral_domain(); ZZ[x].is_integral_domain()
            True
            True
            True
            sage: R = ZZ.quotient(ZZ.ideal(10)); R.is_integral_domain()
            False
        """
        return True

    def is_integrally_closed(self):
        r"""
        Return True if this ring is integrally closed in its field of
        fractions; otherwise return False.

        When no algorithm is implemented for this, then this
        function raises a NotImplementedError.

        Note that ``is_integrally_closed`` has a naive implementation
        in fields. For every field `F`, `F` is its own field of fractions,
        hence every element of `F` is integral over `F`.

        EXAMPLES::

            sage: ZZ.is_integrally_closed()
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: QQ.is_integrally_closed()
            True
            sage: QQbar.is_integrally_closed()
            True
        """
        raise NotImplementedError

    def is_field(self):
        r"""
        Return True if this ring is a field.

        EXAMPLES::

            sage: GF(7).is_field()
            True

        The following examples have their own ``is_field`` implementations::

            sage: ZZ.is_field(); QQ.is_field()
            False
            True
            sage: R.<x> = PolynomialRing(QQ); R.is_field()
            False

        An example where we raise a NotImplementedError::

            sage: R = IntegralDomain(ZZ)
            sage: R.is_field()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if self.is_finite():
            return True
        raise NotImplementedError, "unable to determine whether or not is a field."

cdef class NoetherianRing(CommutativeRing):
    """
    Generic Noetherian ring class.

    A Noetherian ring is a commutative ring in which every ideal is
    finitely generated.

    At present this is not actually used anywhere in the Sage code base
    (largely because of the lack of multiple inheritance for Cython classes).
    """
    def is_noetherian(self):
        """
        Return True since this ring is Noetherian.

        EXAMPLES::

            sage: ZZ.is_noetherian()
            True
            sage: QQ.is_noetherian()
            True
            sage: R.<x> = PolynomialRing(QQ)
            sage: R.is_noetherian()
            True
        """
        return True

cdef class DedekindDomain(IntegralDomain):
    """
    Generic Dedekind domain class.

    A Dedekind domain is a Noetherian integral domain of Krull
    dimension one that is integrally closed in its field of fractions.
    """
    def krull_dimension(self):
        """
        Return 1 since Dedekind domains have Krull dimension 1.

        EXAMPLES:

        The following are examples of Dedekind domains (Noetherian integral
        domains of Krull dimension one that are integrally closed over its
        field of fractions)::

            sage: ZZ.krull_dimension()
            1
            sage: K = NumberField(x^2 + 1, 's')
            sage: OK = K.ring_of_integers()
            sage: OK.krull_dimension()
            1

        The following are not Dedekind domains but have
        a ``krull_dimension`` function::

            sage: QQ.krull_dimension()
            0
            sage: T.<x,y> = PolynomialRing(QQ,2); T
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: T.krull_dimension()
            2
            sage: U.<x,y,z> = PolynomialRing(ZZ,3); U
            Multivariate Polynomial Ring in x, y, z over Integer Ring
            sage: U.krull_dimension()
            4

            sage: K.<i> = QuadraticField(-1)
            sage: R = K.order(2*i); R
            Order in Number Field in i with defining polynomial x^2 + 1
            sage: R.is_maximal()
            False
            sage: R.krull_dimension()
            1
        """
        return 1

    def is_integrally_closed(self):
        """
        Return True since Dedekind domains are integrally closed.

        EXAMPLES:

        The following are examples of Dedekind domains (Noetherian integral
        domains of Krull dimension one that are integrally closed over its
        field of fractions). (Note that the ring of integers does not have
        an implemented ``is_integrally_closed`` function.)

        ::

            sage: ZZ.is_integrally_closed()
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: K = NumberField(x^2 + 1, 's')
            sage: OK = K.ring_of_integers()
            sage: OK.is_integrally_closed()
            True

        These, however, are not Dedekind domains::

            sage: QQ.is_integrally_closed()
            True
            sage: S = ZZ[sqrt(5)]; S.is_integrally_closed()
            False
            sage: T.<x,y> = PolynomialRing(QQ,2); T
            Multivariate Polynomial Ring in x, y over Rational Field
            sage: T.is_integral_domain()
            True
        """
        return True

    def integral_closure(self):
        r"""
        Return self since Dedekind domains are integrally closed.

        EXAMPLES::

            sage: K = NumberField(x^2 + 1, 's')
            sage: OK = K.ring_of_integers()
            sage: OK.integral_closure()
            Maximal Order in Number Field in s with defining polynomial x^2 + 1
            sage: OK.integral_closure() == OK
            True

            sage: QQ.integral_closure() == QQ
            True
        """
        return self

    def is_noetherian(self):
        r"""
        Return True since Dedekind domains are Noetherian.

        EXAMPLES:

        The integers, `\ZZ`, and rings of integers of number
        fields are Dedekind domains::

            sage: ZZ.is_noetherian()
            True
            sage: K = NumberField(x^2 + 1, 's')
            sage: OK = K.ring_of_integers()
            sage: OK.is_noetherian()
            True
            sage: QQ.is_noetherian()
            True
        """
        return True


cdef class PrincipalIdealDomain(IntegralDomain):
    """
    Generic principal ideal domain.
    """
    def class_group(self):
        """
        Return the trivial group, since the class group of a PID is trivial.

        EXAMPLES::

            sage: QQ.class_group()
            Trivial Abelian Group
        """
        from sage.groups.abelian_gps.abelian_group import AbelianGroup
        return AbelianGroup([])

    def gcd(self, x, y, coerce=True):
        r"""
        Return the greatest common divisor of x and y, as elements
        of self.

        EXAMPLES:

        The integers are a principal ideal domain and hence a GCD domain::

            sage: ZZ.gcd(42, 48)
            6
            sage: 42.factor(); 48.factor()
            2 * 3 * 7
            2^4 * 3
            sage: ZZ.gcd(2^4*7^2*11, 2^3*11*13)
            88
            sage: 88.factor()
            2^3 * 11

        In a field, any nonzero element is a GCD of any nonempty set
        of elements.  For concreteness, Sage returns 1 in these cases::

            sage: QQ.gcd(ZZ(42), ZZ(48)); type(QQ.gcd(ZZ(42), ZZ(48)))
            1
            <type 'sage.rings.rational.Rational'>
            sage: QQ.gcd(1/2, 1/3)
            1

        Polynomial rings over fields are GCD domains as well. Here is a simple
        example over the ring of polynomials over the rationals as well as
        over an extension ring. Note that ``gcd`` requires x and y to be
        coercible::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = NumberField(x^2 - 2, 'a')
            sage: f = (x - a)*(x + a); g = (x - a)*(x^2 - 2)
            sage: print f; print g
            x^2 - 2
            x^3 - a*x^2 - 2*x + 2*a
            sage: f in R
            True
            sage: g in R
            False
            sage: R.gcd(f,g)
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce 2*a to a rational
            sage: R.base_extend(S).gcd(f,g)
            x^2 - 2
            sage: R.base_extend(S).gcd(f, (x - a)*(x^2 - 3))
            x - a
        """
        if coerce:
            x = self(x)
            y = self(y)
        return x.gcd(y)

    def content(self, x, y, coerce=True):
        r"""
        Return the content of x and y, i.e. the unique element c of
        self such that x/c and y/c are coprime and integral.

        EXAMPLES::

            sage: QQ.content(ZZ(42), ZZ(48)); type(QQ.content(ZZ(42), ZZ(48)))
            6
            <type 'sage.rings.rational.Rational'>
            sage: QQ.content(1/2, 1/3)
            1/6
            sage: factor(1/2); factor(1/3); factor(1/6)
            2^-1
            3^-1
            2^-1 * 3^-1
            sage: a = (2*3)/(7*11); b = (13*17)/(19*23)
            sage: factor(a); factor(b); factor(QQ.content(a,b))
            2 * 3 * 7^-1 * 11^-1
            13 * 17 * 19^-1 * 23^-1
            7^-1 * 11^-1 * 19^-1 * 23^-1

        Note the changes to the second entry::

            sage: c = (2*3)/(7*11); d = (13*17)/(7*19*23)
            sage: factor(c); factor(d); factor(QQ.content(c,d))
            2 * 3 * 7^-1 * 11^-1
            7^-1 * 13 * 17 * 19^-1 * 23^-1
            7^-1 * 11^-1 * 19^-1 * 23^-1
            sage: e = (2*3)/(7*11); f = (13*17)/(7^3*19*23)
            sage: factor(e); factor(f); factor(QQ.content(e,f))
            2 * 3 * 7^-1 * 11^-1
            7^-3 * 13 * 17 * 19^-1 * 23^-1
            7^-3 * 11^-1 * 19^-1 * 23^-1
        """
        if coerce:
            x = self(x)
            y = self(y)
        return x.content(y)

cdef class EuclideanDomain(PrincipalIdealDomain):
    """
    Generic Euclidean domain class.
    """
    def parameter(self):
        """
        Return an element of degree 1.

        EXAMPLES::

            sage: R.<x>=QQ[]
            sage: R.parameter()
            x
       """
        raise NotImplementedError

def is_Field(x):
    """
    Return True if x is a field.

    EXAMPLES::

        sage: from sage.rings.ring import is_Field
        sage: is_Field(QQ)
        True
        sage: is_Field(ZZ)
        False
        sage: is_Field(pAdicField(2))
        True
        sage: is_Field(5)
        False
    """
    return isinstance(x, Field) or (hasattr(x, 'is_field') and x.is_field())

cdef class Field(PrincipalIdealDomain):
    """
    Generic field
    """
    def category(self):
        """
        Return the category of this field, which is the category
        of fields.

        EXAMPLES:

        Examples with fields::

            sage: QQ.category()
            Category of fields
            sage: RR.category()
            Category of fields
            sage: CC.category()
            Category of fields
            sage: R.<x> = PolynomialRing(ZZ)
            sage: Q = R.fraction_field()
            sage: Q.category()
            Category of fields

        Although fields themselves, number fields belong to the category
        of 'number fields'::

            sage: F = NumberField(x^2 + 1, 'i')
            sage: F.category()
            Category of number fields
        """
        from sage.categories.all import Fields
        return Fields()

    def fraction_field(self):
        """
        Return the fraction field of self.

        EXAMPLES:

        Since fields are their own field of fractions, we simply get the
        original field in return::

            sage: QQ.fraction_field()
            Rational Field
            sage: RR.fraction_field()
            Real Field with 53 bits of precision
            sage: CC.fraction_field()
            Complex Field with 53 bits of precision

            sage: F = NumberField(x^2 + 1, 'i')
            sage: F.fraction_field()
            Number Field in i with defining polynomial x^2 + 1
        """
        return self

    def _pseudo_fraction_field(self):
        """
        The fraction field of self is always available as self.

        EXAMPLES::

            sage: QQ._pseudo_fraction_field()
            Rational Field
            sage: K = GF(5)
            sage: K._pseudo_fraction_field()
            Finite Field of size 5
            sage: K._pseudo_fraction_field() is K
            True
        """
        return self

    def divides(self, x, y, coerce=True):
        """
        Return True if x divides y in this field (usually True in a
        field!).  If ``coerce`` is True (the default), first coerce x and
        y into self.

        EXAMPLES::

            sage: QQ.divides(2, 3/4)
            True
            sage: QQ.divides(0, 5)
            False
        """
        if coerce:
            x = self(x)
            y = self(y)
        if x.is_zero():
            return y.is_zero()
        return True

    def ideal(self, *gens, **kwds):
        """
        Return the ideal generated by gens.

        EXAMPLES::

            sage: QQ.ideal(2)
            Principal ideal (1) of Rational Field
            sage: QQ.ideal(0)
            Principal ideal (0) of Rational Field
        """
        if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
            gens = gens[0]
        if not isinstance(gens, (list, tuple)):
            gens = [gens]
        for x in gens:
            if not self(x).is_zero():
                return self.unit_ideal()
        return self.zero_ideal()

    def integral_closure(self):
        """
        Return this field, since fields are integrally closed in their
        fraction field.

        EXAMPLES::

            sage: QQ.integral_closure()
            Rational Field
            sage: Frac(ZZ['x,y']).integral_closure()
            Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring
        """
        return self

    def is_field(self):
        """
        Return True since this is a field.

        EXAMPLES::

            sage: Frac(ZZ['x,y']).is_field()
            True
        """
        return True

    def is_integrally_closed(self):
        """
        Return True since fields are trivially integrally closed in
        their fraction field (since they are their fraction field).

        EXAMPLES::

            sage: Frac(ZZ['x,y']).is_integrally_closed()
            True
        """
        return True

    def is_noetherian(self):
        """
        Return True since fields are Noetherian rings.

        EXAMPLES::

            sage: QQ.is_noetherian()
            True
        """
        return True

    def krull_dimension(self):
        """
        Return the Krull dimension of this field, which is 0.

        EXAMPLES::

            sage: QQ.krull_dimension()
            0
            sage: Frac(QQ['x,y']).krull_dimension()
            0
        """
        return 0

    def prime_subfield(self):
        """
        Return the prime subfield of self.

        EXAMPLES::

            sage: k = GF(9, 'a')
            sage: k.prime_subfield()
            Finite Field of size 3
        """
        if self.characteristic() == 0:
            import sage.rings.rational_field
            return sage.rings.rational_field.RationalField()
        else:
            import sage.rings.finite_field
            return sage.rings.finite_field.FiniteField(self.characteristic())

cdef class FiniteFieldIterator:
    r"""
    An iterator over a finite field. This should only be used when the field is
    an extension of a smaller field which already has a separate iterator function.
    """
    cdef object iter
    cdef FiniteField parent

    def __init__(self,FiniteField parent):
        r"""
        Standard init function.

        EXAMPLE::

            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = iter(FiniteField_ext_pari(9, 'a')) # indirect doctest
            sage: isinstance(k, sage.rings.ring.FiniteFieldIterator)
            True
        """
        self.parent = parent
        self.iter =iter(self.parent.vector_space())

    def __next__(self):
        r"""
        Return the next element in the iterator.

        EXAMPLE::

            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = iter(FiniteField_ext_pari(9, 'a'))
            sage: k.next() # indirect doctest
            0
        """
        return self.parent(self.iter.next())

cdef class FiniteField(Field):
#    def __init__(self):
#        """
#        EXAMPLES::
#
#            sage: K = GF(7); K
#            Finite Field of size 7
#            sage: loads(K.dumps()) == K
#            True
#            sage: GF(7^10, 'a')
#            Finite Field in a of size 7^10
#            sage: K = GF(7^10, 'a'); K
#            Finite Field in a of size 7^10
#            sage: loads(K.dumps()) == K
#            True
#        """
#        raise NotImplementedError

    def __repr__(self):
        """
        String representation of this finite field.

        EXAMPLES::

            sage: k = GF(127)
            sage: k # indirect doctest
            Finite Field of size 127

            sage: k.<b> = GF(2^8)
            sage: k
            Finite Field in b of size 2^8

            sage: k.<c> = GF(2^20)
            sage: k
            Finite Field in c of size 2^20

            sage: k.<d> = GF(7^20)
            sage: k
            Finite Field in d of size 7^20
        """
        if self.degree()>1:
            return "Finite Field in %s of size %s^%s"%(self.variable_name(),self.characteristic(),self.degree())
        else:
            return "Finite Field of size %s"%(self.characteristic())

    def _latex_(self):
        r"""
        Returns a string denoting the name of the field in LaTeX.

        The :func:`sage.misc.latex.latex` function calls the
        ``_latex_`` attribute when available.

        EXAMPLES:

        The ``latex`` command parses the string::

            sage: GF(81, 'a')._latex_()
            '\\Bold{F}_{3^{4}}'
            sage: latex(GF(81, 'a'))
            \Bold{F}_{3^{4}}
            sage: GF(3)._latex_()
            '\\Bold{F}_{3}'
            sage: latex(GF(3))
            \Bold{F}_{3}
        """
        if self.degree() > 1:
            e = "^{%s}"%self.degree()
        else:
            e = ""
        return "\\Bold{F}_{%s%s}"%(self.characteristic(), e)

    def _gap_init_(self):
        """
        Return string that initializes the GAP version of
        this finite field.

        EXAMPLES::

            sage: GF(9,'a')._gap_init_()
            'GF(9)'
        """
        return 'GF(%s)'%self.order()

    def _magma_init_(self, magma):
        """
        Return string representation of self that Magma can
        understand.

        EXAMPLES::

            sage: GF(97,'a')._magma_init_(magma)              # optional - magma
            'GF(97)'
            sage: GF(9,'a')._magma_init_(magma)               # optional - magma
            'SageCreateWithNames(ext<GF(3)|_sage_[...]![GF(3)!2,GF(3)!2,GF(3)!1]>,["a"])'
            sage: magma(GF(9,'a'))                            # optional - magma
            Finite field of size 3^2
            sage: magma(GF(9,'a')).1                          # optional - magma
            a
        """
        if self.degree() == 1:
            return 'GF(%s)'%self.order()
        B = self.base_ring()
        p = self.polynomial()
        s = "ext<%s|%s>"%(B._magma_init_(magma),p._magma_init_(magma))
        return magma._with_names(s, self.variable_names())

    def _macaulay2_init_(self):
        """
        Returns the string representation of self that Macaulay2 can
        understand.

        EXAMPLES::

            sage: GF(97,'a')._macaulay2_init_()
            'GF 97'

            sage: macaulay2(GF(97, 'a'))       # optional - macaulay2
            ZZ
            --
            97
            sage: macaulay2(GF(49, 'a'))       # optional - macaulay2
            GF 49
        """
        return "GF %s"%(self.order())

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(GF(5), verify=True)
            # Verified
            GF(5)
            sage: sage_input(GF(32, 'a'), verify=True)
            # Verified
            R.<x> = GF(2)[]
            GF(2^5, 'a', x^5 + x^2 + 1)
            sage: K = GF(125, 'b')
            sage: sage_input((K, K), verify=True)
            # Verified
            R.<x> = GF(5)[]
            GF_5_3 = GF(5^3, 'b', x^3 + 3*x + 3)
            (GF_5_3, GF_5_3)
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: GF(81, 'a')._sage_input_(SageInputBuilder(), False)
            {call: {atomic:GF}({binop:** {atomic:3} {atomic:4}}, {atomic:'a'}, {binop:+ {binop:+ {binop:** {gen:x {constr_parent: {subscr: {call: {atomic:GF}({atomic:3})}[{atomic:'x'}]} with gens: ('x',)}} {atomic:4}} {binop:* {atomic:2} {binop:** {gen:x {constr_parent: {subscr: {call: {atomic:GF}({atomic:3})}[{atomic:'x'}]} with gens: ('x',)}} {atomic:3}}}} {atomic:2}})}
        """
        if self.degree() == 1:
            v = sib.name('GF')(sib.int(self.characteristic()))
            name = 'GF_%d' % self.characteristic()
        else:
            v = sib.name('GF')(sib.int(self.characteristic()) ** sib.int(self.degree()),
                               self.variable_name(),
                               self.modulus())
            name = 'GF_%d_%d' % (self.characteristic(), self.degree())
        sib.cache(self, v, name)
        return v

    cdef int _cmp_c_impl(left, Parent right) except -2:
        """
        Compares this finite field with other.

        .. warning::

            The notation of equality of finite fields in Sage is
            currently not stable, i.e., it may change in a future version.

        EXAMPLES::

            sage: FiniteField(3**2, 'c') == FiniteField(3**3, 'c')
            False
            sage: FiniteField(3**2, 'c') == FiniteField(3**2, 'c')
            True

        The variable name is (currently) relevant for comparison of finite fields::

            sage: FiniteField(3**2, 'c') == FiniteField(3**2, 'd')
            False
        """
        if not PY_TYPE_CHECK(right, FiniteField):
            return cmp(type(left), type(right))
        c = cmp(left.characteristic(), right.characteristic())
        if c: return c
        c = cmp(left.degree(), right.degree())
        if c: return c
        # TODO comparing the polynomials themselves would recursively call
        # this cmp...  Also, as mentioned above, we will get rid of this.
        if left.degree() > 1:
            c = cmp(str(left.polynomial()), str(right.polynomial()))
            if c: return c
        return 0

    def __iter__(self):
        """
        Return an iterator over the elements of this finite field. This generic
        implementation uses the fairly simple :class:`FiniteFieldIterator`
        class; derived classes may implement their own more sophisticated
        replacements.

        EXAMPLE::

            sage: from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(8, 'a')
            sage: i = iter(k); i # indirect doctest
            <sage.rings.ring.FiniteFieldIterator object at ...>
            sage: i.next()
            0
            sage: list(k) # indirect doctest
            [0, 1, a, a + 1, a^2, a^2 + 1, a^2 + a, a^2 + a + 1]
        """
        return FiniteFieldIterator(self)

    def _is_valid_homomorphism_(self, codomain, im_gens):
        """
        Return True if the map from self to codomain sending
        self.0 to the unique element of im_gens is a valid field
        homomorphism. Otherwise, return False.

        EXAMPLES::

            sage: k = FiniteField(73^2, 'a')
            sage: K = FiniteField(73^3, 'b') ; b = K.0
            sage: L = FiniteField(73^4, 'c') ; c = L.0
            sage: k.hom([c]) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: images do not define a valid homomorphism

            sage: k.hom([c^(73*73+1)])
            Ring morphism:
            From: Finite Field in a of size 73^2
            To:   Finite Field in c of size 73^4
            Defn: a |--> 7*c^3 + 13*c^2 + 65*c + 71

            sage: k.hom([b])
            Traceback (most recent call last):
            ...
            TypeError: images do not define a valid homomorphism
        """

        if (self.characteristic() != codomain.characteristic()):
            raise ValueError, "no map from %s to %s"%(self, codomain)
        if (len(im_gens) != 1):
            raise ValueError, "only one generator for finite fields."

        return (im_gens[0].charpoly())(self.gen(0)).is_zero()

    def _Hom_(self, codomain, cat=None):
        """
        Return homset of homomorphisms from self to the finite field codomain.
        This function is implicitly called by the Hom method or function.

        The cat option is currently ignored.

        EXAMPLES::

            sage: K.<a> = GF(25); K
            Finite Field in a of size 5^2
            sage: K.Hom(K) # indirect doctest
            Automorphism group of Finite Field in a of size 5^2
        """
        from sage.rings.finite_field_morphism import FiniteFieldHomset
        return FiniteFieldHomset(self, codomain)

    def gen(self):
        r"""
        Return a generator of this field (over its prime field). As this is an abstract base class,
        this just raises a NotImplementedError.

        EXAMPLE::

            sage: K = GF(17)
            sage: sage.rings.ring.FiniteField.gen(K)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def zeta_order(self):
        """
        Return the order of the distinguished root of unity in self.

        EXAMPLES::

            sage: GF(9,'a').zeta_order()
            8
            sage: GF(9,'a').zeta()
            a
            sage: GF(9,'a').zeta().multiplicative_order()
            8
        """
        return self.order() - 1

    def zeta(self, n=None):
        """
        Returns an element of multiplicative order n in this
        finite field, if there is one.  Raises a ValueError if there
        is not.

        EXAMPLES::

            sage: k = GF(7)
            sage: k.zeta()
            3
            sage: k.zeta().multiplicative_order()
            6
            sage: k.zeta(3)
            2
            sage: k.zeta(3).multiplicative_order()
            3
            sage: k = GF(49, 'a')
            sage: k.zeta().multiplicative_order()
            48
            sage: k.zeta(6)
            3

        Even more examples::

            sage: GF(9,'a').zeta_order()
            8
            sage: GF(9,'a').zeta()
            a
            sage: GF(9,'a').zeta(4)
            a + 1
            sage: GF(9,'a').zeta()^2
            a + 1
        """
        z = self.multiplicative_generator()
        if n is None:
            return z
        else:
            import sage.rings.integer
            n = sage.rings.integer.Integer(n)
            m = z.multiplicative_order()
            if m % n != 0:
                raise ValueError, "No %sth root of unity in self"%n
            return z**(m.__floordiv__(n))

    def multiplicative_generator(self):
        """
        Return a generator for the multiplicative group of this field.

        This generator might change from one version of Sage to
        another.

        EXAMPLES::

            sage: k = GF(997)
            sage: k.multiplicative_generator()
            7
            sage: k = GF(11^3, name='a')
            sage: k.multiplicative_generator()
            a
        """
        from sage.rings.arith import primitive_root

        if self.__multiplicative_generator is not None:
            return self.__multiplicative_generator
        else:
            if self.degree() == 1:
                self.__multiplicative_generator = self(primitive_root(self.order()))
                return self.__multiplicative_generator
            n = self.order() - 1
            a = self.gen(0)
            if a.multiplicative_order() == n:
                self.__multiplicative_generator = a
                return a
            for a in self:
                if a == 0:
                    continue
                if a.multiplicative_order() == n:
                    self.__multiplicative_generator = a
                    return a

    def ngens(self):
        """
        The number of generators of the finite field.  Always 1.

        EXAMPLES::

            sage: k = FiniteField(3^4, 'b')
            sage: k.ngens()
            1
        """
        return 1

    def is_field(self):
        """
        Returns whether or not the finite field is a field, i.e.,
        always returns True.

        EXAMPLES::

            sage: k.<a> = FiniteField(3^4)
            sage: k.is_field()
            True
        """
        return True

    def is_finite(self):
        """
        Return True since a finite field is finite.

        EXAMPLES::

            sage: GF(997).is_finite()
            True
        """
        return True

    def order(self):
        """
        Return the order of this finite field.

        EXAMPLES::

            sage: GF(997).order()
            997
        """
        raise NotImplementedError

    def cardinality(self):
        """
        Return the order of this finite field (same as self.order()).

        EXAMPLES::

            sage: GF(997).cardinality()
            997
        """
        return self.order()

    def __len__(self):
        """
        len(k) returns the cardinality of k, i.e., the number of elements in k.

        EXAMPLE::

            sage: len(GF(8,'a'))
            8
        """
        return self.order()

    def is_prime_field(self):
        """
        Return True if self is a prime field, i.e., has degree 1.

        EXAMPLES::

            sage: GF(3^7, 'a').is_prime_field()
            False
            sage: GF(3, 'a').is_prime_field()
            True
        """
        return self.degree()==1

    def modulus(self):
        r"""
        Return the minimal polynomial of the generator of self (over an
        appropriate base field).

        The minimal polynomial of an element `a` in a field is the unique
        irreducible polynomial of smallest degree with coefficients in the base
        field that has `a` as a root. In finite field extensions, `\GF{p^n}`,
        the base field is `\GF{p}`. Here are several examples::

            sage: F.<a> = GF(7^2, 'a'); F
            Finite Field in a of size 7^2
            sage: F.polynomial_ring()
            Univariate Polynomial Ring in a over Finite Field of size 7
            sage: f = F.modulus(); f
            x^2 + 6*x + 3
            sage: f(a)
            0

        Although `f` is irreducible over the base field, we can double-check
        whether or not `f` factors in `F` as follows. The command
        `F[x](f)` coerces `f` as a polynomial with coefficients in `F`.
        (Instead of a polynomial with coefficients over the base field.)

        ::

            sage: f.factor()
            x^2 + 6*x + 3
            sage: F[x](f).factor()
            (x + a + 6) * (x + 6*a)

        Here is an example with a degree 3 extension::

            sage: G.<b> = GF(7^3, 'b'); G
            Finite Field in b of size 7^3
            sage: g = G.modulus(); g
            x^3 + 6*x^2 + 4
            sage: g.degree(); G.degree()
            3
            3
        """
        return self.polynomial_ring("x")(self.polynomial().list())

    def unit_group_exponent(self):
        """
        The exponent of the unit group of the finite field.  For a
        finite field, this is always the order minus 1.

        EXAMPLES::

            sage: k = GF(2^10, 'a')
            sage: k.order()
            1024
            sage: k.unit_group_exponent()
            1023
        """
        return self.order() - 1


    def random_element(self, bound=None):
        r"""
        A random element of the finite field.

        INPUT:

        - ``bound`` -- ignored (exists for consistency with other
          ``random_element`` methods, e.g. for `\ZZ`)

        EXAMPLES::

            sage: k = GF(2^10, 'a')
            sage: k.random_element() # random output
            a^9 + a
        """
        if self.degree() == 1:
            return self(randrange(self.order()))
        v = self.vector_space().random_element()
        return self(v)

    def polynomial(self):
        """
        Return the defining polynomial of this finite field.

        EXAMPLES::

            sage: f = GF(27,'a').polynomial(); f
            a^3 + 2*a + 1
            sage: parent(f)
            Univariate Polynomial Ring in a over Finite Field of size 3
        """
        raise NotImplementedError

    def polynomial_ring(self, variable_name=None):
        """
        Returns the polynomial ring over the prime subfield in the
        same variable as this finite field.

        EXAMPLES::

            sage: k.<alpha> = FiniteField(3^4)
            sage: k.polynomial_ring()
            Univariate Polynomial Ring in alpha over Finite Field of size 3
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.finite_field import FiniteField as GF

        if variable_name is None and self.__polynomial_ring is not None:
            return self.__polynomial_ring
        else:
            if variable_name is None:
                self.__polynomial_ring = PolynomialRing(GF(self.characteristic()), self.variable_name())
                return self.__polynomial_ring
            else:
                return PolynomialRing(GF(self.characteristic()), variable_name)

    def vector_space(self):
        """
        Return the vector space over the prime subfield isomorphic
        to this finite field as a vector space.

        EXAMPLES::

            sage: GF(27,'a').vector_space()
            Vector space of dimension 3 over Finite Field of size 3
        """
        if self.__vector_space is not None:
            return self.__vector_space
        else:
            import sage.modules.all
            V = sage.modules.all.VectorSpace(self.prime_subfield(),self.degree())
            self.__vector_space = V
            return V

    def __hash__(self):
        r"""
        Return a hash of this finite field.

        EXAMPLES::

            sage: hash(GF(17))
            -1709973406 # 32-bit
            9088054599082082 # 64-bit
        """
        return hash("GF") + hash(self.order())

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: A = FiniteField(127)
            sage: A == loads(dumps(A)) # indirect doctest
            True
            sage: B = FiniteField(3^3,'b')
            sage: B == loads(dumps(B))
            True
            sage: C = FiniteField(2^16,'c')
            sage: C == loads(dumps(C))
            True
            sage: D = FiniteField(3^20,'d')
            sage: D == loads(dumps(D))
            True
        """
        return self._factory_data[0].reduce_data(self)

def unpickle_FiniteField_ext(_type, order, variable_name, modulus, kwargs):
    r"""
    Used to unpickle extensions of finite fields. Now superseded (hence no
    doctest), but kept around for backward compatibility.

    EXAMPLE::

        sage: # not tested
    """
    return _type(order, variable_name, modulus, **kwargs)

def unpickle_FiniteField_prm(_type, order, variable_name, kwargs):
    r"""
    Used to unpickle finite prime fields. Now superseded (hence no doctest),
    but kept around for backward compatibility.

    EXAMPLE::

        sage: # not tested
    """
    return _type(order, variable_name, **kwargs)


def is_FiniteField(x):
    """
    Return True if x is of type finite field, and False otherwise.

    EXAMPLES::

        sage: from sage.rings.ring import is_FiniteField
        sage: is_FiniteField(GF(9,'a'))
        True
        sage: is_FiniteField(GF(next_prime(10^10)))
        True

    Note that the integers modulo n are not of type finite field,
    so this function returns False::

        sage: is_FiniteField(Integers(7))
        False
    """
    return IS_INSTANCE(x, FiniteField)

cdef class Algebra(Ring):
    """
    Generic algebra
    """

    # this is a total no-op!
    #def __init__(self, base_ring, names=None, normalize=True):
    #    ParentWithGens.__init__(self, base_ring, names=names, normalize=normalize)

    def characteristic(self):
        r"""
        Return the characteristic of this algebra, which is the same
        as the characteristic of its base ring.

        See objects with the ``base_ring`` attribute for additional examples.
        Here are some examples that explicitly use the ``Algebra`` class.

        EXAMPLES::

            sage: A = Algebra(ZZ); A
            <type 'sage.rings.ring.Algebra'>
            sage: A.characteristic()
            0
            sage: A = Algebra(GF(7^3, 'a'))
            sage: A.characteristic()
            7
        """
        return self.base_ring().characteristic()


cdef class CommutativeAlgebra(CommutativeRing):
    """
    Generic commutative algebra
    """
    def __init__(self, base_ring, names=None, normalize=True):
        r"""
        Standard init function. This just checks that the base is a commutative
        ring and then passes the buck.

        EXAMPLE::

            sage: sage.rings.ring.CommutativeAlgebra(QQ) # indirect doctest
            <type 'sage.rings.ring.CommutativeAlgebra'>

            sage: sage.rings.ring.CommutativeAlgebra(QuaternionAlgebra(QQ,-1,-1)) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: base ring must be a commutative ring
        """
        if not isinstance(base_ring, CommutativeRing):
            raise TypeError, "base ring must be a commutative ring"
        ParentWithGens.__init__(self, base_ring, names=names, normalize=normalize)

    def is_commutative(self):
        """
        Return True since this algebra is commutative.

        EXAMPLES:

        Any commutative ring is a commutative algebra over itself::

            sage: A = sage.rings.ring.CommutativeAlgebra
            sage: A(ZZ).is_commutative()
            True
            sage: A(QQ).is_commutative()
            True

        Trying to create a commutative algebra over a non-commutative ring
        will result in a ``TypeError``.
        """
        return True


def is_Ring(x):
    """
    Return True if x is a ring.

    EXAMPLES::

        sage: from sage.rings.ring import is_Ring
        sage: is_Ring(ZZ)
        True
    """
    return isinstance(x, Ring)


from sage.structure.parent_gens import _certify_names

def gen_name(x, name_chr):
    r"""
    Used to find a name for a generator when rings are created using the
    ``__getitem__`` syntax, e.g. ``ZZ['x']``. If x is a symbolic variable,
    return the name of x; if x is the symbolic square root of a positive
    integer d, return "sqrtd"; else, return a letter of the alphabet and
    increment a counter to avoid that letter being used again.

    EXAMPLES::

        sage: from sage.rings.ring import gen_name
        sage: gen_name(sqrt(5), 1)
        ('sqrt5', 1)
        sage: gen_name(sqrt(-17), 88)
        ('X', 89)
        sage: gen_name(x, 1)
        ('x', 1)
    """
    from sage.symbolic.ring import is_SymbolicVariable
    if is_SymbolicVariable(x):
        return repr(x), name_chr
    name = str(x)
    m = re.match('^sqrt\((\d+)\)$', name)
    if m:
        name = "sqrt%s" % m.groups()[0]
    try:
        _certify_names([name])
    except ValueError, msg:
        name = chr(name_chr)
        name_chr += 1
    return name, name_chr
