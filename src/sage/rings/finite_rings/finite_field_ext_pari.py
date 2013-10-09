"""
Finite Extension Fields implemented via PARI.

AUTHORS:

- William Stein: initial version

- Jeroen Demeyer (2010-12-16): fix formatting of docstrings (:trac:`10487`)
"""
#*****************************************************************************
#       Copyright (C) 2005,2007 William Stein <wstein@gmail.com>
#       Copyright (C) 2010 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import sage.rings.polynomial.polynomial_element as polynomial_element
import sage.rings.polynomial.multi_polynomial_element as multi_polynomial_element

import sage.rings.integer as integer
import sage.rings.rational as rational

import sage.libs.pari.all as pari

import element_ext_pari

from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteField_generic


import sage.interfaces.gap

class FiniteField_ext_pari(FiniteField_generic):
    r"""
    Finite Field of order `q`, where `q` is a prime power (not a prime),
    implemented using PARI ``POLMOD``. This implementation is the default
    implementation for `q \geq 2^{16}`.

    INPUT:

    - ``q`` -- integer, size of the finite field, not prime

    - ``name`` -- variable used for printing element of the finite
      field.  Also, two finite fields are considered equal if they
      have the same variable name, and not otherwise.

    - ``modulus`` -- you may provide a polynomial to use for reduction or
      a string:

      - ``'conway'`` -- force the use of a Conway polynomial, will
        raise a ``RuntimeError`` if none is found in the database
      - ``'random'`` -- use a random irreducible polynomial
      - ``'default'`` -- a Conway polynomial is used if found. Otherwise
        a random polynomial is used

    OUTPUT:

    A finite field of order `q` with the given variable name

    EXAMPLES::

        sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
        sage: k = FiniteField_ext_pari(9, 'a')
        sage: k
        Finite Field in a of size 3^2
        sage: k.is_field()
        True
        sage: k.characteristic()
        3
        sage: a = k.gen()
        sage: a
        a
        sage: a.parent()
        Finite Field in a of size 3^2
        sage: a.charpoly('x')
        x^2 + 2*x + 2
        sage: [a^i for i in range(8)]
        [1, a, a + 1, 2*a + 1, 2, 2*a, 2*a + 2, a + 2]

    Fields can be coerced into sets or list and iterated over::

        sage: list(k)
        [0, 1, 2, a, a + 1, a + 2, 2*a, 2*a + 1, 2*a + 2]

    The following is a native Python set::

        sage: set(k)
        set([0, 1, 2, a, a + 1, a + 2, 2*a, 2*a + 1, 2*a + 2])

    And the following is a Sage set::

        sage: Set(k)
        {0, 1, 2, a, a + 1, a + 2, 2*a, 2*a + 1, 2*a + 2}

        We can also make a list via comprehension:
        sage: [x for x in k]
        [0, 1, 2, a, a + 1, a + 2, 2*a, 2*a + 1, 2*a + 2]

    Next we compute with the finite field of order 16, where
    the name is named ``b``::

        sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
        sage: k16 = FiniteField_ext_pari(16, "b")
        sage: z = k16.gen()
        sage: z
        b
        sage: z.charpoly('x')
        x^4 + x + 1
        sage: k16.is_field()
        True
        sage: k16.characteristic()
        2
        sage: z.multiplicative_order()
        15

    Of course one can also make prime finite fields::

        sage: k = FiniteField(7)

    Note that the generator is 1::

        sage: k.gen()
        1
        sage: k.gen().multiplicative_order()
        1

    Prime finite fields are implemented elsewhere, they cannot be
    constructed using :class:`FiniteField_ext_pari`::

        sage: k = FiniteField_ext_pari(7, 'a')
        Traceback (most recent call last):
        ...
        ValueError: The size of the finite field must not be prime.

    Illustration of dumping and loading::

        sage: K = FiniteField(7)
        sage: loads(K.dumps()) == K
        True
        sage: K = FiniteField(7^10, 'b', impl='pari_mod')
        sage: loads(K.dumps()) == K
        True
        sage: K = FiniteField(7^10, 'a', impl='pari_mod')
        sage: loads(K.dumps()) == K
        True

    In this example `K` is large enough that Conway polynomials are not
    used.  Note that when the field is dumped the defining polynomial `f`
    is also dumped.  Since `f` is determined by a random algorithm, it's
    important that `f` is dumped as part of `K`.  If you quit Sage and
    restart and remake a finite field of the same order (and the order is
    large enough so that there is no Conway polynomial), then defining
    polynomial is probably different.  However, if you load a previously
    saved field, that will have the same defining polynomial. ::

        sage: K = GF(10007^10, 'a', impl='pari_mod')
        sage: loads(K.dumps()) == K
        True

    .. NOTE::

        We do NOT yet define natural consistent inclusion maps
        between different finite fields.
    """
    def __init__(self, q, name, modulus=None):
        """
        Create finite field of order `q` with variable printed as name.

        EXAMPLES::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(9, 'a'); k
            Finite Field in a of size 3^2
        """
        if element_ext_pari.dynamic_FiniteField_ext_pariElement is None: element_ext_pari._late_import()
        from constructor import FiniteField as GF
        q = integer.Integer(q)
        if q < 2:
            raise ArithmeticError, "q must be a prime power"
        from sage.structure.proof.all import arithmetic
        proof = arithmetic()
        if proof:
            F = q.factor()
        else:
            from sage.rings.arith import is_pseudoprime_small_power
            F = is_pseudoprime_small_power(q, get_data=True)
        if len(F) != 1:
            raise ArithmeticError, "q must be a prime power"

        if F[0][1] > 1:
            base_ring = GF(F[0][0])
        else:
            raise ValueError, "The size of the finite field must not be prime."
            #base_ring = self

        FiniteField_generic.__init__(self, base_ring, name, normalize=True)

        self._kwargs = {}
        self.__char = F[0][0]
        self.__pari_one = pari.pari(1).Mod(self.__char)
        self.__degree = integer.Integer(F[0][1])
        self.__order = q
        self.__is_field = True

        if modulus is None or modulus == "default":
            from conway_polynomials import exists_conway_polynomial
            if exists_conway_polynomial(self.__char, self.__degree):
                modulus = "conway"
            else:
                modulus = "random"

        if isinstance(modulus,str):
            if modulus == "conway":
                from conway_polynomials import conway_polynomial
                modulus = conway_polynomial(self.__char, self.__degree)
            elif modulus == "random":
                # The following is fast/deterministic, but has serious problems since
                # it crashes on 64-bit machines, and I can't figure out why:
                #     self.__pari_modulus = pari.pari.finitefield_init(self.__char, self.__degree, self.variable_name())
                # So instead we iterate through random polys until we find an irreducible one.

                R = GF(self.__char)['x']
                while True:
                    modulus = R.random_element(self.__degree)
                    modulus = modulus.monic()
                    if modulus.degree() == self.__degree and modulus.is_irreducible():
                        break
            else:
                raise ValueError("Modulus parameter not understood")

        elif isinstance(modulus, (list, tuple)):
            modulus = GF(self.__char)['x'](modulus)
        elif sage.rings.polynomial.polynomial_element.is_Polynomial(modulus):
            if modulus.parent() is not base_ring:
                modulus = modulus.change_ring(base_ring)
        else:
            raise ValueError("Modulus parameter not understood")

        self.__modulus = modulus
        f = pari.pari(str(modulus))
        self.__pari_modulus = f.subst(modulus.parent().variable_name(), 'a') * self.__pari_one
        self.__gen = element_ext_pari.FiniteField_ext_pariElement(self, pari.pari('a'))

        self._zero_element = self._element_constructor_(0)
        self._one_element = self._element_constructor_(1)

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: k.<b> = GF(5^20, impl='pari_mod'); type(k)
            <class 'sage.rings.finite_rings.finite_field_ext_pari.FiniteField_ext_pari_with_category'>
            sage: k is loads(dumps(k))
            True
        """
        return self._factory_data[0].reduce_data(self)

    def __cmp__(self, other):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: k = GF(7^20, 'a', impl='pari_mod')
            sage: k == loads(dumps(k))
            True
        """
        if not isinstance(other, FiniteField_ext_pari):
            return cmp(type(self), type(other))
        return cmp((self.__order, self.variable_name(), self.__modulus), (other.__order, other.variable_name(), other.__modulus))

    def __richcmp__(left, right, op):
        r"""
        Compare ``left`` with ``right``.

        EXAMPLE::

            sage: k.<a> = GF(2^17, impl='pari_mod')
            sage: j.<b> = GF(2^18, impl='pari_mod')
            sage: k == j
            False

            sage: GF(2^17, 'a', impl='pari_mod') == copy(GF(2^17, 'a', impl='pari_mod'))
            True
        """
        return left._richcmp_helper(right, op)

    def _pari_one(self):
        r"""
        The PARI object ``Mod(1,p)``.  This is implementation specific
        and should be ignored by users.

        EXAMPLES::

            sage: k = GF(7^20, 'a', impl='pari_mod')
            sage: k._pari_one()
            Mod(1, 7)
        """
        return self.__pari_one

    def _pari_modulus(self):
        """
        The polynomial mod `p` that defines the finite field, as a PARI
        object.  This is implementation specific, and some finite fields
        might not be implemented using PARI, so you should avoid using
        this function.

        OUTPUT:

        - ``gen`` -- a PARI polynomial gen

        EXAMPLES::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: FiniteField_ext_pari(19^2, 'a')._pari_modulus()
            Mod(1, 19)*a^2 + Mod(18, 19)*a + Mod(2, 19)

            sage: FiniteField_ext_pari(13^3, 'a')._pari_modulus()
            Mod(1, 13)*a^3 + Mod(2, 13)*a + Mod(11, 13)

        Note that the PARI modulus is always in terms of a, even if
        the field variable isn't.  This is because the specific choice
        of variable name has meaning in PARI, i.e., it can't be
        arbitrary. ::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: FiniteField_ext_pari(2^4, "b")._pari_modulus()
            Mod(1, 2)*a^4 + Mod(1, 2)*a + Mod(1, 2)
        """
        return self.__pari_modulus

    def gen(self, n=0):
        """
        Return a generator of the finite field.

        This generator is a root of the defining polynomial of the finite
        field, and might differ between different runs or different
        architectures.

        .. WARNING::

            This generator is not guaranteed to be a generator
            for the multiplicative group.  To obtain the latter, use
            :meth:`~sage.rings.finite_rings.finite_field_base.FiniteFields.multiplicative_generator()`.

        INPUT:

        - ``n`` -- ignored

        OUTPUT:

        Field generator of finite field

        EXAMPLES::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: FiniteField_ext_pari(2^4, "b").gen()
            b
            sage: k = FiniteField_ext_pari(3^4, "alpha")
            sage: a = k.gen()
            sage: a
            alpha
            sage: a^4
            alpha^3 + 1

        """
        return self.__gen

    def characteristic(self):
        """
        Returns the characteristic of the finite field, which is a
        prime number.

        EXAMPLES::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(3^4, 'a')
            sage: k.characteristic()
            3
        """
        return self.__char

    def modulus(self):
        r"""
        Return the minimal polynomial of the generator of ``self`` in
        ``self.polynomial_ring('x')``.

        EXAMPLES::

            sage: F.<a> = GF(7^20, 'a', impl='pari_mod')
            sage: f = F.modulus(); f
            x^20 + x^12 + 6*x^11 + 2*x^10 + 5*x^9 + 2*x^8 + 3*x^7 + x^6 + 3*x^5 + 3*x^3 + x + 3

            sage: f(a)
            0
        """
        return self.__modulus

    def degree(self):
        """
        Returns the degree of the finite field, which is a positive
        integer.

        EXAMPLES::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: FiniteField_ext_pari(3^20, 'a').degree()
            20
        """
        return self.__degree

    def _element_constructor_(self, x):
        r"""
        Coerce ``x`` into the finite field.

        INPUT:

        - ``x`` -- object

        OUTPUT:

        If possible, makes a finite field element from ``x``.

        EXAMPLES::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(3^4, 'a')
            sage: b = k(5) # indirect doctest
            sage: b.parent()
            Finite Field in a of size 3^4
            sage: a = k.gen()
            sage: k(a + 2)
            a + 2

        Univariate polynomials coerce into finite fields by evaluating
        the polynomial at the field's generator::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: R.<x> = QQ[]
            sage: k, a = FiniteField_ext_pari(5^2, 'a').objgen()
            sage: k(R(2/3))
            4
            sage: k(x^2)
            a + 3
            sage: R.<x> = GF(5)[]
            sage: k(x^3-2*x+1)
            2*a + 4

            sage: x = polygen(QQ)
            sage: k(x^25)
            a

            sage: Q, q = FiniteField_ext_pari(5^7, 'q').objgen()
            sage: L = GF(5)
            sage: LL.<xx> = L[]
            sage: Q(xx^2 + 2*xx + 4)
            q^2 + 2*q + 4

        Multivariate polynomials only coerce if constant::

            sage: R = k['x,y,z']; R
            Multivariate Polynomial Ring in x, y, z over Finite Field in a of size 5^2
            sage: k(R(2))
            2
            sage: R = QQ['x,y,z']
            sage: k(R(1/5))
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce

        Gap elements can also be coerced into finite fields::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: F = FiniteField_ext_pari(8, 'a')
            sage: a = F.multiplicative_generator(); a
            a
            sage: b = gap(a^3); b
            Z(2^3)^3
            sage: F(b)
            a + 1
            sage: a^3
            a + 1

            sage: a = GF(13)(gap('0*Z(13)')); a
            0
            sage: a.parent()
            Finite Field of size 13

            sage: F = GF(16, 'a', impl='pari_mod')
            sage: F(gap('Z(16)^3'))
            a^3
            sage: F(gap('Z(16)^2'))
            a^2

        You can also call a finite extension field with a string
        to produce an element of that field, like this::

            sage: k = GF(2^8, 'a', impl='pari_mod')
            sage: k('a^200')
            a^4 + a^3 + a^2

        This is especially useful for fast conversions from Singular etc.
        to ``FiniteField_ext_pariElements``.

        AUTHORS:

        - David Joyner (2005-11)
        - Martin Albrecht (2006-01-23)
        - Martin Albrecht (2006-03-06): added coercion from string
        """
        if isinstance(x, element_ext_pari.FiniteField_ext_pariElement):
            if x.parent() is self:
                return x
            elif x.parent() == self:
                # canonically isomorphic finite fields
                return element_ext_pari.FiniteField_ext_pariElement(self, x)
            else:
                # This is where we *would* do coercion from one finite field to another...
                raise TypeError, "no coercion defined"

        elif sage.interfaces.gap.is_GapElement(x):
            from sage.interfaces.gap import gfq_gap_to_sage
            try:
                return gfq_gap_to_sage(x, self)
            except (ValueError, IndexError, TypeError):
                raise TypeError, "no coercion defined"

        if isinstance(x, (int, long, integer.Integer, rational.Rational,
                          pari.pari_gen, list)):

            return element_ext_pari.FiniteField_ext_pariElement(self, x)

        elif isinstance(x, multi_polynomial_element.MPolynomial):
            if x.is_constant():
                return self(x.constant_coefficient())
            else:
                raise TypeError, "no coercion defined"

        elif isinstance(x, polynomial_element.Polynomial):
            if x.is_constant():
                return self(x.constant_coefficient())
            else:
                return x.change_ring(self)(self.gen())

        elif isinstance(x, str):
            x = x.replace(self.variable_name(),'a')
            x = pari.pari(x)
            t = x.type()
            if t == 't_POL':
                if (x.variable() == 'a' \
                    and x.polcoeff(0).type()[2] == 'I'): #t_INT and t_INTMOD
                    return self(x)
            if t[2] == 'I': #t_INT and t_INTMOD
                return self(x)
            raise TypeError, "string element does not match this finite field"

        try:
            if x.parent() == self.vector_space():
                x = pari.pari('+'.join(['%s*a^%s'%(x[i], i) for i in range(self.degree())]))
                return element_ext_pari.FiniteField_ext_pariElement(self, x)
        except AttributeError:
            pass
        try:
            return element_ext_pari.FiniteField_ext_pariElement(self, integer.Integer(x))
        except TypeError, msg:
            raise TypeError, "%s\nno coercion defined"%msg

    def _coerce_map_from_(self, R):
        r"""
        Canonical coercion to ``self``.

        EXAMPLES::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: FiniteField_ext_pari(4,'a')._coerce_(GF(2)(1)) # indirect doctest
            1
            sage: k = FiniteField_ext_pari(4,'a')
            sage: k._coerce_(k.0)
            a
            sage: FiniteField_ext_pari(4,'a')._coerce_(3)
            1
            sage: FiniteField_ext_pari(4,'a')._coerce_(2/3)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Finite Field in a of size 2^2
            sage: FiniteField_ext_pari(8,'a')._coerce_(FiniteField_ext_pari(4,'a').0)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Finite Field in a of size 2^2 to Finite Field in a of size 2^3
            sage: FiniteField_ext_pari(16,'a')._coerce_(FiniteField_ext_pari(4,'a').0)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Finite Field in a of size 2^2 to Finite Field in a of size 2^4
            sage: k = FiniteField_ext_pari(8,'a')
            sage: k._coerce_(FiniteField(7,'a')(2))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Finite Field of size 7 to Finite Field in a of size 2^3
        """
        from sage.rings.integer_ring import ZZ
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
        if R is int or R is long or R is ZZ:
            return True
        if isinstance(R, FiniteField_ext_pari):
            if R is self:
                return True
            if R.characteristic() == self.characteristic():
                if R.degree() == 1:
                    return True
                elif self.degree() % R.degree() == 0:
                    # TODO: This is where we *would* do coercion from one nontrivial finite field to another...
                    return False
        from sage.rings.residue_field import ResidueField_generic
        if isinstance(R, IntegerModRing_generic) and R.characteristic() == self.characteristic() and not isinstance(R, ResidueField_generic):
            return True

    def __len__(self):
        """
        The number of elements of the finite field.

        EXAMPLES::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(2^10, 'a')
            sage: k
            Finite Field in a of size 2^10
            sage: len(k)
            1024
        """
        return self.__order

    def order(self):
        """
        The number of elements of the finite field.

        EXAMPLES::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(2^10,'a')
            sage: k
            Finite Field in a of size 2^10
            sage: k.order()
            1024
        """
        return self.__order

    def polynomial(self, name=None):
        """
        Return the irreducible characteristic polynomial of the
        generator of this finite field, i.e., the polynomial `f(x)` so
        elements of the finite field as elements modulo `f`.

        EXAMPLES::

            sage: k = FiniteField(17)
            sage: k.polynomial('x')
            x
            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(9,'a')
            sage: k.polynomial('x')
            x^2 + 2*x + 2
        """
        if name is None:
            name = self.variable_name()
        try:
            return self.__polynomial[name]
        except (AttributeError, KeyError):
            from constructor import FiniteField as GF
            R = GF(self.characteristic())[name]
            f = R(self._pari_modulus())
            try:
                self.__polynomial[name] = f
            except (KeyError, AttributeError):
                self.__polynomial = {}
                self.__polynomial[name] = f
            return f

    def __hash__(self):
        """
        Return the hash of this field.

        EXAMPLES::

            sage: {GF(9, 'a', impl='pari_mod'): 1} # indirect doctest
            {Finite Field in a of size 3^2: 1}
            sage: {GF(9, 'b', impl='pari_mod'): 1} # indirect doctest
            {Finite Field in b of size 3^2: 1}
        """
        try:
            return self.__hash
        except AttributeError:
            self.__hash = hash((self.__order, self.variable_name(), self.__modulus))
            return self.__hash
