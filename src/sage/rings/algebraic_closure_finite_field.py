r"""
Algebraic closures of finite fields

Let `\Bold{F}` be a finite field, and let `\overline{\Bold{F}}` be an
algebraic closure of `\Bold{F}`; this is unique up to (non-canonical)
isomorphism.  For every `n\ge 1`, there is a unique subfield
`\Bold{F}_n` of `\overline{\Bold{F}}` such that
`\Bold{F}\subset\Bold{F}_n` and `[\Bold{F}_n:\Bold{F}]=n`.

In Sage, algebraic closures of finite fields are implemented using
compatible systems of finite fields.  The resulting Sage object keeps
track of a finite lattice of the subfields `\Bold{F}_n` and the
embeddings between them.  This lattice is extended as necessary.

The Sage class corresponding to `\overline{\Bold{F}}` can be
constructed from the finite field `\Bold{F}` by using the
:meth:`~sage.rings.finite_rings.finite_field_base.FiniteField.algebraic_closure`
method.

The Sage class for elements of `\overline{\Bold{F}}` is
:class:`AlgebraicClosureFiniteFieldElement`.  Such an element is
represented as an element of one of the `\Bold{F}_n`.  This means that
each element `x\in\Bold{F}` has infinitely many different
representations, one for each `n` such that `x` is in `\Bold{F}_n`.

.. NOTE::

    Only prime finite fields are currently accepted as base fields for
    algebraic closures.  To obtain an algebraic closure of a non-prime
    finite field `\Bold{F}`, take an algebraic closure of the prime
    field of `\Bold{F}` and embed `\Bold{F}` into this.

    Algebraic closures of finite fields are currently implemented
    using (pseudo-)Conway polynomials; see
    :class:`AlgebraicClosureFiniteField_pseudo_conway` and the module
    :mod:`~sage.rings.finite_rings.conway_polynomials`.  Other
    implementations may be added by creating appropriate subclasses of
    :class:`AlgebraicClosureFiniteField_generic`.

    In the current implementation, algebraic closures do not satisfy
    the unique parent condition.  Moreover, there is no coercion map
    between different algebraic closures of the same finite field.
    There is a conceptual reason for this, namely that the definition
    of pseudo-Conway polynomials only determines an algebraic closure
    up to *non-unique* isomorphism.  This means in particular that
    different algebraic closures, and their respective elements, never
    compare equal.

AUTHORS:

- Peter Bruin (August 2013): initial version

- Vincent Delecroix (November 2013): additional methods

"""

from sage.misc.abstract_method import abstract_method
from sage.misc.fast_methods import WithEqualityById

from sage.rings.finite_rings.element_base import is_FiniteFieldElement
from sage.rings.finite_rings.finite_field_base import is_FiniteField
from sage.rings.ring import Field
from sage.structure.element import FieldElement

class AlgebraicClosureFiniteFieldElement(FieldElement):
    """
    Element of an algebraic closure of a finite field.

    EXAMPLES::

        sage: F = GF(3).algebraic_closure()
        sage: F.gen(2)
        z2
        sage: type(F.gen(2))
        <class 'sage.rings.algebraic_closure_finite_field.AlgebraicClosureFiniteField_pseudo_conway_with_category.element_class'>

    """
    def __init__(self, parent, value):
        """
        TEST::

            sage: F = GF(3).algebraic_closure()
            sage: TestSuite(F.gen(2)).run(skip=['_test_pickling'])

        .. NOTE::

            The ``_test_pickling`` test has to be skipped because
            there is no coercion map between the parents of ``x``
            and ``loads(dumps(x))``.

        """
        if is_FiniteFieldElement(value):
            n = value.parent().degree()
        else:
            from sage.rings.integer import Integer
            n = Integer(1)
        self._value = parent._subfield(n).coerce(value)
        self._level = n
        FieldElement.__init__(self, parent)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure()
            sage: F._repr_()
            'Algebraic closure of Finite Field of size 3'

        """
        return self._value._repr_()

    def __cmp__(self, right):
        """
        Compare ``self`` with ``right``.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure()
            sage: F.gen(2) == F.gen(3)
            False

        """
        x, y = self.parent()._to_common_subfield(self, right)
        return cmp(x, y)

    def _add_(self, right):
        """
        Return ``self`` + ``right``.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure()
            sage: F.gen(2) + F.gen(3)
            z6^5 + 2*z6^4 + 2*z6^3 + z6^2 + 2*z6 + 1

        """
        F = self.parent()
        x, y = F._to_common_subfield(self, right)
        return self.__class__(F, x + y)

    def _sub_(self, right):
        """
        Return ``self`` - ``right``.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure()
            sage: F.gen(2) - F.gen(3)
            z6^4 + 2*z6^3 + z6^2 + 2*z6

        """
        F = self.parent()
        x, y = F._to_common_subfield(self, right)
        return self.__class__(F, x - y)

    def _mul_(self, right):
        """
        Return ``self`` * ``right``.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure()
            sage: F.gen(2) * F.gen(3)
            z6^5 + 2*z6^4 + z6^2 + 2

        """
        F = self.parent()
        x, y = F._to_common_subfield(self, right)
        return self.__class__(F, x * y)

    def _div_(self, right):
        """
        Return ``self`` / ``right``.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure()
            sage: F.gen(2) / F.gen(3)
            z6^5 + 2*z6^4 + z6^3 + 1

        """
        F = self.parent()
        x, y = F._to_common_subfield(self, right)
        return self.__class__(F, x / y)

    def change_level(self, n):
        """
        Return a representation of ``self`` as an element of the
        subfield of degree `n` of the parent, if possible.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure()
            sage: z = F.gen(4)
            sage: (z^10).change_level(6)
            2*z6^5 + 2*z6^3 + z6^2 + 2*z6 + 2
            sage: z.change_level(6)
            Traceback (most recent call last):
            ...
            ValueError: z4 is not in the image of Ring morphism:
              From: Finite Field in z2 of size 3^2
              To:   Finite Field in z4 of size 3^4
              Defn: z2 |--> 2*z4^3 + 2*z4^2 + 1

            sage: a = F(1).change_level(3); a
            1
            sage: a.change_level(2)
            1
            sage: F.gen(3).change_level(1)
            Traceback (most recent call last):
            ...
            TypeError: not in prime subfield

        """
        F = self.parent()
        l = self._level
        m = l.gcd(n)
        xl = self._value
        if m == 1:
            # Sections of canonical coercion maps from F_p to other
            # finite fields of characteristic p are not implemented.
            xm = F.base_ring()(xl)
        else:
            xm = F.inclusion(m, l).section()(xl)
        xn = F.inclusion(m, n)(xm)
        return self.__class__(F, xn)

    def _latex_(self):
        """
        Return a LaTeX representation of ``self``.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure()
            sage: s = F.gen(1) + F.gen(2) + F.gen(3)
            sage: s
            z6^5 + 2*z6^4 + 2*z6^3 + z6^2 + 2*z6 + 2
            sage: latex(s)
            z_{6}^{5} + 2 z_{6}^{4} + 2 z_{6}^{3} + z_{6}^{2} + 2 z_{6} + 2

        """
        return self._value._latex_()

    def minpoly(self):
        """
        Return the minimal polynomial of ``self`` over the prime
        field.

        EXAMPLES::

            sage: F = GF(11).algebraic_closure()
            sage: F.gen(3).minpoly()
            x^3 + 2*x + 9

        """
        return self._value.minpoly()

    minimal_polynomial = minpoly

    def is_square(self):
        """
        Return ``True`` if ``self`` is a square.

        This always returns ``True``.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure()
            sage: F.gen(2).is_square()
            True

        """
        return True

    def sqrt(self):
        """
        Return a square root of ``self``.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure()
            sage: F.gen(2).sqrt()
            z4^3 + z4 + 1

        """
        F = self.parent()
        x = self._value
        if x.is_square():
            return self.__class__(F, x.sqrt(extend=False))
        else:
            l = self._level
            x = F.inclusion(l, 2*l)(x)
            return self.__class__(F, x.sqrt(extend=False))

    def nth_root(self, n):
        """
        Return an `n`-th root of ``self``.

        EXAMPLES::

            sage: F = GF(5).algebraic_closure()
            sage: t = F.gen(2) + 1
            sage: s = t.nth_root(15); s
            4*z6^5 + 3*z6^4 + 2*z6^3 + 2*z6^2 + 4
            sage: s**15 == t
            True

        .. TODO::

            This function could probably be made faster.

        """
        from sage.rings.integer import Integer
        F = self.parent()
        x = self._value
        n = Integer(n)
        l = self._level
        # In order to be smart we look for the smallest subfield that
        # actually contains the root.
        for d in n.divisors():
            xx = F.inclusion(l, d*l)(x)
            try:
                y = xx.nth_root(n, extend=False)
            except ValueError:
                continue
            return self.__class__(F, y)

        raise AssertionError('cannot find n-th root in algebraic closure of finite field')

    def multiplicative_order(self):
        """
        Return the multiplicative order of ``self``.

        EXAMPLES::

            sage: K = GF(7).algebraic_closure()
            sage: K.gen(5).multiplicative_order()
            16806
            sage: (K.gen(1) + K.gen(2) + K.gen(3)).multiplicative_order()
            7353

        """
        return self._value.multiplicative_order()

    def pth_power(self, k=1):
        """
        Return the `p^k`-th power of ``self``, where `p` is the
        characteristic of ``self.parent()``.

        EXAMPLES::

            sage: K = GF(13).algebraic_closure('t')
            sage: t3 = K.gen(3)
            sage: s = 1 + t3 + t3**2
            sage: s.pth_power()
            10*t3^2 + 6*t3
            sage: s.pth_power(2)
            2*t3^2 + 6*t3 + 11
            sage: s.pth_power(3)
            t3^2 + t3 + 1
            sage: s.pth_power(3).parent() is K
            True

        """
        return self.__class__(self.parent(), self._value.pth_power(k))

    def pth_root(self, k=1):
        """
        Return the unique `p^k`-th root of ``self``, where `p` is the
        characteristic of ``self.parent()``.

        EXAMPLES::

            sage: K = GF(13).algebraic_closure('t')
            sage: t3 = K.gen(3)
            sage: s = 1 + t3 + t3**2
            sage: s.pth_root()
            2*t3^2 + 6*t3 + 11
            sage: s.pth_root(2)
            10*t3^2 + 6*t3
            sage: s.pth_root(3)
            t3^2 + t3 + 1
            sage: s.pth_root(2).parent() is K
            True

        """
        return self.__class__(self.parent(), self._value.pth_root(k))

    def as_finite_field_element(self, minimal=False):
        """
        Return ``self`` as a finite field element.

        INPUT:

        - ``minimal`` -- boolean (default: ``False``).  If ``True``,
          always return the smallest subfield containing ``self``.

        OUTPUT:

        - a triple (``field``, ``element``, ``morphism``) where
          ``field`` is a finite field, ``element`` an element of
          ``field`` and ``morphism`` a morphism from ``field`` to
          ``self.parent()``.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure('t')
            sage: t = F.gen(5)
            sage: t.as_finite_field_element()
            (Finite Field in t5 of size 3^5,
             t5,
             Ring morphism:
              From: Finite Field in t5 of size 3^5
              To:   Algebraic closure of Finite Field of size 3
              Defn: t5 |--> t5)

        By default, ``field`` is not necessarily minimal.  We can
        force it to be minimal using the ``minimal`` option::

            sage: s = t + 1 - t
            sage: s.as_finite_field_element()[0]
            Finite Field in t5 of size 3^5
            sage: s.as_finite_field_element(minimal=True)[0]
            Finite Field of size 3

        """
        if not minimal:
            l = self._level
        else:
            l = self._value.minpoly().degree()

        F, phi = self.parent().subfield(l)
        return (F, F(self._value), phi)


class AlgebraicClosureFiniteField_generic(Field):
    """
    Algebraic closure of a finite field.

    """
    def __init__(self, base_ring, name, category=None):
        """
        TEST::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_generic
            sage: F = AlgebraicClosureFiniteField_generic(GF(5), 'z')
            sage: F
            Algebraic closure of Finite Field of size 5

        """
        Field.__init__(self, base_ring=base_ring, names=name,
                       normalize=False, category=category)

    def __cmp__(self, other):
        """
        Compare ``self`` with ``other``.

        TEST::

            sage: F3 = GF(3).algebraic_closure()
            sage: F3 == F3
            True
            sage: F5 = GF(5).algebraic_closure()
            sage: F3 == F5
            False

        """
        if self is other:
            return 0
        c = cmp(type(self), type(other))
        if c != 0:
            return c
        return cmp((self.base_ring(), self.variable_name(), self.category()),
                   (other.base_ring(), other.variable_name(), other.category()))

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        This always returns ``+Infinity``.

        .. TODO::

            When :trac:`10963` is merged we should remove that method and set the
            category to infinite fields (i.e. ``Fields().Infinite()``).

        EXAMPLES::

            sage: F = GF(3).algebraic_closure()
            sage: F.cardinality()
            +Infinity

        """
        from sage.rings.infinity import Infinity
        return Infinity

    def is_finite(self):
        """
        Returns ``False`` as an algebraically closed field is always infinite.

        .. TODO::

            When :trac:`10963` is merged we should remove that method and set the
            category to infinite fields (i.e. ``Fields().Infinite()``).

        EXAMPLES::

            sage: GF(3).algebraic_closure().is_finite()
            False

        """
        return False

    def characteristic(self):
        """
        Return the characteristic of ``self``.

        EXAMPLES::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: p = next_prime(1000)
            sage: F = AlgebraicClosureFiniteField(GF(p), 'z')
            sage: F.characteristic() == p
            True

        """
        return self.base_ring().characteristic()

    Element = AlgebraicClosureFiniteFieldElement

    def _element_constructor_(self, x):
        """
        Construct an element of ``self``.

        TEST::

            sage: F = GF(5).algebraic_closure()
            sage: type(F(3))
            <class 'sage.rings.algebraic_closure_finite_field.AlgebraicClosureFiniteField_pseudo_conway_with_category.element_class'>

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: F1 = AlgebraicClosureFiniteField(GF(3), 'z')
            sage: F2 = AlgebraicClosureFiniteField(GF(3), 'z')
            sage: F1(F2.gen(1))
            Traceback (most recent call last):
            ...
            ValueError: no conversion defined between different algebraic closures

        """
        if isinstance(x, self.element_class):
            if x.parent() is not self:
                raise ValueError('no conversion defined between different algebraic closures')
            return x
        else:
            return self.element_class(self, x)

    def _coerce_map_from_(self, other):
        """
        Return ``True`` if elements of ``other`` can be coerced into
        ``self``.

        EXAMPLES::

            sage: F = GF(7).algebraic_closure()
            sage: F.has_coerce_map_from(Integers())
            True

        """
        if other is self:
            return True
        elif is_FiniteField(other) and self._subfield(other.degree()) is other:
            return True
        elif self._subfield(1).has_coerce_map_from(other):
            return True

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: F = AlgebraicClosureFiniteField(GF(5), 'z')
            sage: F._repr_()
            'Algebraic closure of Finite Field of size 5'

        """
        return 'Algebraic closure of %s' % self.base_ring()

    def _to_common_subfield(self, x, y):
        """
        Coerce `x` and `y` to a common subfield of ``self``.

        TEST::

            sage: F = GF(3).algebraic_closure()
            sage: x, y = F._to_common_subfield(F.gen(2), F.gen(3))
            sage: x.parent()
            Finite Field in z6 of size 3^6
            sage: y.parent()
            Finite Field in z6 of size 3^6

        """
        if x._level == y._level:
            return x._value, y._value
        n = x._level.lcm(y._level)
        mx = self.inclusion(x._level, n)
        my = self.inclusion(y._level, n)
        return mx(x._value), my(y._value)

    @abstract_method
    def _get_polynomial(self, n):
        """
        Return the polynomial defining the unique subfield of degree
        `n` of ``self``.

        This must be implemented by subclasses.

        EXAMPLES::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_generic
            sage: F = AlgebraicClosureFiniteField_generic(GF(5), 'z')
            sage: F._get_polynomial(1)
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method _get_polynomial at ...>
        """

    @abstract_method
    def _get_im_gen(self, m, n):
        """
        Return the image of ``self.gen(m)`` under the canonical
        inclusion into ``self.subfield(n)``.

        This must be implemented by subclasses.

        EXAMPLES::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_generic
            sage: F = AlgebraicClosureFiniteField_generic(GF(5), 'z')
            sage: F._get_im_gen(2, 4)
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method _get_im_gen at ...>

        """

    def _subfield(self, n):
        """
        Return the unique subfield of degree `n` of ``self``.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure()
            sage: F._subfield(4)
            Finite Field in z4 of size 3^4

        """
        if n == 1:
            return self.base_ring()
        else:
            from sage.rings.finite_rings.constructor import FiniteField
            return FiniteField(self.base_ring().cardinality() ** n,
                               name=self.variable_name() + str(n),
                               modulus=self._get_polynomial(n),
                               check_irreducible=False)

    def subfield(self, n):
        """
        Return the unique subfield of degree `n` of ``self``
        together with its canonical embedding into ``self``.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure()
            sage: F.subfield(1)
            (Finite Field of size 3,
             Ring morphism:
               From: Finite Field of size 3
               To:   Algebraic closure of Finite Field of size 3
               Defn: 1 |--> 1)
            sage: F.subfield(4)
            (Finite Field in z4 of size 3^4,
             Ring morphism:
               From: Finite Field in z4 of size 3^4
               To:   Algebraic closure of Finite Field of size 3
               Defn: z4 |--> z4)

        """
        Fn = self._subfield(n)
        return Fn, Fn.hom((self.gen(n),))

    def inclusion(self, m, n):
        """
        Return the canonical inclusion map from the subfield
        of degree `m` to the subfield of degree `n`.

        EXAMPLES::

            sage: F = GF(3).algebraic_closure()
            sage: F.inclusion(1, 2)
            Ring Coercion morphism:
              From: Finite Field of size 3
              To:   Finite Field in z2 of size 3^2
            sage: F.inclusion(2, 4)
            Ring morphism:
              From: Finite Field in z2 of size 3^2
              To:   Finite Field in z4 of size 3^4
              Defn: z2 |--> 2*z4^3 + 2*z4^2 + 1

        """
        if m == 1:
            return self.base_ring().hom(self._subfield(n))
        elif m.divides(n):
            return self._subfield(m).hom((self._get_im_gen(m, n),))
        else:
            raise ValueError("subfield of degree %s not contained in subfield of degree %s" % (m, n))

    def ngens(self):
        """
        Return the number of generators of ``self``, which is
        infinity.

        EXAMPLES::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: AlgebraicClosureFiniteField(GF(5), 'z').ngens()
            +Infinity

        """
        from sage.rings.infinity import Infinity
        return Infinity

    def gen(self, n):
        """
        Return the `n`-th generator of ``self``.

        EXAMPLES::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: F = AlgebraicClosureFiniteField(GF(5), 'z')
            sage: F.gen(2)
            z2

        """
        F = self._subfield(n)
        return self(F.gen())

    def gens(self):
        """
        Return a family of generators of ``self``.

        OUTPUT:

        - a :class:`~sage.sets.family.Family`, indexed by the positive
          integers, whose `n`-th element is ``self.gen(n)``.

        EXAMPLES::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: F = AlgebraicClosureFiniteField(GF(5), 'z')
            sage: g = F.gens()
            sage: g
            Lazy family (<lambda>(i))_{i in Positive integers}
            sage: g[3]
            z3

        """
        from sage.sets.family import Family
        from sage.sets.positive_integers import PositiveIntegers

        return Family(PositiveIntegers(), lambda n: self.gen(n))

    def _first_ngens(self, n):
        """
        Return the first `n` generators of ``self``.

        EXAMPLES::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: F = AlgebraicClosureFiniteField(GF(5), 'z')
            sage: F._first_ngens(3)
            (1, z2, z3)

        """
        return tuple([self.gen(i + 1) for i in xrange(n)])

    def algebraic_closure(self):
        """
        Return an algebraic closure of ``self``.

        This always returns ``self``.

        EXAMPLES::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: F = AlgebraicClosureFiniteField(GF(5), 'z')
            sage: F.algebraic_closure() is F
            True

        """
        return self

    def _an_element_(self):
        """
        TEST::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
            sage: F = AlgebraicClosureFiniteField(GF(5), 'w')
            sage: F.an_element()  # indirect doctest
            w2

        """
        return self.gen(2)

    def some_elements(self):
        r"""
        Return some elements of this field.

        EXAMPLES::

            sage: F = GF(7).algebraic_closure()
            sage: F.some_elements()
            (1, z2, z3 + 1)

        """
        return (self(1), self.gen(2), 1+self.gen(3))


class AlgebraicClosureFiniteField_pseudo_conway(AlgebraicClosureFiniteField_generic, WithEqualityById):
    """
    Algebraic closure of a finite field, constructed using
    pseudo-Conway polynomials.

    EXAMPLES::

        sage: F = GF(5).algebraic_closure(implementation='pseudo_conway')
        sage: F.cardinality()
        +Infinity
        sage: F.algebraic_closure() is F
        True
        sage: x = F(3).nth_root(12); x
        z4^3 + z4^2 + 4*z4
        sage: x**12
        3

    TESTS::

        sage: F3 = GF(3).algebraic_closure()
        sage: F3 == F3
        True
        sage: F5 = GF(5).algebraic_closure()
        sage: F3 == F5
        False

    """
    def __init__(self, base_ring, name, category=None, lattice=None, use_database=True):
        """
        INPUT:

        - ``base_ring`` -- the finite field of which to construct an
          algebraic closure.  Currently only prime fields are
          accepted.

        - ``name`` -- prefix to use for generators of the finite
          subfields.

        - ``category`` -- if provided, specifies the category in which
          this algebraic closure will be placed.

        - ``lattice`` -- :class:`~sage.rings.finite_rings.conway_polynomials.PseudoConwayPolynomialLattice`
          (default: None).  If provided, use this pseudo-Conway
          polynonomial lattice to construct an algebraic closure.

        - ``use_database`` -- boolean.  If True (default), use actual
          Conway polynomials whenever they are available in the
          database.  If False, always compute pseudo-Conway
          polynomials from scratch.

        TESTS::

            sage: F = GF(5).algebraic_closure(implementation='pseudo_conway')
            sage: print F.__class__.__name__
            AlgebraicClosureFiniteField_pseudo_conway_with_category
            sage: TestSuite(F).run(skip=['_test_elements', '_test_pickling'])

            sage: from sage.rings.finite_rings.conway_polynomials import PseudoConwayLattice
            sage: L = PseudoConwayLattice(11, use_database=False)
            sage: F = GF(7).algebraic_closure(lattice=L)
            Traceback (most recent call last):
            ...
            TypeError: lattice must be a pseudo-Conway lattice with characteristic 7
            sage: F = GF(11).algebraic_closure(lattice=L)
            sage: F.gen(2).minimal_polynomial()
            x^2 + 4*x + 2

            sage: F = GF(11).algebraic_closure(use_database=True)
            sage: F.gen(2).minimal_polynomial()
            x^2 + 7*x + 2

        .. NOTE::

            In the test suite, ``_test_pickling`` has to be skipped
            because ``F`` and ``loads(dumps(F))`` cannot consistently
            be made to compare equal, and ``_test_elements`` has to be
            skipped for the reason described in
            :meth:`AlgebraicClosureFiniteFieldElement.__init__`.

        """
        if not (is_FiniteField(base_ring) and base_ring.is_prime_field()):
            raise NotImplementedError('algebraic closures of finite fields are only implemented for prime fields')
        from sage.rings.finite_rings.conway_polynomials import PseudoConwayLattice
        p = base_ring.characteristic()
        if lattice is None:
            lattice = PseudoConwayLattice(p, use_database)
        elif not isinstance(lattice, PseudoConwayLattice) or lattice.p != p:
            raise TypeError('lattice must be a pseudo-Conway lattice with characteristic %s' % p)
        self._pseudo_conway_lattice = lattice
        AlgebraicClosureFiniteField_generic.__init__(self, base_ring, name, category)

    def _get_polynomial(self, n):
        """
        Return the defining polynomial of the unique subfield of
        degree `n` of ``self``.

        EXAMPLES::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_pseudo_conway
            sage: F = AlgebraicClosureFiniteField_pseudo_conway(GF(5), 'z')
            sage: F._get_polynomial(1)
            x + 3

        """
        return self._pseudo_conway_lattice.polynomial(n)

    def _get_im_gen(self, m, n):
        """
        Return the image of ``self.gen(m)`` under the canonical
        inclusion into ``self.subfield(n)``.

        EXAMPLES::

            sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField_pseudo_conway
            sage: F = AlgebraicClosureFiniteField_pseudo_conway(GF(5), 'z')
            sage: F._get_im_gen(2, 4)
            z4^3 + z4^2 + z4 + 3

        """
        p = self.characteristic()
        return self._subfield(n).gen() ** ((p**n - 1)//(p**m - 1))


def AlgebraicClosureFiniteField(base_ring, name, category=None, implementation=None, **kwds):
    """
    Construct an algebraic closure of a finite field.

    The recommended way to use this functionality is by calling the
    :meth:`~sage.rings.finite_rings.finite_field_base.FiniteField.algebraic_closure`
    method of the finite field.

    .. NOTE::

        Algebraic closures of finite fields in Sage do not have the
        unique representation property, because they are not
        determined up to unique isomorphism by their defining data.

    EXAMPLES::

        sage: from sage.rings.algebraic_closure_finite_field import AlgebraicClosureFiniteField
        sage: F = GF(2).algebraic_closure()
        sage: F1 = AlgebraicClosureFiniteField(GF(2), 'z')
        sage: F1 is F
        False

    In the pseudo-Conway implementation, non-identical instances never
    compare equal::

        sage: F1 == F
        False
        sage: loads(dumps(F)) == F
        False

    This is to ensure that the result of comparing two instances
    cannot change with time.

    """
    if category is None:
        from sage.categories.fields import Fields
        category = Fields()
    if implementation is None:
        implementation = 'pseudo_conway'

    if implementation == 'pseudo_conway':
        return AlgebraicClosureFiniteField_pseudo_conway(base_ring, name, category, **kwds)
    else:
        raise ValueError('unknown implementation for algebraic closure of finite field: %s'
                         % implementation)
