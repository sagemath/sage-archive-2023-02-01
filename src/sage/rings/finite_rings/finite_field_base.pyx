"""
Base Classes for Finite Fields

TESTS::

    sage: K.<a> = NumberField(x^2 + 1)
    sage: F = K.factor(3)[0][0].residue_field()
    sage: loads(dumps(F)) == F
    True
"""
include "sage/ext/stdsage.pxi"

from sage.structure.parent cimport Parent
from sage.misc.cachefunc import cached_method
from sage.misc.prandom import randrange
from sage.rings.finite_rings.homset import FiniteFieldHomset

cdef class FiniteFieldIterator:
    r"""
    An iterator over a finite field. This should only be used when the field
    is an extension of a smaller field which already has a separate iterator
    function.
    """
    cdef object iter
    cdef FiniteField parent

    def __init__(self,FiniteField parent):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = iter(FiniteField_ext_pari(9, 'a')) # indirect doctest
            sage: isinstance(k, sage.rings.finite_rings.finite_field_base.FiniteFieldIterator)
            True
        """
        self.parent = parent
        self.iter =iter(self.parent.vector_space())

    def __next__(self):
        r"""
        Return the next element in the iterator.

        EXAMPLE::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = iter(FiniteField_ext_pari(9, 'a'))
            sage: k.next() # indirect doctest
            0
        """
        return self.parent(self.iter.next())

    def __iter__(self):
        """
        Return ``self`` since this is an interator class.

        EXAMPLES::

            sage: K.<a> = GF(7^4)
            sage: K.list()[:7]
            [0, a, a^2, a^3, 2*a^2 + 3*a + 4, 2*a^3 + 3*a^2 + 4*a, 3*a^3 + a^2 + 6*a + 1]
            sage: K.<a> = GF(5^9)
            sage: for x in K:
            ...       if x == a+3: break
            ...       print x
            0
            1
            2
            3
            4
            a
            a + 1
            a + 2
        """
        return self

from sage.categories.finite_fields import FiniteFields
_FiniteFields = FiniteFields()
cdef class FiniteField(Field):
    """
    Abstract base class for finite fields.
    """
    def __init__(self, base, names, normalize):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: K = GF(7); K
            Finite Field of size 7
            sage: loads(K.dumps()) == K
            True
            sage: GF(7^10, 'a')
            Finite Field in a of size 7^10
            sage: K = GF(7^10, 'a'); K
            Finite Field in a of size 7^10
            sage: loads(K.dumps()) == K
            True
        """
        Field.__init__(self, base, names, normalize, category=_FiniteFields)

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

        The :func:`~sage.misc.latex.latex` function calls the
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
        Return string representation of ``self`` that Magma can
        understand.

        EXAMPLES::

            sage: GF(97,'a')._magma_init_(magma)            # optional - magma
            'GF(97)'
            sage: GF(9,'a')._magma_init_(magma)             # optional - magma
            'SageCreateWithNames(ext<GF(3)|_sage_[...]![GF(3)!2,GF(3)!2,GF(3)!1]>,["a"])'
            sage: magma(GF(9,'a'))                          # optional - magma
            Finite field of size 3^2
            sage: magma(GF(9,'a')).1                        # optional - magma
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
        Returns the string representation of ``self`` that Macaulay2 can
        understand.

        EXAMPLES::

            sage: GF(97,'a')._macaulay2_init_()
            'GF 97'

            sage: macaulay2(GF(97, 'a'))       # optional - macaulay2
            GF 97
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

    def _cmp_(left, Parent right):
        """
        Compares this finite field with other.

        .. WARNING::

            The notation of equality of finite fields in Sage is
            currently not stable, i.e., it may change in a future version.

        EXAMPLES::

            sage: FiniteField(3**2, 'c') == FiniteField(3**3, 'c') # indirect doctest
            False
            sage: FiniteField(3**2, 'c') == FiniteField(3**2, 'c')
            True

        The variable name is (currently) relevant for comparison of finite
        fields::

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

        EXAMPLES::

            sage: from sage.rings.finite_rings.finite_field_ext_pari import FiniteField_ext_pari
            sage: k = FiniteField_ext_pari(8, 'a')
            sage: i = iter(k); i # indirect doctest
            <sage.rings.finite_rings.finite_field_base.FiniteFieldIterator object at ...>
            sage: i.next()
            0
            sage: list(k) # indirect doctest
            [0, 1, a, a + 1, a^2, a^2 + 1, a^2 + a, a^2 + a + 1]
        """
        return FiniteFieldIterator(self)

    def _is_valid_homomorphism_(self, codomain, im_gens):
        """
        Return ``True`` if the map from self to codomain sending
        ``self.0`` to the unique element of ``im_gens`` is a valid field
        homomorphism. Otherwise, return ``False``.

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
        Return homset of homomorphisms from ``self`` to the finite field
        codomain. This function is implicitly called by the Hom method or
        function.

        The ``cat`` option is currently ignored.

        EXAMPLES::

            sage: K.<a> = GF(25); K
            Finite Field in a of size 5^2
            sage: K.Hom(K) # indirect doctest
            Automorphism group of Finite Field in a of size 5^2
        """
        return FiniteFieldHomset(self, codomain)

    def gen(self):
        r"""
        Return a generator of this field (over its prime field). As this is an
        abstract base class, this just raises a ``NotImplementedError``.

        EXAMPLES::

            sage: K = GF(17)
            sage: sage.rings.finite_rings.finite_field_base.FiniteField.gen(K)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def zeta_order(self):
        """
        Return the order of the distinguished root of unity in ``self``.

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
        Returns an element of multiplicative order ``n`` in this
        finite field, if there is one.  Raises a ``ValueError`` if there
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
        Return a primitive element of this finite field, i.e. a generator
        of the multiplicative group.

        You can use :meth:`multiplicative_generator()` or
        :meth:`primitive_element()`, these mean the same thing.

        .. WARNING::

           This generator might change from one version of Sage to another.

        EXAMPLES::

            sage: k = GF(997)
            sage: k.multiplicative_generator()
            7
            sage: k.<a> = GF(11^3)
            sage: k.primitive_element()
            a
            sage: k.<b> = GF(19^32)
            sage: k.multiplicative_generator()
            b + 4

        TESTS:

        Check that large characteristics work (:trac:`11946`)::

            sage: p = 10^20 + 39
            sage: x = polygen(GF(p))
            sage: K.<a> = GF(p^2, modulus=x^2+1)
            sage: K.multiplicative_generator()
            a + 12
        """
        from sage.rings.arith import primitive_root

        if self.__multiplicative_generator is not None:
            return self.__multiplicative_generator
        if self.degree() == 1:
            self.__multiplicative_generator = self(primitive_root(self.order()))
            return self.__multiplicative_generator
        n = self.order() - 1
        g = self.gen(0)
        # We check whether x+g is a multiplicative generator, where
        # x runs through the finite field.
        # This has the advantage that g is the first element we try,
        # so we always get g as generator if possible.  Second, the
        # PARI finite field iterator gives all the constant elements
        # first, so we try g+(constant) before anything else.
        for x in self:
            a = g+x
            if a != 0 and a.multiplicative_order() == n:
                self.__multiplicative_generator = a
                return a

    primitive_element = multiplicative_generator

    def ngens(self):
        """
        The number of generators of the finite field.  Always 1.

        EXAMPLES::

            sage: k = FiniteField(3^4, 'b')
            sage: k.ngens()
            1
        """
        return 1

    def is_field(self, proof = True):
        """
        Returns whether or not the finite field is a field, i.e.,
        always returns ``True``.

        EXAMPLES::

            sage: k.<a> = FiniteField(3^4)
            sage: k.is_field()
            True
        """
        return True

    def is_finite(self):
        """
        Return ``True`` since a finite field is finite.

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

    def factored_order(self):
        """
        Returns the factored order of this field.  For compatibility with
        :mod:`~sage.rings.finite_rings.integer_mod_ring`.

        EXAMPLES::

            sage: GF(7^2,'a').factored_order()
            7^2
        """
        from sage.structure.factorization import Factorization
        return Factorization([(self.characteristic(), self.degree())])

    def factored_unit_order(self):
        """
        Returns the factorization of ``self.order()-1``, as a 1-element list.

        The format is for compatibility with
        :mod:`~sage.rings.finite_rings.integer_mod_ring`.

        EXAMPLES::

            sage: GF(7^2,'a').factored_unit_order()
            [2^4 * 3]
        """
        if self.__factored_unit_order is None:
            if self.characteristic() in []: # want to be [2,3,5,7,11] once #7240 is finished.
                from sage.rings.factorint import factor_cunningham
                self.__factored_unit_order = [factor_cunningham(self.order()-1)]
            else:
                self.__factored_unit_order = [(self.order()-1).factor()]
        return self.__factored_unit_order

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        Same as :meth:`order`.

        EXAMPLES::

            sage: GF(997).cardinality()
            997
        """
        return self.order()

    __len__ = cardinality

    def is_prime_field(self):
        """
        Return ``True`` if ``self`` is a prime field, i.e., has degree 1.

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
        ``F[x](f)`` coerces `f` as a polynomial with coefficients in `F`.
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


    def random_element(self, *args, **kwds):
        r"""
        A random element of the finite field.  Passes arguments to
        ``random_element()`` function of underlying vector space.

        EXAMPLES::

            sage: k = GF(19^4, 'a')
            sage: k.random_element()
            a^3 + 3*a^2 + 6*a + 9

        Passes extra positional or keyword arguments through::

            sage: k.random_element(prob=0)
            0

        """
        if self.degree() == 1:
            return self(randrange(self.order()))
        v = self.vector_space().random_element(*args, **kwds)
        return self(v)

    def some_elements(self):
        """
        Returns a collection of elements of this finite field for use in unit
        testing.

        EXAMPLES::

            sage: k = GF(2^8,'a')
            sage: k.some_elements() # random output
            [a^4 + a^3 + 1, a^6 + a^4 + a^3, a^5 + a^4 + a, a^2 + a]
        """
        return [self.random_element() for i in range(4)]

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
        from sage.rings.finite_rings.constructor import FiniteField as GF

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

    def algebraic_closure(self):
        """
        Return the algebraic closure of ``self`` (not implemented).

        .. NOTE::

           This is not yet implemented for finite fields.

        EXAMPLES::

            sage: GF(5).algebraic_closure()
            Traceback (most recent call last):
            ...
            NotImplementedError: Algebraic closures of finite fields not implemented.
        """
        raise NotImplementedError, "Algebraic closures of finite fields not implemented."

    @cached_method
    def is_conway(self):
        """
        Return ``True`` if self is defined by a Conway polynomial.

        EXAMPLES:

            sage: GF(5^3, 'a').is_conway()
            True
            sage: GF(5^3, 'a', modulus='adleman-lenstra').is_conway()
            False
            sage: GF(next_prime(2^16, 2), 'a').is_conway()
            False
        """
        from conway_polynomials import conway_polynomial, exists_conway_polynomial
        p = self.characteristic()
        n = self.degree()
        return (exists_conway_polynomial(p, n)
                and self.polynomial() == self.polynomial_ring()(conway_polynomial(p, n)))


def unpickle_FiniteField_ext(_type, order, variable_name, modulus, kwargs):
    r"""
    Used to unpickle extensions of finite fields. Now superseded (hence no
    doctest), but kept around for backward compatibility.

    EXAMPLES::

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
    Return ``True`` if ``x`` is of type finite field, and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.rings.finite_rings.finite_field_base import is_FiniteField
        sage: is_FiniteField(GF(9,'a'))
        True
        sage: is_FiniteField(GF(next_prime(10^10)))
        True

    Note that the integers modulo n are not of type finite field,
    so this function returns ``False``::

        sage: is_FiniteField(Integers(7))
        False
    """
    return IS_INSTANCE(x, FiniteField)

