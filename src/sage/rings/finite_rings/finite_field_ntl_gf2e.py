"""
Finite Fields of Characteristic 2
"""

from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.libs.pari.all import pari
from finite_field_ext_pari import FiniteField_ext_pari
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic

def late_import():
    """
    Imports various modules after startup.

    EXAMPLES::

       sage: sage.rings.finite_rings.finite_field_ntl_gf2e.late_import()
       sage: sage.rings.finite_rings.finite_field_ntl_gf2e.GF2 is None # indirect doctest
       False
    """
    if "GF2" in globals():
        return
    global ResidueField_generic, is_FiniteField, exists_conway_polynomial, conway_polynomial, Cache_ntl_gf2e, GF, GF2, is_Polynomial
    import sage.rings.residue_field
    ResidueField_generic = sage.rings.residue_field.ResidueField_generic

    import sage.rings.finite_rings.finite_field_base
    is_FiniteField = sage.rings.finite_rings.finite_field_base.is_FiniteField

    import sage.rings.finite_rings.conway_polynomials
    exists_conway_polynomial = sage.rings.finite_rings.conway_polynomials.exists_conway_polynomial
    conway_polynomial = sage.rings.finite_rings.conway_polynomials.conway_polynomial

    import sage.rings.finite_rings.element_ntl_gf2e
    Cache_ntl_gf2e = sage.rings.finite_rings.element_ntl_gf2e.Cache_ntl_gf2e

    import sage.rings.finite_rings.constructor
    GF = sage.rings.finite_rings.constructor.GF
    GF2 = GF(2)

    import sage.rings.polynomial.polynomial_element
    is_Polynomial = sage.rings.polynomial.polynomial_element.is_Polynomial

class FiniteField_ntl_gf2e(FiniteField):
    """
    Finite Field of characteristic 2 and order `2^n`.

    INPUT:

    - ``q`` -- `2^n` (must be 2 power)

    - ``names`` -- variable used for poly_repr (default: ``'a'``)

    - ``modulus`` -- you may provide a polynomial to use for reduction or
      a string:

      - ``'conway'`` -- force the use of a Conway polynomial, will
        raise a ``RuntimeError`` if ``None`` is found in the database;
      - ``'minimal_weight'`` -- use a minimal weight polynomial, should
        result in faster arithmetic;
      - ``'random'`` -- use a random irreducible polynomial.
      - ``'default'`` -- a Conway polynomial is used if found. Otherwise
        a sparse polynomial is used.

    - ``repr`` -- controls the way elements are printed to the user:
                 (default: ``'poly'``)

      - ``'poly'``: polynomial representation

    OUTPUT:

    Finite field with characteristic 2 and cardinality `2^n`.

    EXAMPLES::

        sage: k.<a> = GF(2^16)
        sage: type(k)
        <class 'sage.rings.finite_rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e_with_category'>
        sage: k.<a> = GF(2^1024)
        sage: k.modulus()
        x^1024 + x^19 + x^6 + x + 1
        sage: set_random_seed(0)
        sage: k.<a> = GF(2^17, modulus='random')
        sage: k.modulus()
        x^17 + x^16 + x^15 + x^10 + x^8 + x^6 + x^4 + x^3 + x^2 + x + 1
        sage: k.modulus().is_irreducible()
        True
        sage: k.<a> = GF(2^211, modulus='minimal_weight')
        sage: k.modulus()
        x^211 + x^11 + x^10 + x^8 + 1
        sage: k.<a> = GF(2^211, modulus='conway')
        sage: k.modulus()
        x^211 + x^9 + x^6 + x^5 + x^3 + x + 1
        sage: k.<a> = GF(2^23, modulus='conway')
        sage: a.multiplicative_order() == k.order() - 1
        True
    """

    def __init__(self, q, names="a",  modulus=None, repr="poly"):
        """
        Initialize ``self``.

        TESTS::

            sage: k.<a> = GF(2^100, modulus='strangeinput')
            Traceback (most recent call last):
            ...
            ValueError: no such algorithm for finding an irreducible polynomial: strangeinput
            sage: k.<a> = GF(2^20) ; type(k)
            <class 'sage.rings.finite_rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e_with_category'>
            sage: loads(dumps(k)) is k
            True
            sage: k1.<a> = GF(2^16)
            sage: k2.<a> = GF(2^17)
            sage: k1 == k2
            False
            sage: k3 = k1._finite_field_ext_pari_()
            sage: k1 == k3
            False

            sage: TestSuite(k).run()

            sage: k.<a> = GF(2^64)
            sage: k._repr_option('element_is_atomic')
            False
            sage: P.<x> = PolynomialRing(k)
            sage: (a+1)*x # indirect doctest
            (a + 1)*x
        """
        late_import()
        q = Integer(q)
        if q < 2:
            raise ValueError("q must be a 2-power")
        k = q.exact_log(2)
        if q != 1 << k:
            raise ValueError("q must be a 2-power")
        p = Integer(2)
        FiniteField.__init__(self, GF(p), names, normalize=True)

        self._kwargs = {'repr':repr}

        if modulus is None or modulus == 'default':
            if exists_conway_polynomial(p, k):
                modulus = "conway"
            else:
                modulus = "minimal_weight"
        if modulus == "conway":
            modulus = conway_polynomial(p, k)
        if is_Polynomial(modulus):
            modulus = modulus.list()
        self._cache = Cache_ntl_gf2e(self, k, modulus)
        self._polynomial = {}

    def characteristic(self):
        """
        Return the characteristic of ``self`` which is 2.

        EXAMPLES::

            sage: k.<a> = GF(2^16,modulus='random')
            sage: k.characteristic()
            2
        """
        return Integer(2)

    def order(self):
        """
        Return the cardinality of this field.

        EXAMPLES::

            sage: k.<a> = GF(2^64)
            sage: k.order()
            18446744073709551616
        """
        return self._cache.order()

    def degree(self):
        r"""
        If this field has cardinality `2^n` this method returns `n`.

        EXAMPLES::

            sage: k.<a> = GF(2^64)
            sage: k.degree()
            64
        """
        return self._cache.degree()

    def _element_constructor_(self, e):
        """
        Coerces several data types to ``self``.

        INPUT:

        - ``e`` -- data to coerce

        EXAMPLES::

            sage: k.<a> = GF(2^20)
            sage: k(1) # indirect doctest
            1
            sage: k(int(2))
            0

            sage: k('a+1')
            a + 1
            sage: k('b+1')
            Traceback (most recent call last):
            ...
            NameError: name 'b' is not defined

            sage: R.<x>=GF(2)[]
            sage: k(1+x+x^10+x^55)
            a^19 + a^17 + a^16 + a^15 + a^12 + a^11 + a^8 + a^6 + a^4 + a^2 + 1

            sage: V = k.vector_space()
            sage: v = V.random_element(); v
            (1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1)
            sage: k(v)
            a^19 + a^15 + a^14 + a^13 + a^11 + a^10 + a^9 + a^6 + a^5 + a^4 + 1
            sage: vector(k(v)) == v
            True

            sage: k(pari('Mod(1,2)*a^20'))
            a^10 + a^9 + a^7 + a^6 + a^5 + a^4 + a + 1
        """
        return self._cache.import_data(e)

    def gen(self, ignored=None):
        r"""
        Return a generator of ``self``.

        EXAMPLES::

            sage: k.<a> = GF(2^19)
            sage: k.gen() == a
            True
            sage: a
            a
        """
        return self._cache._gen

    def prime_subfield(self):
        r"""
        Return the prime subfield `\GF{p}` of ``self`` if ``self`` is
        `\GF{p^n}`.

        EXAMPLES::

            sage: F.<a> = GF(2^16)
            sage: F.prime_subfield()
            Finite Field of size 2
        """
        return GF2

    def fetch_int(self, number):
        r"""
        Given an integer `n` less than :meth:`cardinality` with base `2`
        representation `a_0 + 2 \cdot a_1 + \cdots + 2^k a_k`, returns
        `a_0 + a_1 \cdot x + \cdots + a_k x^k`, where `x` is the
        generator of this finite field.

        INPUT:

        - ``number`` -- an integer

        EXAMPLES::

            sage: k.<a> = GF(2^48)
            sage: k.fetch_int(2^43 + 2^15 + 1)
            a^43 + a^15 + 1
            sage: k.fetch_int(33793)
            a^15 + a^10 + 1
            sage: 33793.digits(2) # little endian
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        """
        return self._cache.fetch_int(number)

    def polynomial(self, name=None):
        """
        Return the defining polynomial of this field as an element of
        :class:`PolynomialRing`.

        This is the same as the characteristic polynomial of the
        generator of ``self``.

        INPUT:

        ``name`` -- optional variable name

        EXAMPLES::

            sage: k.<a> = GF(2^20)
            sage: k.polynomial()
            a^20 + a^10 + a^9 + a^7 + a^6 + a^5 + a^4 + a + 1
            sage: k.polynomial('FOO')
            FOO^20 + FOO^10 + FOO^9 + FOO^7 + FOO^6 + FOO^5 + FOO^4 + FOO + 1
            sage: a^20
            a^10 + a^9 + a^7 + a^6 + a^5 + a^4 + a + 1

        """
        try:
            return self._polynomial[name]
        except KeyError:
            R = self.polynomial_ring(name)
            f = R(self._cache.polynomial())
            self._polynomial[name] = f
            return f

    def _finite_field_ext_pari_(self):
        """
        Return a :class:`FiniteField_ext_pari` isomorphic to ``self`` with
        the same defining polynomial.

        .. NOTE::

            This method will vanish eventually because that implementation of
            finite fields will be deprecated.

        EXAMPLES::

            sage: k.<a> = GF(2^20)
            sage: kP = k._finite_field_ext_pari_()
            sage: kP
            Finite Field in a of size 2^20
            sage: type(kP)
            <class 'sage.rings.finite_rings.finite_field_ext_pari.FiniteField_ext_pari_with_category'>
        """
        f = self.polynomial()
        return FiniteField_ext_pari(self.order(), self.variable_name(), f)

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: k1.<a> = GF(2^16)
            sage: {k1:1} # indirect doctest
            {Finite Field in a of size 2^16: 1}
        """
        try:
            return self._hash
        except AttributeError:
            self._hash = hash((self.characteristic(),self.polynomial(),self.variable_name(),"ntl_gf2e"))
            return self._hash

    def _pari_modulus(self):
        """
        Return PARI object which is equivalent to the
        polynomial/modulus of ``self``.

        EXAMPLES::

            sage: k1.<a> = GF(2^16)
            sage: k1._pari_modulus()
            Mod(1, 2)*a^16 + Mod(1, 2)*a^5 + Mod(1, 2)*a^3 + Mod(1, 2)*a^2 + Mod(1, 2)
        """
        f = pari(str(self.modulus()))
        return f.subst('x', 'a') * pari("Mod(1,%s)"%self.characteristic())
