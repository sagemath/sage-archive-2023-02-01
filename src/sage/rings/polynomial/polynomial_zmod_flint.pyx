"""
Dense univariate polynomials over `\ZZ/n\ZZ`, implemented using FLINT.

This module gives a fast implementation of `(\ZZ/n\ZZ)[x]` whenever `n` is at
most ``sys.maxsize``. We use it by default in preference to NTL when the modulus
is small, falling back to NTL if the modulus is too large, as in the example
below.

EXAMPLES::

    sage: R.<a> = PolynomialRing(Integers(100))
    sage: type(a)
    <type 'sage.rings.polynomial.polynomial_zmod_flint.Polynomial_zmod_flint'>
    sage: R.<a> = PolynomialRing(Integers(5*2^64))
    sage: type(a)
    <type 'sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_modn_ntl_ZZ'>
    sage: R.<a> = PolynomialRing(Integers(5*2^64), implementation="FLINT")
    Traceback (most recent call last):
    ...
    ValueError: FLINT does not support modulus 92233720368547758080

AUTHORS:

- Burcin Erocal (2008-11) initial implementation
- Martin Albrecht (2009-01) another initial implementation
"""
#*****************************************************************************
#       Copyright (C) 2009-2010 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2009 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.ntl.ntl_lzz_pX import ntl_zz_pX
from sage.structure.factorization import Factorization
from sage.structure.element import coerce_binop, parent
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint

# We need to define this stuff before including the templating stuff
# to make sure the function get_cparent is found since it is used in
# 'polynomial_template.pxi'.

cdef inline cparent get_cparent(parent) except? 0:
    try:
        return <unsigned long>(parent.modulus())
    except AttributeError:
        return 0

# first we include the definitions
include "sage/libs/flint/nmod_poly_linkage.pxi"

# and then the interface
include "polynomial_template.pxi"

cdef extern from "zn_poly/zn_poly.h":
    ctypedef struct zn_mod_struct:
        pass
    cdef void zn_mod_init(zn_mod_struct *mod, unsigned long m)
    cdef void zn_mod_clear(zn_mod_struct *mod)
    cdef void zn_array_mul(unsigned long* res, unsigned long* op1, size_t n1, unsigned long* op2, size_t n2, zn_mod_struct *mod)

from sage.libs.flint.fmpz_poly cimport *

cdef class Polynomial_zmod_flint(Polynomial_template):
    r"""
    Polynomial on `\ZZ/n\ZZ` implemented via FLINT.

    TESTS::

        sage: f = Integers(4)['x'].random_element()
        sage: from sage.rings.polynomial.polynomial_zmod_flint import Polynomial_zmod_flint
        sage: isinstance(f, Polynomial_zmod_flint)
        True

    .. automethod:: _add_
    .. automethod:: _sub_
    .. automethod:: _lmul_
    .. automethod:: _rmul_
    .. automethod:: _mul_
    .. automethod:: _mul_trunc_
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        EXAMPLE::

            sage: P.<x> = GF(32003)[]
            sage: f = 24998*x^2 + 29761*x + 2252
        """
        cdef long nlen

        if isinstance(x, list) or isinstance(x, tuple):
            k = parent._base
            if check:
                lst = [k(i) for i in x]
            else:
                lst = x
            # remove trailing zeroes
            nlen = len(lst)
            while nlen and lst[nlen-1] == 0:
                nlen -= 1
            lst = lst[:nlen]
            Polynomial_template.__init__(self, parent, 0, check, is_gen, construct)
            self._set_list(lst)
            return
        elif isinstance(x, Polynomial_integer_dense_flint):
            Polynomial_template.__init__(self, parent, 0, check, is_gen, construct)
            self._set_fmpz_poly((<Polynomial_integer_dense_flint>x).__poly)
            return
        else:
            if isinstance(x, ntl_zz_pX):
                x = x.list()
            try:
                if x.parent() is parent.base_ring() or x.parent() == parent.base_ring():
                    x = int(x) % parent.modulus()
            except AttributeError:
                pass
        Polynomial_template.__init__(self, parent, x, check, is_gen, construct)

    cdef Polynomial_template _new(self):
        """
        EXAMPLES::

            sage: P.<x> = GF(5)[]
            sage: (2*x+1).monic() #indirect doctest
            x + 3
        """
        cdef type t = type(self)
        cdef Polynomial_template e = <Polynomial_template>t.__new__(t)
        nmod_poly_init(&e.x, self._parent.modulus())
        e._parent = self._parent
        e._cparent = self._cparent
        return e

    cpdef Polynomial _new_constant_poly(self, x, Parent P):
        r"""
        Quickly creates a new constant polynomial with value x in parent P.

        ASSUMPTION:

        x must convertible to an int.

        The modulus of P must coincide with the modulus of this element.
        That assumption is not verified!

        EXAMPLE::

            sage: R.<x> = GF(3)[]
            sage: x._new_constant_poly(4,R)
            1
            sage: x._new_constant_poly('4',R)
            1
            sage: x._new_constant_poly('4.1',R)
            Traceback (most recent call last):
            ...
            ValueError: invalid literal for int() with base 10: '4.1'

        """
        cdef type t = type(self)
        cdef Polynomial_template r = <Polynomial_template>t.__new__(t)
        r._parent = P
        r._cparent = get_cparent(P)
        nmod_poly_init(&r.x, nmod_poly_modulus(&self.x))
        celement_set_si(&r.x, int(x), (<Polynomial_template>self)._cparent)
        return r

    cdef int _set_list(self, x) except -1:
        """
        Set the coefficients of ``self`` from a list of coefficients.

        INPUT:

        - ``x`` - a list of coefficients - the coefficients are assumed to be
          reduced already and the list contains no trailing zeroes.


        EXAMPLES::

            sage: P.<a>=GF(7)[]
            sage: P([2^60,0,1])
            a^2 + 1
            sage: P([])
            0
            sage: P(range(15))
            6*a^13 + 5*a^12 + 4*a^11 + 3*a^10 + 2*a^9 + a^8 + 6*a^6 + 5*a^5 + 4*a^4 + 3*a^3 + 2*a^2 + a
        """
        cdef list l_in = x
        cdef unsigned long length = len(l_in)
        cdef unsigned long modulus = nmod_poly_modulus(&self.x)
        cdef int i
        if length == 0:
            nmod_poly_zero(&self.x)
            return 0

        # resize to length of list
        sig_on()
        nmod_poly_realloc(&self.x, length)
        sig_off()

        sig_on()
        # The following depends on the internals of FLINT
        for i from 0 <= i < length:
            self.x.coeffs[i] = l_in[i]
        self.x.length = length
        sig_off()
        return 0

    cdef int _set_fmpz_poly(self, fmpz_poly_t x) except -1:
        """
        Set the coefficients of ``self`` from the coefficients of an ``fmpz_poly_t`` element.

        INPUT:

        - ``x`` - an ``fmpz_poly_t`` element

        EXAMPLES::

            sage: a = ZZ['x'](range(17))
            sage: R = Integers(7)['x']
            sage: R(a)
            2*x^16 + x^15 + 6*x^13 + 5*x^12 + 4*x^11 + 3*x^10 + 2*x^9 + x^8 + 6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x

        TESTS:

        The following test from :trac:`12173` used to be horribly slow::

            sage: a = ZZ['x'](range(100000))
            sage: R = Integers(3)['x']
            sage: R(a)  # long time (7s on sage.math, 2013)
            2*x^99998 + ... + x
        """
        sig_on()
        fmpz_poly_get_nmod_poly(&self.x, x)
        sig_off()
        return 0

    cdef get_unsafe(self, Py_ssize_t i):
        """
        Return the `i`-th coefficient of ``self``.

        EXAMPLES::

            sage: P.<x> = GF(32003)[]
            sage: f = 24998*x^2 + 29761*x + 2252
            sage: f[100]
            0
            sage: f[1]
            29761
            sage: f[0]
            2252
            sage: f[-1]
            0
            sage: f[:2]
            29761*x + 2252
            sage: f[:50] == f
            True
        """
        cdef unsigned long c = nmod_poly_get_coeff_ui(&self.x, i)
        return self._parent.base_ring()(c)

    def __call__(self, *x, **kwds):
        """
        Evaluate polynomial at x=a.

        INPUT: **either**

        - a -- ring element; need not be in the coefficient ring of the
          polynomial.
        - a dictionary for kwds:value pairs.  If the variable name of the
          polynomial is a keyword it is substituted in; otherwise this
          polynomial is returned unchanged.

        EXAMPLE::

            sage: P.<x> = PolynomialRing(GF(7))
            sage: f = x^2 + 1
            sage: f(0)
            1
            sage: f(2)
            5
            sage: f(3)
            3

            sage: f(x+1)
            x^2 + 2*x + 2

        Test some simple (but important) optimizations::

            sage: f(2) == f(P(2))
            True
            sage: f(x) is f
            True
            sage: f(1/x)
            (x^2 + 1)/x^2
        """
        cdef Polynomial_zmod_flint t, y
        cdef long c
        K = self._parent.base_ring()
        if not kwds and len(x) == 1:
            P = parent(x[0])
            if K.has_coerce_map_from(P):
                x = K(x[0])
                return K(nmod_poly_evaluate_nmod(&self.x, x))
            elif self._parent.has_coerce_map_from(P):
                y = <Polynomial_zmod_flint>self._parent(x[0])
                t = self._new()
                if nmod_poly_degree(&y.x) == 0:
                    c = nmod_poly_evaluate_nmod(&self.x, nmod_poly_get_coeff_ui(&y.x, 0))
                    nmod_poly_set_coeff_ui(&t.x, 0, c)
                elif nmod_poly_degree(&y.x) == 1 and nmod_poly_get_coeff_ui(&y.x, 0) == 0:
                    c = nmod_poly_get_coeff_ui(&y.x, 1)
                    if c == 1:
                        return self
                nmod_poly_compose(&t.x, &self.x, &y.x)
                return t
        return Polynomial.__call__(self, *x, **kwds)

    @coerce_binop
    def resultant(self, Polynomial_zmod_flint other):
        """
        Returns the resultant of self and other, which must lie in the same
        polynomial ring.

        INPUT:

        - other -- a polynomial

        OUTPUT: an element of the base ring of the polynomial ring

        EXAMPLES::

            sage: R.<x> = GF(19)['x']
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: r = f.resultant(g); r
            11
            sage: r.parent() is GF(19)
            True

        The following example shows that :trac:`11782` has been fixed::

            sage: R.<x> = ZZ.quo(9)['x']
            sage: f = 2*x^3 + x^2 + x;  g = 6*x^2 + 2*x + 1
            sage: f.resultant(g)
            5
        """
        # As of version 1.6 of FLINT, the base ring must be a field to compute
        # resultants correctly. (see http://www.flintlib.org/flint-1.6.pdf p.58)
        # If it is not a field we fall back to direct computation through the
        # Sylvester matrix.
        if self.base_ring().is_field():
            res = nmod_poly_resultant(&(<Polynomial_template>self).x,
                                      &(<Polynomial_template>other).x)
            return self.parent().base_ring()(res)
        else:
            return self.sylvester_matrix(other).determinant()

    def small_roots(self, *args, **kwds):
        r"""
        See :func:`sage.rings.polynomial.polynomial_modn_dense_ntl.small_roots`
        for the documentation of this function.

        EXAMPLE::

            sage: N = 10001
            sage: K = Zmod(10001)
            sage: P.<x> = PolynomialRing(K)
            sage: f = x^3 + 10*x^2 + 5000*x - 222
            sage: f.small_roots()
            [4]
        """
        from sage.rings.polynomial.polynomial_modn_dense_ntl import small_roots
        return small_roots(self, *args, **kwds)

    def _unsafe_mutate(self, n, value):
        r"""
        Never use this unless you really know what you are doing.

        INPUT:

        - n -- degree
        - value -- coefficient

        .. warning::

            This could easily introduce subtle bugs, since Sage assumes
            everywhere that polynomials are immutable.  It's OK to use this if
            you really know what you're doing.

        EXAMPLES::

            sage: R.<x> = GF(7)[]
            sage: f = (1+2*x)^2; f
            4*x^2 + 4*x + 1
            sage: f._unsafe_mutate(1, -5)
            sage: f
            4*x^2 + 2*x + 1
        """
        n = int(n)
        value = self.base_ring()(value)
        if n >= 0:
            nmod_poly_set_coeff_ui(&self.x, n, int(value)%nmod_poly_modulus(&self.x))
        else:
            raise IndexError("Polynomial coefficient index must be nonnegative.")

    def _mul_zn_poly(self, other):
        r"""
        Returns the product of two polynomials using the zn_poly library.

        See http://www.math.harvard.edu/~dmharvey/zn_poly/ for details
        on zn_poly.

        INPUT:

        - self: Polynomial
        - right: Polynomial (over same base ring as self)

        OUTPUT: (Polynomial) the product self*right.


        EXAMPLE::

            sage: P.<x> = PolynomialRing(GF(next_prime(2^30)))
            sage: f = P.random_element(1000)
            sage: g = P.random_element(1000)
            sage: f*g == f._mul_zn_poly(g)
            True

            sage: P.<x> = PolynomialRing(Integers(100))
            sage: P
            Univariate Polynomial Ring in x over Ring of integers modulo 100
            sage: r = (10*x)._mul_zn_poly(10*x); r
            0
            sage: r.degree()
            -1

        ALGORITHM:

           uses David Harvey's zn_poly library.

        NOTE: This function is a technology preview. It might
        disappear or be replaced without a deprecation warning.
        """
        cdef Polynomial_zmod_flint _other = <Polynomial_zmod_flint>self._parent._coerce_(other)

        cdef type t = type(self)
        cdef Polynomial_zmod_flint r = <Polynomial_zmod_flint>t.__new__(t)
        r._parent = (<Polynomial_zmod_flint>self)._parent
        r._cparent = (<Polynomial_zmod_flint>self)._cparent

        cdef unsigned long p = self._parent.modulus()
        cdef unsigned long n1 = self.x.length
        cdef unsigned long n2 = _other.x.length

        cdef zn_mod_struct zn_mod

        nmod_poly_init2(&r.x, p, n1 + n2 -1 )

        zn_mod_init(&zn_mod, p)
        zn_array_mul(<unsigned long *> r.x.coeffs, <unsigned long *> self.x.coeffs, n1, <unsigned long*> _other.x.coeffs, n2, &zn_mod)
        r.x.length = n1 + n2 -1
        _nmod_poly_normalise(&r.x)
        zn_mod_clear(&zn_mod)
        return r

    cpdef Polynomial _mul_trunc_(self, Polynomial right, long n):
        """
        Return the product of this polynomial and other truncated to the
        given length `n`.

        This function is usually more efficient than simply doing the
        multiplication and then truncating. The function is tuned for length
        `n` about half the length of a full product.


        EXAMPLES::

            sage: P.<a>=GF(7)[]
            sage: a = P(range(10)); b = P(range(5, 15))
            sage: a._mul_trunc_(b, 5)
            4*a^4 + 6*a^3 + 2*a^2 + 5*a

        TESTS::

            sage: a._mul_trunc_(b, 0)
            Traceback (most recent call last):
            ...
            ValueError: length must be > 0
        """
        if n <= 0:
            raise ValueError("length must be > 0")
        cdef Polynomial_zmod_flint op2 = <Polynomial_zmod_flint> right
        cdef Polynomial_zmod_flint res = self._new()
        nmod_poly_mullow(&res.x, &self.x, &op2.x, n)
        return res

    _mul_short = _mul_trunc_

    cpdef _mul_trunc_opposite(self, Polynomial_zmod_flint other, n):
        """
        Return the product of this polynomial and other ignoring the least
        significant `n` terms of the result which may be set to anything.

        This function is more efficient than doing the full multiplication if
        the operands are relatively short. It is tuned for `n` about half the
        length of a full product.

        EXAMPLES::

            sage: P.<a>=GF(7)[]
            sage: b = P(range(10)); c = P(range(5, 15))
            sage: b._mul_trunc_opposite(c, 10)
            5*a^17 + 2*a^16 + 6*a^15 + 4*a^14 + 4*a^13 + 5*a^10 + 2*a^9 + 5*a^8 + 4*a^5 + 4*a^4 + 6*a^3 + 2*a^2 + 5*a
            sage: list(b._mul_trunc_opposite(c, 10))[10:18]
            [5, 0, 0, 4, 4, 6, 2, 5]
            sage: list(b*c)[10:18]
            [5, 0, 0, 4, 4, 6, 2, 5]
            sage: list(b._mul_trunc_opposite(c, 18))[18:]
            []

        TESTS::

            sage: a._mul_trunc_opposite(b, -1)
            Traceback (most recent call last):
            ...
            ValueError: length must be >= 0
        """
        cdef Polynomial_zmod_flint res = self._new()
        if n < 0:
            raise ValueError("length must be >= 0")
        nmod_poly_mulhigh(&res.x, &self.x, &other.x, n)
        return res

    _mul_short_opposite = _mul_trunc_opposite

    cpdef rational_reconstruct(self, m, n_deg=0, d_deg=0):
        """
        Construct a rational function n/d such that `p*d` is equivalent to `n`
        modulo `m` where `p` is this polynomial.

        EXAMPLES::

            sage: P.<x> = GF(5)[]
            sage: p = 4*x^5 + 3*x^4 + 2*x^3 + 2*x^2 + 4*x + 2
            sage: n, d = p.rational_reconstruct(x^9, 4, 4); n, d
            (3*x^4 + 2*x^3 + x^2 + 2*x, x^4 + 3*x^3 + x^2 + x)
            sage: (p*d % x^9) == n
            True
        """
        if n_deg < 0 or d_deg < 0:
            raise ValueError("The degree bounds n_deg and d_deg should be positive.")

        if n_deg == 0:
            n_deg = (m.degree() - 1)//2
        if d_deg == 0:
            d_deg = (m.degree() - 1)//2
        P = self._parent

        cdef Polynomial_zmod_flint s0 = self._new()
        cdef Polynomial_zmod_flint t0 = P.one()
        cdef Polynomial_zmod_flint s1 = m
        cdef Polynomial_zmod_flint t1 = self%m

        cdef Polynomial_zmod_flint q
        cdef Polynomial_zmod_flint r0
        cdef Polynomial_zmod_flint r1

        while nmod_poly_length(&t1.x) != 0 and n_deg < nmod_poly_degree(&t1.x):
            q = self._new()
            r1 = self._new()
            nmod_poly_divrem(&q.x, &r1.x, &s1.x, &t1.x)
            r0 = s0 - q*t0
            s0 = t0
            s1 = t1
            t0 = r0
            t1 = r1

        assert(t0 != 0)
        if d_deg < nmod_poly_degree(&t0.x):
            raise ValueError("could not complete rational reconstruction")

        # make the denominator monic
        c = t0.leading_coefficient()
        t0 = t0.monic()
        t1 = t1/c

        return t1, t0

    def is_irreducible(self):
        """
        Return True if this polynomial is irreducible.

        EXAMPLES::

            sage: R.<x> = GF(5)[]
            sage: (x^2 + 1).is_irreducible()
            False
            sage: (x^3 + x + 1).is_irreducible()
            True

        TESTS::

            sage: R(0).is_irreducible()
            False
            sage: R(1).is_irreducible()
            False
            sage: R(2).is_irreducible()
            False

            sage: S.<s> = Zmod(10)[]
            sage: (s^2).is_irreducible()
            Traceback (most recent call last):
            ...
            NotImplementedError: checking irreducibility of polynomials over rings with composite characteristic is not implemented
            sage: S(1).is_irreducible()
            False
            sage: S(2).is_irreducible()
            Traceback (most recent call last):
            ...
            NotImplementedError: checking irreducibility of polynomials over rings with composite characteristic is not implemented
        """
        if not self:
            return False
        if self.is_unit():
            return False

        if not self.base_ring().is_field():
            raise NotImplementedError("checking irreducibility of polynomials over rings with composite characteristic is not implemented")

        sig_on()
        if 1 == nmod_poly_is_irreducible(&self.x):
            sig_off()
            return True
        else:
            sig_off()
            return False

    def squarefree_decomposition(self):
        """
        Returns the squarefree decomposition of this polynomial.

        EXAMPLES::

            sage: R.<x> = GF(5)[]
            sage: ((x+1)*(x^2+1)^2*x^3).squarefree_decomposition()
            (x + 1) * (x^2 + 1)^2 * x^3

        TESTS::

            sage: (2*x*(x+1)^2).squarefree_decomposition()
            (2) * x * (x + 1)^2
            sage: P.<x> = Zmod(10)[]
            sage: (x^2).squarefree_decomposition()
            Traceback (most recent call last):
            ...
            NotImplementedError: square free factorization of polynomials over rings with composite characteristic is not implemented

        """
        if not self.base_ring().is_field():
            raise NotImplementedError("square free factorization of polynomials over rings with composite characteristic is not implemented")

        return factor_helper(self, True)

    def factor(self):
        """
        Returns the factorization of the polynomial.

        EXAMPLES::

            sage: R.<x> = GF(5)[]
            sage: (x^2 + 1).factor()
            (x + 2) * (x + 3)

        TESTS::

            sage: (2*x^2 + 2).factor()
            (2) * (x + 2) * (x + 3)
            sage: P.<x> = Zmod(10)[]
            sage: (x^2).factor()
            Traceback (most recent call last):
            ...
            NotImplementedError: factorization of polynomials over rings with composite characteristic is not implemented

        """
        if not self.base_ring().is_field():
            raise NotImplementedError("factorization of polynomials over rings with composite characteristic is not implemented")

        return factor_helper(self)

    def monic(self):
        """
        Return this polynomial divided by its leading coefficient.

        Raises ValueError if the leading cofficient is not invertible in the
        base ring.

        EXAMPLES::

            sage: R.<x> = GF(5)[]
            sage: (2*x^2+1).monic()
            x^2 + 3

        TESTS::

            sage: R.<x> = Zmod(10)[]
            sage: (5*x).monic()
            Traceback (most recent call last):
            ...
            ValueError: leading coefficient must be invertible
        """
        if self.base_ring().characteristic().gcd(\
                self.leading_coefficient().lift()) != 1:
            raise ValueError("leading coefficient must be invertible")
        cdef Polynomial_zmod_flint res = self._new()
        nmod_poly_make_monic(&res.x, &self.x)
        return res

    def reverse(self, degree=None):
        """
        Return a polynomial with the coefficients of this polynomial reversed.

        If an optional degree argument is given the coefficient list will be
        truncated or zero padded as necessary and the reverse polynomial will
        have the specified degree.

        EXAMPLES::

            sage: R.<x> = GF(5)[]
            sage: p = R([1,2,3,4]); p
            4*x^3 + 3*x^2 + 2*x + 1
            sage: p.reverse()
            x^3 + 2*x^2 + 3*x + 4
            sage: p.reverse(degree=6)
            x^6 + 2*x^5 + 3*x^4 + 4*x^3
            sage: p.reverse(degree=2)
            x^2 + 2*x + 3

            sage: R.<x> = GF(101)[]
            sage: f = x^3 - x + 2; f
            x^3 + 100*x + 2
            sage: f.reverse()
            2*x^3 + 100*x^2 + 1
            sage: f.reverse() == f(1/x) * x^f.degree()
            True

        Note that if `f` has zero constant coefficient, its reverse will
        have lower degree.

        ::

            sage: f = x^3 + 2*x
            sage: f.reverse()
            2*x^2 + 1

        In this case, reverse is not an involution unless we explicitly
        specify a degree.

        ::

            sage: f
            x^3 + 2*x
            sage: f.reverse().reverse()
            x^2 + 2
            sage: f.reverse(5).reverse(5)
            x^3 + 2*x

        TESTS::

            sage: p.reverse(degree=1.5r)
            Traceback (most recent call last):
            ...
            ValueError: degree argument must be a non-negative integer, got 1.5
        """
        cdef Polynomial_zmod_flint res = self._new()
        cdef unsigned long d
        if degree:
            d = degree
            if d != degree:
                raise ValueError("degree argument must be a non-negative integer, got %s"%(degree))
            nmod_poly_reverse(&res.x, &self.x, d+1) # FLINT expects length
        else:
            nmod_poly_reverse(&res.x, &self.x, nmod_poly_length(&self.x))
        return res

    def revert_series(self, n):
        r"""
        Return a polynomial `f` such that `f(self(x)) = self(f(x)) = x mod x^n`.

        EXAMPLES::

            sage: R.<t> =  GF(5)[]
            sage: f = t + 2*t^2 - t^3 - 3*t^4
            sage: f.revert_series(5)
            3*t^4 + 4*t^3 + 3*t^2 + t

            sage: f.revert_series(-1)
            Traceback (most recent call last):
            ...
            ValueError: argument n must be a non-negative integer, got -1

            sage: g = - t^3 + t^5
            sage: g.revert_series(6)
            Traceback (most recent call last):
            ...
            ValueError: self must have constant coefficient 0 and a unit for coefficient t^1

            sage: g = t + 2*t^2 - t^3 -3*t^4 + t^5
            sage: g.revert_series(6)
            Traceback (most recent call last):
            ...
            ValueError: the integers 1 up to n=5 are required to be invertible over the base field
        """
        cdef Polynomial_zmod_flint res = self._new()
        cdef unsigned long m
        if n < 0:
            raise ValueError("argument n must be a non-negative integer, got {}".format(n))
        m = n
        if not self[0].is_zero() or not self[1].is_unit():
            raise ValueError("self must have constant coefficient 0 and a unit for coefficient {}^1".format(self.parent().gen()))
        if not all((self.base_ring())(i) != 0 for i in range(1,n)):
            raise ValueError("the integers 1 up to n={} are required to be invertible over the base field".format(n-1))

        sig_on()
        nmod_poly_revert_series(&res.x, &self.x, m)
        sig_off()

        return res
