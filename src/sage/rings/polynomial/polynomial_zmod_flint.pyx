"""
Univariate Polynomials over Z/nZ for n <= sys.maxint via FLINT.

AUTHOR:
    -- Martin Albrecht (2009-01) another initial implementation
    -- Burcin Erocal (2008-11) initial implementation
"""

from sage.libs.ntl.ntl_lzz_pX import ntl_zz_pX

# We need to define this stuff before including the templating stuff
# to make sure the function get_cparent is found since it is used in
# 'polynomial_template.pxi'.

ctypedef unsigned long cparent

cdef inline cparent get_cparent(parent):
    if parent is None:
        return 0
    try:
        m = parent.modulus()
    except:
        return 0
    return <unsigned long>(parent.modulus())

# first we include the definitions
include "../../libs/flint/zmod_poly_linkage.pxi"

# and then the interface
include "polynomial_template.pxi"

from sage.libs.all import pari

cdef extern from "zn_poly/zn_poly.h":
    ctypedef struct zn_mod_struct:
        pass
    cdef void zn_mod_init(zn_mod_struct *mod, unsigned long m)
    cdef void zn_mod_clear(zn_mod_struct *mod)
    cdef void zn_array_mul(unsigned long* res, unsigned long* op1, size_t n1, unsigned long* op2, size_t n2, zn_mod_struct *mod)

cdef class Polynomial_zmod_flint(Polynomial_template):
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        EXAMPLE:
            sage: P.<x> = GF(32003)[]
            sage: f = 24998*x^2 + 29761*x + 2252
        """
        cdef long nlen

        if PY_TYPE_CHECK(x, list) or PY_TYPE_CHECK(x, tuple):
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
        else:
            if PY_TYPE_CHECK(x, ntl_zz_pX):
                x = x.list()
            try:
                if x.parent() is parent.base_ring() or x.parent() == parent.base_ring():
                    x = int(x) % parent.modulus()
            except AttributeError:
                pass
        Polynomial_template.__init__(self, parent, x, check, is_gen, construct)

    cdef _set_list(self, x):
        """
        EXAMPLES:
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
        cdef unsigned long modulus = zmod_poly_modulus(&self.x)
        cdef int i
        if length == 0:
            zmod_poly_zero(&self.x)
            return

        # resize to length of list
        _sig_on
        zmod_poly_realloc(&self.x, length)
        _sig_off

        for i from 0 <= i < length:
            _zmod_poly_set_coeff_ui(&self.x, i, l_in[i])
        # we need to set the length manually, we used _zmod_poly_set_coeff_ui
        self.x.length = length

    def __getitem__(self, i):
        """
        EXAMPLE:
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
        """
        cdef unsigned long c = 0
        if 0 <= i < zmod_poly_length(&self.x):
            c = zmod_poly_get_coeff_ui(&self.x, i)
        return self._parent.base_ring()(c)

    def __call__(self, *x, **kwds):
        """
        Evaluate polynomial at x=a.

        INPUT:
            a -- ring element a; need not be in the coefficient
                 ring of the polynomial.
          -- or --
            a dictionary for kwds:value pairs.  If the variable
            name of the polynomial is a kwds it is substituted in;
            otherwise this polynomial is returned unchanged.

        EXAMPLE:
            sage: P.<x> = PolynomialRing(GF(7))
            sage: f= x^2 + 1
            sage: f(0)
            1
            sage: f(2)
            5
            sage: f(3)
            3
        """
        K = self._parent.base_ring()
        if len(kwds) == 0 and len(x) == 1:
            try:
                x = K._coerce_(x[0])
                return K(zmod_poly_evaluate_horner(&self.x, int(x)))
            except TypeError:
                pass
        return Polynomial.__call__(self, *x, **kwds)

    def resultant(self, Polynomial_zmod_flint other):
        """
        Returns the resultant of self and other, which must lie in the same
        polynomial ring.

        INPUT:
            other -- a polynomial

        OUTPUT:
            an element of the base ring of the polynomial ring

        EXAMPLES:
            sage: R.<x> = GF(19)['x']
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: r = f.resultant(g); r
            11
            sage: r.parent() is GF(19)
            True
        """
        other = self.parent()._coerce_(other)
        res = zmod_poly_resultant(&(<Polynomial_template>self).x, &(<Polynomial_template>other).x)
        return self.parent().base_ring()(res)

    def small_roots(self, *args, **kwds):
        r"""
        See \code{sage.rings.polynomial.polynomial_modn_dense_ntl.small_roots}
        for the documentation of this function.

        EXAMPLE:
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
            n -- degree
            value -- coefficient

        WARNING: This could easily introduce subtle bugs, since \Sage
        assumes everywhere that polynomials are immutable.  It's OK to
        use this if you really know what you're doing.

        EXAMPLES:
            sage: R.<x> = GF(7)[]
            sage: f = (1+2*x)^2; f
            4*x^2 + 4*x + 1
            sage: f._unsafe_mutate(1, -5)
            sage: f
            4*x^2 + 2*x + 1
        """
        cdef cparent _parent = get_cparent((<Polynomial_template>self)._parent)

        n = int(n)
        value = self.base_ring()(value)
        if n >= 0:
            zmod_poly_set_coeff_ui(&self.x, n, int(value)%zmod_poly_modulus(&self.x))
        else:
            raise IndexError, "Polynomial coefficient index must be nonnegative."

    def _mul_zn_poly(self, other):
        r"""
        Returns the product of two polynomials using the zn_poly library.

        See \url{http://www.math.harvard.edu/~dmharvey/zn_poly/} for
        details on zn_poly.

        INPUT:
           self: Polynomial
           right: Polynomial (over same base ring as self)

        OUTPUT: Polynomial
           The product self*right.


        EXAMPLE:
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

        cdef Polynomial_zmod_flint r = <Polynomial_zmod_flint>PY_NEW(self.__class__)
        r._parent = (<Polynomial_zmod_flint>self)._parent

        cdef unsigned long p = self._parent.modulus()
        cdef unsigned long n1 = self.x.length
        cdef unsigned long n2 = _other.x.length

        cdef zn_mod_struct zn_mod

        zmod_poly_init2(&r.x, p, n1 + n2 -1 )

        zn_mod_init(&zn_mod, p)
        zn_array_mul(r.x.coeffs, self.x.coeffs, n1, _other.x.coeffs, n2, &zn_mod)
        r.x.length = n1 + n2 -1
        __zmod_poly_normalise(&r.x)
        zn_mod_clear(&zn_mod)
        return r

