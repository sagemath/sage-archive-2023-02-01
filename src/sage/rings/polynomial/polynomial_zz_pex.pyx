"""
Univariate Polynomials over GF(p^e) via NTL's ZZ_pEX.

AUTHOR:
- Yann Laigle-Chapuy (2010-01) initial implementation
"""

from sage.rings.integer_ring import ZZ
from sage.rings.integer_ring cimport IntegerRing_class

from sage.libs.ntl.ntl_ZZ_pEContext cimport ntl_ZZ_pEContext_class
from sage.libs.ntl.ntl_ZZ_pEContext_decl cimport ZZ_pEContext_c
from sage.libs.ntl.ntl_ZZ_pE_decl cimport ZZ_pE_to_PyString
from sage.libs.ntl.ntl_ZZ_pE_decl cimport ZZ_pE_to_ZZ_pX
from sage.libs.ntl.ntl_ZZ_pX_decl cimport ZZ_pX_to_PyString
from sage.libs.ntl.ntl_ZZ_pX_decl cimport ZZ_pX_deg, ZZ_pX_coeff
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX
from sage.libs.ntl.ntl_ZZ_p_decl cimport ZZ_p_to_PyString
from sage.libs.ntl.ntl_ZZ_p_decl cimport ZZ_p_rep
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class

# We need to define this stuff before including the templating stuff
# to make sure the function get_cparent is found since it is used in
# 'polynomial_template.pxi'.

cdef cparent get_cparent(parent) except? NULL:
    if parent is None:
        return NULL
    cdef ntl_ZZ_pEContext_class pec
    try:
        pec = parent._modulus
    except AttributeError:
        return NULL
    return &(pec.ptrs)

# first we include the definitions
include "sage/libs/ntl/ntl_ZZ_pEX_linkage.pxi"

# and then the interface
include "polynomial_template.pxi"

from sage.libs.all import pari
from sage.libs.ntl.ntl_ZZ_pE cimport ntl_ZZ_pE

cdef inline ZZ_pE_c_to_list(ZZ_pE_c x):
    cdef list L = []
    cdef ZZ_pX_c c_pX
    cdef ZZ_p_c c_p
    cdef ZZ_c c_c

    c_pX = ZZ_pE_to_ZZ_pX(x)
    d = ZZ_pX_deg(c_pX)
    if d>=0:
        for 0 <= j <= d:
            c_p = ZZ_pX_coeff(c_pX, j)
            c_c = ZZ_p_rep(c_p)
            L.append((<IntegerRing_class>ZZ)._coerce_ZZ(&c_c))
    return L


cdef class Polynomial_ZZ_pEX(Polynomial_template):
    """
    Univariate Polynomials over GF(p^n) via NTL's ZZ_pEX.

    EXAMPLE::

        sage: K.<a>=GF(next_prime(2**60)**3)
        sage: R.<x> = PolynomialRing(K,implementation='NTL')
        sage: (x^3 + a*x^2 + 1) * (x + a)
        x^4 + 2*a*x^3 + a^2*x^2 + x + a
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        Create a new univariate polynomials over GF(p^n).

        EXAMPLE::

            sage: K.<a>=GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K,implementation='NTL')
            sage: x^2+a
            x^2 + a

        TESTS:

        The following tests against a bug that was fixed in :trac:`9944`.
        With the ring definition above, we now have::

            sage: R([3,'1234'])
            1234*x + 3
            sage: R([3,'12e34'])
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x (=12e34) to an integer
            sage: R([3,x])
            Traceback (most recent call last):
            ...
            TypeError: not a constant polynomial

        Check that NTL contexts are correctly restored and that
        :trac:`9524` has been fixed::

            sage: x = polygen(GF(9, 'a'))
            sage: x = polygen(GF(49, 'a'))
            sage: -x
            6*x
            sage: 5*x
            5*x

        Check that :trac:`11239` is fixed::

            sage: Fq.<a> = GF(2^4); Fqq.<b> = GF(3^7)
            sage: PFq.<x> = Fq[]; PFqq.<y> = Fqq[]
            sage: f = x^3 + (a^3 + 1)*x
            sage: sage.rings.polynomial.polynomial_zz_pex.Polynomial_ZZ_pEX(PFqq, f)
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce from a finite field other than the prime subfield
        """
        cdef ntl_ZZ_pE d
        try:
            if (x.parent() is parent.base_ring()) or (x.parent() == parent.base_ring()):
                Polynomial.__init__(self, parent, is_gen=is_gen)
                (<Polynomial_template>self)._cparent = get_cparent(parent)
                celement_construct(&self.x, (<Polynomial_template>self)._cparent)
                d = parent._modulus.ZZ_pE(list(x.polynomial()))
                ZZ_pEX_SetCoeff(self.x, 0, d.x)
                return
        except AttributeError:
            pass

        if isinstance(x, Polynomial):
            x = x.list()

        if PY_TYPE_CHECK(x, list) or PY_TYPE_CHECK(x, tuple):
            Polynomial.__init__(self, parent, is_gen=is_gen)
            (<Polynomial_template>self)._cparent = get_cparent(parent)
            celement_construct(&self.x, (<Polynomial_template>self)._cparent)
            K = parent.base_ring()
            for i,e in enumerate(x):
                # self(x) is supposed to be a conversion,
                # not necessarily a coercion. So, we must
                # not do K.coerce(e) but K(e).
                e = K(e)
                d = parent._modulus.ZZ_pE(list(e.polynomial()))
                ZZ_pEX_SetCoeff(self.x, i, d.x)
            return

        Polynomial_template.__init__(self, parent, x, check, is_gen, construct)

    def __getitem__(self,i):
        """
        EXAMPLES::

            sage: K.<a>=GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K,implementation='NTL')
            sage: f = x^3+(2*a+1)*x+a
            sage: f[0]
            a
            sage: f[1]
            2*a + 1
            sage: f[2]
            0
            sage: f[1:4]
            x^3 + (2*a + 1)*x
            sage: f[-5:50] == f
            True
        """
        cdef ZZ_pE_c c_pE
        cdef Polynomial_template r
        if isinstance(i, slice):
            start, stop = i.start, i.stop
            if start < 0:
                start = 0
            if stop > celement_len(&self.x, (<Polynomial_template>self)._cparent) or stop is None:
                stop = celement_len(&self.x, (<Polynomial_template>self)._cparent)
            x = (<Polynomial_template>self)._parent.gen()
            v = [self[t] for t from start <= t < stop]

            r = <Polynomial_template>PY_NEW(self.__class__)
            Polynomial_template.__init__(r, (<Polynomial_template>self)._parent, v)
            return r << start
        else:
            self._parent._modulus.restore()
            c_pE = ZZ_pEX_coeff(self.x, i)

            K = self._parent.base_ring()
            return K(ZZ_pE_c_to_list(c_pE))

    def list(self):
        """
        Returs the list of coefficients.

        EXAMPLE::

            sage: K.<a> = GF(5^3)
            sage: P = PolynomialRing(K, 'x')
            sage: f = P.random_element(100)
            sage: f.list() == [f[i] for i in  range(f.degree()+1)]
            True
            sage: P.0.list()
            [0, 1]

        """
        cdef Py_ssize_t i

        self._parent._modulus.restore()

        K = self._parent.base_ring()
        return [K(ZZ_pE_c_to_list(ZZ_pEX_coeff(self.x, i))) for i in range(celement_len(&self.x, (<Polynomial_template>self)._cparent))]

    cpdef ModuleElement _rmul_(self, RingElement left):
        """
        EXAMPLE::
            sage: K.<a>=GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K,implementation='NTL')
            sage: (2*a+1)*x # indirect doctest
            (2*a + 1)*x
        """
        cdef ntl_ZZ_pE d
        cdef Polynomial_ZZ_pEX r
        r = PY_NEW(Polynomial_ZZ_pEX)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent
        d = self._parent._modulus.ZZ_pE(list(left.polynomial()))
        ZZ_pEX_mul_ZZ_pE(r.x, self.x, d.x)
        return r

    def __call__(self, *x, **kwds):
        """
        Evaluate polynomial at `a`.

        EXAMPLE::

            sage: K.<u>=GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K,implementation='NTL')
            sage: P = (x-u)*(x+u+1)
            sage: P(u)
            0
            sage: P(u+1)
            2*u + 2

        TESTS:

        The following was fixed at trac ticket #10475::

            sage: F.<x> = GF(4)
            sage: P.<y> = F[]
            sage: p = y^4 + x*y^3 + y^2 + (x + 1)*y + x + 1
            sage: SR(p)      #indirect doctest
            (((y + x)*y + 1)*y + x + 1)*y + x + 1
            sage: p(2)
            x + 1
            sage: p(y=2)
            x + 1
            sage: p(3,y=2)
            Traceback (most recent call last):
            ...
            TypeError: <type 'sage.rings.polynomial.polynomial_zz_pex.Polynomial_ZZ_pEX'>__call__() takes exactly 1 argument
            sage: p(x=2)
            Traceback (most recent call last):
            ...
            TypeError: <type 'sage.rings.polynomial.polynomial_zz_pex.Polynomial_ZZ_pEX'>__call__() accepts no named argument except 'y'

        Check that polynomial evaluation works when using logarithmic
        representation of finite field elements (:trac:`16383`)::

            sage: for i in xrange(10):
            ....:     F = FiniteField(random_prime(15) ** ZZ.random_element(2, 5), 'a', repr='log')
            ....:     b = F.random_element()
            ....:     P = PolynomialRing(F, 'x')
            ....:     f = P.random_element(8)
            ....:     assert f(b) == sum(c * b^i for i, c in enumerate(f))

        """
        cdef ntl_ZZ_pE _a
        cdef ZZ_pE_c c_b

        K = self._parent.base_ring()

        if kwds:
            if x:
                raise TypeError("%s__call__() takes exactly 1 argument"%type(self))
            try:
                x = [kwds.pop(self.variable_name())]
            except KeyError:
                pass
        if kwds:
            raise TypeError("%s__call__() accepts no named argument except '%s'"%(type(self),self.variable_name()))
        if len(x)!=1:
            raise TypeError("%s__call__() takes exactly 1 positional argument"%type(self))

        a = x[0]
        try:
            if a.parent() is not K:
                a = K.coerce(a)
        except (TypeError, AttributeError, NotImplementedError):
            return Polynomial.__call__(self, a)

        _a = self._parent._modulus.ZZ_pE(list(a.polynomial()))
        ZZ_pEX_eval(c_b, self.x, _a.x)
        return K(ZZ_pE_c_to_list(c_b))

    def resultant(self, other):
        """
        Returns the resultant of self and other, which must lie in the same
        polynomial ring.

        INPUT:

        :argument other: a polynomial

        OUTPUT: an element of the base ring of the polynomial ring

        EXAMPLES::

            sage: K.<a>=GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K,implementation='NTL')
            sage: f=(x-a)*(x-a**2)*(x+1)
            sage: g=(x-a**3)*(x-a**4)*(x+a)
            sage: r = f.resultant(g)
            sage: r == prod(u-v for (u,eu) in f.roots() for (v,ev) in g.roots())
            True
        """
        cdef ZZ_pE_c r
        self._parent._modulus.restore()

        if other.parent() is not self._parent:
            other = self._parent.coerce(other)

        ZZ_pEX_resultant(r, self.x, (<Polynomial_ZZ_pEX>other).x)

        K = self._parent.base_ring()
        return K(K.polynomial_ring()(ZZ_pE_c_to_list(r)))

    def is_irreducible(self, algorithm="fast_when_false", iter=1):
        """
        Returns `True` precisely when self is irreducible over its base ring.

        INPUT:

        :argument algorithm: a string (default "fast_when_false"),
            there are 3 available algorithms:
            "fast_when_true", "fast_when_false" and "probabilistic".
        :argument iter: (default: 1) if the algorithm is "probabilistic"
            defines the number of iterations. The error probability is bounded
            by `q**-iter` for polynomials in `GF(q)[x]`.

        EXAMPLES::

            sage: K.<a>=GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K,implementation='NTL')
            sage: P = x^3+(2-a)*x+1
            sage: P.is_irreducible(algorithm="fast_when_false")
            True
            sage: P.is_irreducible(algorithm="fast_when_true")
            True
            sage: P.is_irreducible(algorithm="probabilistic")
            True
            sage: Q = (x^2+a)*(x+a^3)
            sage: Q.is_irreducible(algorithm="fast_when_false")
            False
            sage: Q.is_irreducible(algorithm="fast_when_true")
            False
            sage: Q.is_irreducible(algorithm="probabilistic")
            False
            """
        self._parent._modulus.restore()
        if algorithm=="fast_when_false":
            sig_on()
            res = ZZ_pEX_IterIrredTest(self.x)
            sig_off()
        elif algorithm=="fast_when_true":
            sig_on()
            res = ZZ_pEX_DetIrredTest(self.x)
            sig_off()
        elif algorithm=="probabilistic":
            sig_on()
            res = ZZ_pEX_ProbIrredTest(self.x, iter)
            sig_off()
        else:
            raise ValueError("unknown algorithm")
        return res != 0

    cdef int _cmp_c_impl(left,Element right) except -2:
        """
        EXAMPLE::

            sage: K.<a>=GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K,implementation='NTL')
            sage: P1 = (a**2+a+1)*x^2+a*x+1
            sage: P2 = (     a+1)*x^2+a*x+1
            sage: P1 < P2 # indirect doctests
            False

        TEST::

            sage: P3 = (a**2+a+1)*x^2+  x+1
            sage: P4 =                  x+1
            sage: P1 < P3
            False
            sage: P1 < P4
            False
            sage: P1 > P2
            True
            sage: P1 > P3
            True
            sage: P1 > P4
            True
        """
        cdef long ld, rd, i

        left._parent._modulus.restore()
        ld = left.degree()
        rd = right.degree()
        if ld < rd: return -1
        if ld > rd: return 1
        # degrees are equal
        for i in range(ld,-1,-1):
            li = left[i]
            ri = right[i]
            t = li.__cmp__(ri)
            if t != 0:
                return t
        return 0

    def shift(self, int n):
        """
        EXAMPLE::

            sage: K.<a>=GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K,implementation='NTL')
            sage: f = x^3 + x^2 + 1
            sage: f.shift(1)
            x^4 + x^3 + x
            sage: f.shift(-1)
            x^2 + x
        """
        self._parent._modulus.restore()
        cdef Polynomial_ZZ_pEX r
        r = PY_NEW(Polynomial_ZZ_pEX)
        celement_construct(&r.x, (<Polynomial_template>self)._cparent)
        r._parent = (<Polynomial_template>self)._parent
        r._cparent = (<Polynomial_template>self)._cparent
        ZZ_pEX_LeftShift(r.x, self.x, n)
        return r

    def __lshift__(self, int n):
        """
        EXAMPLE::

            sage: K.<a>=GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K,implementation='NTL')
            sage: f = x^3 + x^2 + 1
            sage: f << 1
            x^4 + x^3 + x
            sage: f << -1
            x^2 + x
        """
        return self.shift(n)

    def __rshift__(self, int n):
        """
        EXAMPLE::

            sage: K.<a>=GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K,implementation='NTL')
            sage: f = x^3 + x^2 + 1
            sage: f >> 1
            x^2 + x
            sage: f >> -1
            x^4 + x^3 + x
        """
        return self.shift(-n)

