
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

# We need to define this stuff before including the templating stuff
# to make sure the function get_cparent is found since it is used in
# 'polynomial_template.pxi'.

ctypedef ZZ_pEContext_c *cparent

cdef cparent get_cparent(parent):
    if parent is None:
        return NULL
    cdef ntl_ZZ_pEContext_class c
    c = parent._PolynomialRing_field__modulus
    return &(c.x)

# first we include the definitions
include "../../libs/ntl/ntl_ZZ_pEX_linkage.pxi"

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
        Create a new univariate polynomials over GF(2).

        EXAMPLE::

            sage: K.<a>=GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K,implementation='NTL')
            sage: x^2+a
            x^2 + a
        """
        cdef cparent _parent
        cdef ntl_ZZ_pE d
        try:
            if (x.parent() is parent.base_ring()) or (x.parent() == parent.base_ring()):
                _parent = get_cparent(parent)
                Polynomial.__init__(self, parent, is_gen=is_gen)
                celement_construct(&self.x, _parent)
                d = parent._PolynomialRing_field__modulus.ZZ_pE(list(x.polynomial()))
                ZZ_pEX_SetCoeff(self.x, 0, d.x)
                return
        except AttributeError:
            pass

        if PY_TYPE_CHECK(x, list) or PY_TYPE_CHECK(x, tuple):
            _parent = get_cparent(parent)
            Polynomial.__init__(self, parent, is_gen=is_gen)
            celement_construct(&self.x, _parent)
            K = parent.base_ring()
            for i,e in enumerate(x):
                if not hasattr(e,'polynomial'):
                    try:
                        e = K._coerce_(e)
                    except:
                        TypeError("unable to coerce this value to the base ring")
                d = parent._PolynomialRing_field__modulus.ZZ_pE(list(e.polynomial()))
                ZZ_pEX_SetCoeff(self.x, i, d.x)
            return

        Polynomial_template.__init__(self, parent, x, check, is_gen, construct)

    def __getitem__(self,i):
        """
        EXAMPLE::

            sage: K.<a>=GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K,implementation='NTL')
            sage: f = x^3+(2*a+1)*x+a
            sage: f[0]
            a
            sage: f[1]
            2*a + 1
            sage: f[2]
            0
        """
        cdef list L = []
        cdef ZZ_pE_c c_pE
        cdef cparent _parent

        _parent = get_cparent(self._parent)
        _parent[0].restore()
        c_pE = ZZ_pEX_coeff(self.x, i)

        K = self._parent.base_ring()
        return K(K.polynomial_ring()(ZZ_pE_c_to_list(c_pE)))
        return ZZ_pE_to_PyString(&c_pE)

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
        celement_construct(&r.x, get_cparent(self._parent))
        r._parent = self._parent
        d = self._parent._PolynomialRing_field__modulus.ZZ_pE(list(left.polynomial()))
        ZZ_pEX_mul_ZZ_pE(r.x, self.x, d.x)
        return r

    def __call__(self, a):
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
        """
        cdef ntl_ZZ_pE _a
        cdef ZZ_pE_c c_b

        K = self._parent.base_ring()

        try:
            if a.parent() is not K:
                a = coerce(K,a)
        except:
            return Polynomial.__call__(self, a)

        _a = self._parent._PolynomialRing_field__modulus.ZZ_pE(list(a.polynomial()))
        ZZ_pEX_eval(c_b, self.x, _a.x)
        g = K.gen()
        res = K(0)
        m = K(1)
        for c in ZZ_pE_c_to_list(c_b):
            res += m*c
            m *= g
        return res

    def resultant(self, Polynomial_ZZ_pEX other):
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
        self._parent._PolynomialRing_field__modulus.restore()

        if other._parent is not self._parent:
            other = self._parent._coerce_(other)

        ZZ_pEX_resultant(r, self.x, other.x)

        K = self._parent.base_ring()
        return K(K.polynomial_ring()(ZZ_pE_c_to_list(r)))

    def is_irreducible(self, algorithm="fast_when_false"):
        """
        Returns `True` precisely when self is irreducible over its base ring.

        INPUT:

        :argument algorithm: a string (default "fast_when_false"),
            there are 3 available algorithms:
            "fast_when_true", "fast_when_false" and "probabilistic".

        EXAMPLES::

            sage: K.<a>=GF(next_prime(2**60)**3)
            sage: R.<x> = PolynomialRing(K,implementation='NTL')
            sage: P = x^3+(a+3)*x+1
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
        self._parent._PolynomialRing_field__modulus.restore()
        if algorithm=="fast_when_false":
            _sig_on
            res = ZZ_pEX_IterIrredTest(self.x)
            _sig_off
        elif algorithm=="fast_when_true":
            _sig_on
            res = ZZ_pEX_DetIrredTest(self.x)
            _sig_off
        elif algorithm=="probabilistic":
            _sig_on
            res = ZZ_pEX_ProbIrredTest(self.x, 1)
            _sig_off
        else:
            raise ValueError("unknown algorithm")
        return res != 0

    cdef int _cmp_c_impl(left,Element right) except -2:
        left._parent._PolynomialRing_field__modulus.restore()
        ld = left.degree()
        rd = right.degree()
        if ld < rd: return -1
        if ld > rd: return 1
        # degrees are equal
        cdef int i
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
        self._parent._PolynomialRing_field__modulus.restore()
        cdef Polynomial_ZZ_pEX r
        r = PY_NEW(Polynomial_ZZ_pEX)
        celement_construct(&r.x, get_cparent(self._parent))
        r._parent = self._parent
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

