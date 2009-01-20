"""
Univariate Polynomials over GF(2) via NTL's GF2X.

AUTHOR:
    -- Martin Albrecht (2008-10) initial implementation
"""


# We need to define this stuff before including the templating stuff
# to make sure the function get_cparent is found since it is used in
# 'polynomial_template.pxi'.

ctypedef long cparent
cdef cparent get_cparent(parent):
    return 0

# first we include the definitions
include "../../libs/ntl/ntl_GF2X_linkage.pxi"

# and then the interface
include "polynomial_template.pxi"

from sage.libs.all import pari

from sage.matrix.matrix_mod2_dense cimport mzd_write_bit, mzd_read_bit, Matrix_mod2_dense, word

cdef class Polynomial_GF2X(Polynomial_template):
    """
    Univariate Polynomials over GF(2) via NTL's GF2X.

    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: x^3 + x^2 + 1
        x^3 + x^2 + 1
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        Create a new univariate polynomials over GF(2).

        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: x^3 + x^2 + 1
            x^3 + x^2 + 1
        """
        try:
            if x.parent() is parent.base_ring() or x.parent() == parent.base_ring():
                x = int(x)
        except AttributeError:
            pass
        Polynomial_template.__init__(self, parent, x, check, is_gen, construct)

    def __getitem__(self, int i):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: f = x^3 + x^2 + 1; f
            x^3 + x^2 + 1
            sage: f[0]
            1
            sage: f[1]
            0
        """
        cdef long c = 0
        if 0 <= i < GF2X_NumBits(self.x):
            c = GF2_conv_to_long(GF2X_coeff(self.x, i))
        return self._parent.base_ring()(c)

    def _pari_(self, variable=None):
        """
        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: f = x^3 + x^2 + 1
            sage: pari(f)
            Mod(1, 2)*x^3 + Mod(1, 2)*x^2 + Mod(1, 2)
        """
        #TODO: put this in a superclass
        parent = self._parent
        if variable is None:
            variable = parent.variable_name()
        return pari(self.list()).Polrev(variable) * pari(1).Mod(2)

    def modular_composition(Polynomial_GF2X self, Polynomial_GF2X g, Polynomial_GF2X h, algorithm=None):
        """
        Compute f(g) % h.

        Both implementations Use Brent-Kung's Algorithm 2.1 (Fast
        Algorithms for Manipulation Formal Power Series, JACM 1978)

        INPUT:
            g -- a polynomial
            h -- a polynomial
            algorithm -- either 'native' or 'ntl' (default: 'native')

        EXAMPLE:
            sage: P.<x> = GF(2)[]
            sage: r = 279
            sage: f = x^r + x +1
            sage: g = x^r
            sage: g.modular_composition(g, f) == g(g) % f
            True

        AUTHORS:
             -- Paul Zimmermann (2008-10) initial implementation
             -- Martin Albrecht (2008-10) performance improvements
        """
        if g.parent() is not self.parent() or h.parent() is not self.parent():
            raise TypeError("Parents of the first three parameters must match.")

        from sage.misc.misc import verbose, cputime
        from sage.calculus.calculus import ceil
        from sage.matrix.constructor import Matrix
        from sage.rings.all import FiniteField as GF

        cdef Polynomial_GF2X res
        cdef GF2XModulus_c modulus
        GF2XModulus_build(modulus, (<Polynomial_GF2X>h).x)

        res = <Polynomial_GF2X>PY_NEW(Polynomial_GF2X)
        res._parent = self._parent

        if algorithm == "ntl":
            t = cputime()
            _sig_on
            GF2X_CompMod(res.x, self.x, g.x, modulus)
            _sig_off
            verbose("NTL %5.3f s"%cputime(t),level=1)
            return res

        cdef Py_ssize_t i, j, k, l, n, maxlength
        cdef Matrix_mod2_dense F, G, H

        cdef GF2X_c _f = (<Polynomial_GF2X>self).x
        cdef GF2X_c _g = (<Polynomial_GF2X>g).x
        cdef GF2X_c gpow, g2, tt
        GF2X_conv_long(gpow, 1)

        maxlength = GF2X_NumBits(_f)

        t = cputime()

        n = h.degree()

        k = ceil(Integer(n+1).sqrt_approx())
        l = ceil((self.degree() + 1) / k)

        # we store all matrices transposed for performance reasons
        G = <Matrix_mod2_dense>Matrix(GF(2), k, n)

        # first compute g^j mod h, 2 <= j < k
        # first deal with j=0
        for i from 0 <= i < GF2X_NumBits(gpow):
            mzd_write_bit(G._entries, 0, i, GF2_conv_to_long(GF2X_coeff(gpow, i)))
        # precompute g^2
        GF2X_SqrMod_pre(g2, _g, modulus)
        gpow = _g
        for j in range(1, k, 2):
            if j > 1:
                GF2X_MulMod_pre(gpow, gpow, g2, modulus) # gpow = g^j
            for i from 0 <= i < GF2X_NumBits(gpow):
                mzd_write_bit(G._entries, j, i, GF2_conv_to_long(GF2X_coeff(gpow, i)))
            # we now process 2j, 4j, 8j, ... by squaring each time
            if 2*j < k:
                tt = gpow
                jj = j
                while 2*jj < k:
                   GF2X_SqrMod_pre(tt, tt, modulus)
                   jj = 2*jj
                   for i from 0 <= i < GF2X_NumBits(tt):
                       mzd_write_bit(G._entries, jj, i, GF2_conv_to_long(GF2X_coeff(tt, i)))
        # we need that gpow = g^k at the end
        if k % 2 == 1: # k is odd, last j is k-2
            GF2X_MulMod_pre(gpow, gpow, g2, modulus)
        else:          # k is even, last j is k-1
            GF2X_MulMod_pre(gpow, gpow, _g, modulus)
        verbose("G %d x %d %5.3f s"%(G.nrows(), G.ncols(),cputime(t)),level=1)

        # split f in chunks of degree < k
        t = cputime()
        F = <Matrix_mod2_dense>Matrix(GF(2), l, k)
        for j in range(0, l):
            if j*k+k <= maxlength:
                for i from j*k <= i < j*k+k:
                    mzd_write_bit(F._entries, j, i-j*k, GF2_conv_to_long(GF2X_coeff(_f, i)))
            else:
                for i from j*k <= i < maxlength:
                    mzd_write_bit(F._entries, j, i-j*k, GF2_conv_to_long(GF2X_coeff(_f, i)))

        verbose("F %d x %d %5.3f s"%(F.nrows(), F.ncols(), cputime(t)),level=1)

        t = cputime()
        H = <Matrix_mod2_dense>(F * G)
        verbose("H %d x %d %5.3f s"%(H.nrows(), H.ncols(), cputime(t)),level=1)

        t = cputime()
        # H is a n x l matrix now H[i,j] = sum(G[i,m]*F[m,j],
        # m=0..k-1) = sum(g^m[i] * f[j*k+m], m=0..k-1) where g^m[i] is
        # the coefficient of degree i in g^m and f[j*k+m] is the
        # coefficient of degree j*k+m in f thus f[j*k+m]*g^m[i] should
        # be multiplied by g^(j*k) gpow = (g^k) % h

        GF2X_conv_long(res.x, 0)
        j = l - 1
        while j >= 0:
            #res = (res * gpow) % h
            GF2X_MulMod_pre(res.x, res.x, gpow, modulus)

            # res = res + parent([H[j,i] for i in range(0,n)])
            GF2X_conv_long(tt, 0)
            for i from 0<= i < n:
                GF2X_SetCoeff_long(tt, i, mzd_read_bit(H._entries, j, i))
            GF2X_add(res.x, res.x, tt)
            j = j - 1

        verbose("Res %5.3f s"%cputime(t),level=1)
        return res
