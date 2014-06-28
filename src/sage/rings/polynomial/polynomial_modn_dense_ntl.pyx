"""
Dense univariate polynomials over  `\ZZ/n\ZZ`, implemented using NTL.

This implementation is generally slower than the FLINT implementation in
:mod:`~sage.rings.polynomial.polynomial_zmod_flint`, so we use FLINT by
default when the modulus is small enough; but NTL does not require that `n` be
``int``-sized, so we use it as default when `n` is too large for FLINT.

Note that the classes :class:`Polynomial_dense_modn_ntl_zz` and
:class:`Polynomial_dense_modn_ntl_ZZ` are different; the former is limited to
moduli less than a certain bound, while the latter supports arbitrarily large
moduli.

AUTHORS:

- Robert Bradshaw: Split off from polynomial_element_generic.py (2007-09)
- Robert Bradshaw: Major rewrite to use NTL directly (2007-09)
"""

################################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
################################################################################

from sage.rings.polynomial.polynomial_element import is_Polynomial, Polynomial_generic_dense

from sage.libs.all import pari, pari_gen

from sage.libs.ntl.all import ZZ as ntl_ZZ, ZZX, zero_ZZX, ZZ_p, ZZ_pX
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod import IntegerMod_abstract

from sage.rings.fraction_field_element import FractionFieldElement
import sage.rings.polynomial.polynomial_ring

from sage.rings.infinity import infinity

import polynomial_singular_interface
from sage.interfaces.all import singular as singular_default

from sage.structure.element import generic_power, canonical_coercion, bin_op, coerce_binop

from sage.libs.ntl.ntl_ZZ_p_decl cimport *, ZZ_p_c
from sage.libs.ntl.ntl_lzz_p_decl cimport *, zz_p_c
from sage.libs.ntl.ntl_lzz_pX_decl cimport *, zz_pX_c, zz_pX_Modulus_c
from sage.libs.ntl.ntl_ZZ_pX_decl cimport *, ZZ_pX_c, ZZ_pX_Modulus_c

def make_element(parent, args):
    return parent(*args)

include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"
include "sage/ext/cdefs.pxi"

zz_p_max = NTL_SP_BOUND

cdef class Polynomial_dense_mod_n(Polynomial):
    """
    A dense polynomial over the integers modulo n, where n is composite, with
    the underlying arithmetic done using NTL.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(Integers(16))
        sage: f = x^3 - x + 17
        sage: f^2
        x^6 + 14*x^4 + 2*x^3 + x^2 + 14*x + 1

        sage: loads(f.dumps()) == f
        True

        sage: R.<x> = Integers(100)[]
        sage: p = 3*x
        sage: q = 7*x
        sage: p+q
        10*x
        sage: R.<x> = Integers(8)[]
        sage: parent(p)
        Univariate Polynomial Ring in x over Ring of integers modulo 100
        sage: p + q
        10*x
        sage: R({10:-1})
        7*x^10

    """
    def __init__(self, parent, x=None, check=True,
                 is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)
        cdef Polynomial_dense_mod_n numer, denom

        if construct:
            if isinstance(x, ZZ_pX):
                self.__poly = x
                return
            self.__poly = ZZ_pX(x, parent.modulus())
            return

        self.__poly = ZZ_pX([], parent.modulus())

        if x is None:
            return         # leave initialized to 0 polynomial.

        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                self.__poly = (<Polynomial_dense_modn_ntl_zz>x).__poly.__copy__()
                return
            else:
                R = parent.base_ring()
                x = [ZZ(R(a)) for a in x.list()]
                check = False

        elif isinstance(x, dict):
            R = parent.base_ring()
            x = self._dict_to_list(x, R(0))


        elif isinstance(x, ZZX):
            self.__poly = x.copy()
            return

        elif isinstance(x, pari_gen):
            x = [ZZ(w) for w in x.list()]
            check = False

        elif isinstance(x, FractionFieldElement) and \
                 isinstance(x.numerator(), Polynomial_dense_mod_n):
            if x.denominator().is_unit():
                numer = x.numerator()
                denom = x.denominator().inverse_of_unit()
                x = numer.__poly * denom.__poly
                check = False
            else:
                raise TypeError("Denominator not a unit.")

        elif not isinstance(x, list) and not isinstance(x, tuple):
            x = [x]   # constant polynomials

        if check:
            R = parent.base_ring()
            x = [ZZ(R(a)) for a in x]

        self.__poly = ZZ_pX(x, parent.modulus())


    def __reduce__(self):
        return make_element, (self.parent(), (self.list(), False, self.is_gen()))

    def int_list(self):
        return eval(str(self.__poly).replace(' ',','))

    def _pari_(self, variable=None):
        """
        EXAMPLES::

            sage: t = PolynomialRing(IntegerModRing(17),"t").gen()
            sage: f = t^3 + 3*t - 17
            sage: pari(f)
            Mod(1, 17)*t^3 + Mod(3, 17)*t
        """
        if variable is None:
            variable = self.parent().variable_name()
        return pari(self.int_list()).Polrev(variable) * \
               pari(1).Mod(self.parent().base_ring().order())

    def ntl_ZZ_pX(self):
        r"""
        Return underlying NTL representation of this polynomial.  Additional
        ''bonus'' functionality is available through this function.

        .. warning::

            You must call ``ntl.set_modulus(ntl.ZZ(n))`` before doing
            arithmetic with this object!
        """
        return self.__poly

    def __getitem__(self, n):
        """
        Returns coefficient of the monomial of degree `n` if `n` is an integer,
        returns the monomials of self of degree in slice `n` if `n` is a slice.

        EXAMPLES::

            sage: R.<x> = Integers(100)[]
            sage: from sage.rings.polynomial.polynomial_modn_dense_ntl import Polynomial_dense_mod_n
            sage: f = Polynomial_dense_mod_n(R,[5,10,13,1,4]); f
            4*x^4 + x^3 + 13*x^2 + 10*x + 5
            sage: f[2]
            13
            sage: f[1:3]
            13*x^2 + 10*x
        """
        if isinstance(n, slice):
            start, stop = n.start, n.stop
            R = self.base_ring()
            if start < 0:
                start = 0
            if stop > self.__poly.degree()+1 or stop is None:
                stop = self.__poly.degree()+1
            v = [R(self.__poly[k]._sage_()) for k in range(start,stop)]
            return self.parent()([0]*int(start) + v)
        else:
            return self.parent().base_ring()(self.__poly[n]._sage_())

    def _unsafe_mutate(self, n, value):
        n = int(n)
        if n < 0:
            raise IndexError("n must be >= 0")
        self.__poly[n] = int(value)

    def _pow(self, n):
        n = int(n)

        if self.degree() <= 0:
            return self.parent()(self[0]**n)
        if n < 0:
            return (~self)**(-n)
        return self.parent()(self.__poly**n, construct=True)

    cpdef ModuleElement _add_(self, ModuleElement right):
        return self.parent()(self.__poly + (<Polynomial_dense_mod_n>right).__poly, construct=True)

    cpdef RingElement _mul_(self, RingElement right):
        """
        EXAMPLES::

            sage: x = PolynomialRing(Integers(100), 'x').0
            sage: (x - 2)*(x^2 - 8*x + 16)
            x^3 + 90*x^2 + 32*x + 68
        """
        return self.parent()(self.__poly * (<Polynomial_dense_mod_n>right).__poly, construct=True)

    cpdef ModuleElement _rmul_(self, RingElement c):
        try:
            return self.parent()(ZZ_pX([c], self.parent().modulus()) * self.__poly, construct=True)
        except RuntimeError as msg: # should this really be a TypeError
            raise TypeError(msg)

    cpdef ModuleElement _lmul_(self, RingElement c):
        try:
            return self.parent()(ZZ_pX([c], self.parent().modulus()) * self.__poly, construct=True)
        except RuntimeError as msg: # should this really be a TypeError
            raise TypeError(msg)

    @coerce_binop
    def quo_rem(self, right):
        """
        Returns a tuple (quotient, remainder) where self = quotient*other +
        remainder.
        """
        v = self.__poly.quo_rem((<Polynomial_dense_mod_n>right).__poly)
        P = self.parent()
        return (P(v[0], construct=True), P(v[1], construct=True) )

    def shift(self, n):
        r"""
        Returns this polynomial multiplied by the power `x^n`. If `n` is negative,
        terms below `x^n` will be discarded. Does not change this polynomial.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(12345678901234567890))
            sage: p = x^2 + 2*x + 4
            sage: p.shift(0)
             x^2 + 2*x + 4
            sage: p.shift(-1)
             x + 2
            sage: p.shift(-5)
             0
            sage: p.shift(2)
             x^4 + 2*x^3 + 4*x^2

        TESTS::

            sage: p = R(0)
            sage: p.shift(3).is_zero()
            True
            sage: p.shift(-3).is_zero()
            True

        AUTHOR:

        - David Harvey (2006-08-06)
        """
        if n == 0 or self.degree() < 0:
            return self
        return self.parent()(self.__poly.left_shift(n),
                             construct=True)

    cpdef ModuleElement _sub_(self, ModuleElement right):
        return self.parent()(self.__poly - (<Polynomial_dense_mod_n>right).__poly, construct=True)

    def __floordiv__(self, right):
        q, _ = self.quo_rem(right)
        return q


    def degree(self, gen=None):
        """
        Return the degree of this polynomial.  The zero polynomial
        has degree -1.
        """
        return max(self.__poly.degree(), -1)

    def list(self):
        """
        Return a new copy of the list of the underlying
        elements of self.

        EXAMPLES::

            sage: _.<x> = Integers(100)[]
            sage: f = x^3 + 3*x - 17
            sage: f.list()
            [83, 3, 0, 1]
        """
        R = self.base_ring()
        return [R(x) for x in self.int_list()]

    def ntl_set_directly(self, v):
        r"""
        Set the value of this polynomial directly from a vector or string.

        Polynomials over the integers modulo n are stored internally using
        NTL's ``ZZ_pX`` class.  Use this function to set the value of this
        polynomial using the NTL constructor, which is potentially *very* fast.
        The input v is either a vector of ints or a string of the form ``[ n1
        n2 n3 ... ]`` where the ni are integers and there are no commas between
        them. The optimal input format is the string format, since that's what
        NTL uses by default.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(100))
            sage: from sage.rings.polynomial.polynomial_modn_dense_ntl import Polynomial_dense_mod_n as poly_modn_dense
            sage: poly_modn_dense(R, ([1,-2,3]))
            3*x^2 + 98*x + 1
            sage: f = poly_modn_dense(R, 0)
            sage: f.ntl_set_directly([1,-2,3])
            sage: f
            3*x^2 + 98*x + 1
            sage: f.ntl_set_directly('[1 -2 3 4]')
            sage: f
            4*x^3 + 3*x^2 + 98*x + 1
        """
        if self.is_gen():
            raise TypeError("Cannot change the value of the generator.")
        self.__poly = ZZ_pX(v, self.parent().modulus())

    # Polynomial_singular_repr stuff, copied due to lack of multiple inheritance
    def _singular_(self, singular=singular_default, have_ring=False, force=False):
        if not have_ring:
            self.parent()._singular_(singular,force=force).set_ring() # this is expensive
        if self.__singular is not None:
            try:
                self.__singular._check_valid()
                if self.__singular.parent() is singular:
                    return self.__singular
            except (AttributeError, ValueError):
                pass
        return self._singular_init_(singular, have_ring=have_ring)

    def _singular_init_(self, singular=singular_default, have_ring=False, force=False):
        if not have_ring:
            self.parent()._singular_(singular,force=force).set_ring() # this is expensive
        self.__singular = singular(str(self))
        return self.__singular

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
        return small_roots(self, *args, **kwds)

def small_roots(self, X=None, beta=1.0, epsilon=None, **kwds):
    r"""
    Let `N` be the characteristic of the base ring this polynomial
    is defined over: ``N = self.base_ring().characteristic()``.
    This method returns small roots of this polynomial modulo some
    factor `b` of `N` with the constraint that `b >= N^\beta`.
    Small in this context means that if `x` is a root of `f`
    modulo `b` then `|x| < X`. This `X` is either provided by the
    user or the maximum `X` is chosen such that this algorithm
    terminates in polynomial time. If `X` is chosen automatically
    it is `X = ceil(1/2 N^{\beta^2/\delta - \epsilon})`.
    The algorithm may also return some roots which are larger than `X`.
    'This algorithm' in this context means Coppersmith's algorithm for finding
    small roots using the LLL algorithm. The implementation of this algorithm
    follows Alexander May's PhD thesis referenced below.

    INPUT:

    - ``X`` -- an absolute bound for the root (default: see above)
    - ``beta`` -- compute a root mod `b` where `b` is a factor of `N` and `b
      \ge N^\beta`. (Default: 1.0, so `b = N`.)
    - ``epsilon`` -- the parameter `\epsilon` described above. (Default: `\beta/8`)
    - ``**kwds`` -- passed through to method :meth:`Matrix_integer_dense.LLL() <sage.matrix.matrix_integer_dense.Matrix_integer_dense.LLL>`.

    EXAMPLES:

    First consider a small example::

        sage: N = 10001
        sage: K = Zmod(10001)
        sage: P.<x> = PolynomialRing(K)
        sage: f = x^3 + 10*x^2 + 5000*x - 222

    This polynomial has no roots without modular reduction (i.e. over `\ZZ`)::

        sage: f.change_ring(ZZ).roots()
        []

    To compute its roots we need to factor the modulus `N` and use the Chinese
    remainder theorem::

        sage: p,q = N.prime_divisors()
        sage: f.change_ring(GF(p)).roots()
        [(4, 1)]
        sage: f.change_ring(GF(q)).roots()
        [(4, 1)]

        sage: crt(4, 4, p, q)
        4

    This root is quite small compared to `N`, so we can attempt to
    recover it without factoring `N` using Coppersmith's small root
    method::

        sage: f.small_roots()
        [4]

    An application of this method is to consider RSA. We are using 512-bit RSA
    with public exponent `e=3` to encrypt a 56-bit DES key. Because it would be
    easy to attack this setting if no padding was used we pad the key `K` with
    1s to get a large number::

        sage: Nbits, Kbits = 512, 56
        sage: e = 3

    We choose two primes of size 256-bit each::

        sage: p = 2^256 + 2^8 + 2^5 + 2^3 + 1
        sage: q = 2^256 + 2^8 + 2^5 + 2^3 + 2^2 + 1
        sage: N = p*q
        sage: ZmodN = Zmod( N )

    We choose a random key::

        sage: K = ZZ.random_element(0, 2^Kbits)

    and pad it with 512-56=456 1s::

        sage: Kdigits = K.digits(2)
        sage: M = [0]*Kbits + [1]*(Nbits-Kbits)
        sage: for i in range(len(Kdigits)): M[i] = Kdigits[i]

        sage: M = ZZ(M, 2)

    Now we encrypt the resulting message::

        sage: C = ZmodN(M)^e

    To recover `K` we consider the following polynomial modulo `N`::

        sage: P.<x> = PolynomialRing(ZmodN)
        sage: f = (2^Nbits - 2^Kbits + x)^e - C

    and recover its small roots::

        sage: Kbar = f.small_roots()[0]
        sage: K == Kbar
        True

    The same algorithm can be used to factor `N = pq` if partial
    knowledge about `q` is available. This example is from the
    Magma handbook:

    First, we set up `p`, `q` and `N`::

        sage: length = 512
        sage: hidden = 110
        sage: p = next_prime(2^int(round(length/2)))
        sage: q = next_prime( round(pi.n()*p) )
        sage: N = p*q

    Now we disturb the low 110 bits of `q`::

        sage: qbar = q + ZZ.random_element(0,2^hidden-1)

    And try to recover `q` from it::

        sage: F.<x> = PolynomialRing(Zmod(N))
        sage: f = x - qbar

    We know that the error is `\le 2^{\text{hidden}}-1` and that the modulus
    we are looking for is `\ge \sqrt{N}`::

        sage: set_verbose(2)
        sage: d = f.small_roots(X=2^hidden-1, beta=0.5)[0] # time random
        verbose 2 (<module>) m = 4
        verbose 2 (<module>) t = 4
        verbose 2 (<module>) X = 1298074214633706907132624082305023
        verbose 1 (<module>) LLL of 8x8 matrix (algorithm fpLLL:wrapper)
        verbose 1 (<module>) LLL finished (time = 0.006998)
        sage: q == qbar - d
        True

    REFERENCES:

    Don Coppersmith. *Finding a small root of a univariate modular equation.*
    In Advances in Cryptology, EuroCrypt 1996, volume 1070 of Lecture
    Notes in Computer Science, p. 155--165. Springer, 1996.
    http://cr.yp.to/bib/2001/coppersmith.pdf

    Alexander May. *New RSA Vulnerabilities Using Lattice Reduction Methods.*
    PhD thesis, University of Paderborn, 2003.
    http://www.cs.uni-paderborn.de/uploads/tx_sibibtex/bp.pdf
    """
    from sage.misc.misc import verbose
    from sage.matrix.constructor import Matrix
    from sage.rings.all import RR

    N = self.parent().characteristic()

    if not self.is_monic():
        raise ArithmeticError("Polynomial must be monic.")

    beta = RR(beta)
    if beta <= 0.0 or beta > 1.0:
        raise ValueError("0.0 < beta <= 1.0 not satisfied.")

    f = self.change_ring(ZZ)

    P,(x,) = f.parent().objgens()

    delta = f.degree()

    if epsilon is None:
        epsilon = beta/8
    verbose("epsilon = %d"%epsilon, level=2)

    m = max(beta**2/(delta * epsilon), 7*beta/delta).ceil()
    verbose("m = %d"%m, level=2)

    t = int( ( delta*m*(1/beta -1) ).floor() )
    verbose("t = %d"%t, level=2)

    if X is None:
        X = (0.5 * N**(beta**2/delta - epsilon)).ceil()
    verbose("X = %s"%X, level=2)

    # we could do this much faster, but this is a cheap step
    # compared to LLL
    g  = [x**j * N**(m-i) * f**i for i in range(m) for j in range(delta) ]
    g.extend([x**i * f**m for i in range(t)]) # h

    B = Matrix(ZZ, len(g), delta*m + max(delta,t) )
    for i in range(B.nrows()):
        for j in range( g[i].degree()+1 ):
            B[i,j] = g[i][j]*X**j

    B =  B.LLL(**kwds)

    f = sum([ZZ(B[0,i]//X**i)*x**i for i in range(B.ncols())])
    R = f.roots()

    ZmodN = self.base_ring()
    roots = set([ZmodN(r) for r,m in R if abs(r) <= X])
    Nbeta = N**beta
    return [root for root in roots if N.gcd(ZZ(self(root))) >= Nbeta]

cdef class Polynomial_dense_modn_ntl_zz(Polynomial_dense_mod_n):

    def __init__(self, parent, v=None, check=True, is_gen=False, construct=False):
        r"""
        EXAMPLES::

            sage: R = Integers(5**21) ; S.<x> = R[]
            sage: S(1/4)
            357627868652344
        """
        if isinstance(v, Polynomial):
            if (<Element>v)._parent == parent:
                Polynomial.__init__(self, parent, is_gen=is_gen)
                self.x = (<Polynomial_dense_modn_ntl_zz>v).x
                self.c = (<Polynomial_dense_modn_ntl_zz>v).c
                return

        Polynomial_dense_mod_n.__init__(self, parent, v, check=check, is_gen=is_gen, construct=construct)
        v = [a for a in self.__poly.list()]
        self.__poly = None # this will eventually go away
        cdef ntl_zz_pX ntl = ntl_zz_pX(v, parent.modulus()) # let it handle the hard work
        self.x = ntl.x
        self.c = ntl.c

    def __dealloc__(self):
        if <object>self.c is not None:
            self.c.restore_c()
        zz_pX_destruct(&self.x)

    def ntl_set_directly(self, v):
        # TODO: Get rid of this
        Polynomial_dense_mod_n.ntl_set_directly(self, v)
        # verbatim from __init__
        v = [int(a) for a in self.__poly.list()]
        self.__poly = None # this will eventually go away
        cdef ntl_zz_pX ntl = ntl_zz_pX(v, self._parent.modulus()) # let it handle the hard work
        self.x = ntl.x
        self.c = ntl.c

    cdef Polynomial_dense_modn_ntl_zz _new(self):
        cdef Polynomial_dense_modn_ntl_zz y = <Polynomial_dense_modn_ntl_zz>PY_NEW(Polynomial_dense_modn_ntl_zz)
        y.c = self.c
        y._parent = self._parent
        return y

    def int_list(self):
        """
        Returns the coefficients of self as efficiently as possible as a
        list of python ints.

        EXAMPLES::

            sage: R.<x> = Integers(100)[]
            sage: from sage.rings.polynomial.polynomial_modn_dense_ntl import Polynomial_dense_mod_n as poly_modn_dense
            sage: f = poly_modn_dense(R,[5,0,0,1])
            sage: f.int_list()
            [5, 0, 0, 1]
            sage: [type(a) for a in f.int_list()]
            [<type 'int'>, <type 'int'>, <type 'int'>, <type 'int'>]
        """
        cdef long i
        return [ zz_p_rep(zz_pX_GetCoeff(self.x, i)) for i from 0 <= i <= zz_pX_deg(self.x) ]

    def __getitem__(self, n):
        """
        Returns coefficient of the monomial of degree `n` if `n` is an integer,
        returns the monomials of self of degree in slice `n` if `n` is a slice.

        EXAMPLES::

            sage: R.<x> = Integers(100)[]
            sage: from sage.rings.polynomial.polynomial_modn_dense_ntl import Polynomial_dense_modn_ntl_zz
            sage: f = Polynomial_dense_modn_ntl_zz(R,[2, 1])^7
            sage: f[3]
            60
            sage: f[3:6]
            84*x^5 + 80*x^4 + 60*x^3
            sage: f[-5:50] == f
            True
            sage: f[6:]
            x^7 + 14*x^6
        """
        if isinstance(n, slice):
            start, stop = n.start, n.stop
            R = self.base_ring()
            if start < 0:
                start = 0
            if stop > zz_pX_deg(self.x)+1 or stop is None:
                stop = zz_pX_deg(self.x)+1
            v = [ zz_p_rep(zz_pX_GetCoeff(self.x, t)) for t from start <= t < stop ]
            return Polynomial_dense_modn_ntl_zz(self._parent, v, check=False) << start
        else:
            R = self._parent._base
            if n < 0 or n > zz_pX_deg(self.x):
                return R(0)
            else:
                return R(zz_p_rep(zz_pX_GetCoeff(self.x, n)))

    def _unsafe_mutate(self, n, value):
        self.c.restore_c()
        zz_pX_SetCoeff_long(self.x, n, value)

    cpdef ModuleElement _add_(self, ModuleElement _right):
        """
        TESTS::

            sage: R.<x> = Integers(100)[]
            sage: (x+5) + (x^2 - 6)
            x^2 + x + 99
        """
        cdef Polynomial_dense_modn_ntl_zz right = <Polynomial_dense_modn_ntl_zz>_right
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) + zz_pX_deg(right.x) > 1000000
        if do_sig: sig_on()
        self.c.restore_c()
        zz_pX_add(r.x, self.x, right.x)
        if do_sig: sig_off()
        return r

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        """
        TESTS::

            sage: R.<x> = Integers(100)[]
            sage: (x+5) - (x^2 - 6)
            99*x^2 + x + 11
        """
        cdef Polynomial_dense_modn_ntl_zz right = <Polynomial_dense_modn_ntl_zz>_right
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) + zz_pX_deg(right.x) > 1000000
        if do_sig: sig_on()
        self.c.restore_c()
        zz_pX_sub(r.x, self.x, right.x)
        if do_sig: sig_off()
        return r

    cpdef RingElement _mul_(self, RingElement _right):
        """
        TESTS::

            sage: R.<x> = Integers(100)[]
            sage: (x+5) * (x^2 - 1)
            x^3 + 5*x^2 + 99*x + 95
        """
        cdef Polynomial_dense_modn_ntl_zz right = <Polynomial_dense_modn_ntl_zz>_right
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) + zz_pX_deg(right.x) > 10000
        if do_sig: sig_on()
        self.c.restore_c()
        if self is right:
            zz_pX_sqr(r.x, self.x)
        else:
            zz_pX_mul(r.x, self.x, right.x)
        if do_sig: sig_off()
        return r

    cpdef Polynomial_dense_modn_ntl_zz _mul_trunc(self, Polynomial_dense_modn_ntl_zz right, long n):
        r"""
        Return the product of ``self`` and ``right`` truncated to the
        given length `n`, only return terms of degree less than `n`.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(100), implementation="NTL")
            sage: f = x - 2
            sage: g = x^2 - 8*x + 16
            sage: f*g
            x^3 + 90*x^2 + 32*x + 68
            sage: f._mul_trunc(g, 42)
            x^3 + 90*x^2 + 32*x + 68
            sage: f._mul_trunc(g, 3)
            90*x^2 + 32*x + 68
            sage: f._mul_trunc(g, 2)
            32*x + 68
            sage: f._mul_trunc(g, 1)
            68
            sage: f._mul_trunc(g, 0)
            0
            sage: f = x^2 - 8*x + 16
            sage: f._mul_trunc(f, 2)
            44*x + 56
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) + zz_pX_deg(right.x) > 10000
        if do_sig: sig_on()
        self.c.restore_c()
        if self is right:
            zz_pX_SqrTrunc(r.x, self.x, n)
        else:
            zz_pX_MulTrunc(r.x, self.x, right.x, n)
        if do_sig: sig_off()
        return r

    cpdef ModuleElement _rmul_(self, RingElement c):
        """
        TESTS::

            sage: R.<x> = Integers(100)[]
            sage: (x+5) * 3
            3*x + 15
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) > 100000
        if do_sig: sig_on()
        self.c.restore_c()
        zz_pX_rmul(r.x, self.x, c)
        if do_sig: sig_off()
        return r

    cpdef ModuleElement _lmul_(self, RingElement c):
        """
        TESTS::

            sage: R.<x> = Integers(100)[]
            sage: 3 * (x+5)
            3*x + 15
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef bint do_sig = zz_pX_deg(self.x) > 100000
        if do_sig: sig_on()
        self.c.restore_c()
        zz_pX_lmul(r.x, c, self.x)
        if do_sig: sig_off()
        return r

    def __pow__(Polynomial_dense_modn_ntl_zz self, ee, modulus):
        """
        TESTS::

            sage: R.<x> = Integers(100)[]
            sage: (x-1)^5
            x^5 + 95*x^4 + 10*x^3 + 90*x^2 + 5*x + 99

            sage: R.<x> = Integers(101)[]
            sage: (x-1)^(-5)
            1/(x^5 + 96*x^4 + 10*x^3 + 91*x^2 + 5*x + 100)

        We define ``0^0`` to be unity, :trac:`13895`::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: R(0)^0
            1

        The value returned from ``0^0`` should belong to our ring::

            sage: R.<x> = PolynomialRing(Integers(100), implementation='NTL')
            sage: type(R(0)^0) == type(R(0))
            True

        """
        cdef bint recip = 0, do_sig
        cdef long e = ee
        if e != ee:
            raise TypeError("Only integral powers defined.")
        elif e < 0:
            recip = 1
            e = -e
        if self == 0 and e == 0:
            return self.parent(1)

        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef zz_pX_Modulus_c *mod

        self.c.restore_c()

        if modulus is None:
            if zz_pX_IsX(self.x):
                zz_pX_LeftShift(r.x, self.x, e-1)
            else:
                do_sig = zz_pX_deg(self.x) *e > 1000
                if do_sig: sig_on()
                zz_pX_power(r.x, self.x, e)
                if do_sig: sig_off()
        else:
            if not PY_TYPE_CHECK(modulus, Polynomial_dense_modn_ntl_zz):
                modulus = self.parent()._coerce_(modulus)
            zz_pX_Modulus_construct(mod)
            zz_pX_Modulus_build(mod[0], (<Polynomial_dense_modn_ntl_zz>modulus).x)

            do_sig = zz_pX_deg(self.x) * e * self.c.p_bits > 1e5
            if do_sig: sig_on()
            zz_pX_PowerMod_long_pre(r.x, self.x, e, mod[0])
            if do_sig: sig_off()
            zz_pX_Modulus_destruct(mod)

        if recip:
            return ~r
        else:
            return r

    @coerce_binop
    def quo_rem(self, right):
        """
        Returns `q` and `r`, with the degree of `r` less than the degree of `right`,
        such that `q * right + r = self`.

        EXAMPLES::

            sage: R.<x> = Integers(125)[]
            sage: f = x^5+1; g = (x+1)^2
            sage: q, r = f.quo_rem(g)
            sage: q
            x^3 + 123*x^2 + 3*x + 121
            sage: r
            5*x + 5
            sage: q*g + r
            x^5 + 1
        """
        cdef Polynomial_dense_modn_ntl_zz q = self._new()
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        cdef Polynomial_dense_modn_ntl_zz denom = <Polynomial_dense_modn_ntl_zz>right
        sig_on()
        self.c.restore_c()
        zz_pX_divrem(q.x, r.x, self.x, denom.x)
        sig_off()
        return q, r

    def __floordiv__(self, right):
        """
        Returns the whole part of self/right, without remainder.

        For q = n // d, we have deg(n - q*d) < deg(d)

        EXAMPLES::

            sage: R.<x> = Integers(25)[]
            sage: f = x^7 + 1; g = x^2 - 1
            sage: q = f // g; q
            x^5 + x^3 + x
            sage: f - q*g
            x + 1
        """
        if PY_TYPE(self) != PY_TYPE(right) or (<Element>self)._parent is not (<Element>right)._parent:
            self, right = canonical_coercion(self, right)
            return self // right
        cdef Polynomial_dense_modn_ntl_zz numer = <Polynomial_dense_modn_ntl_zz>self
        cdef Polynomial_dense_modn_ntl_zz denom = <Polynomial_dense_modn_ntl_zz>right
        cdef Polynomial_dense_modn_ntl_zz q = numer._new()
        sig_on()
        numer.c.restore_c()
        zz_pX_div(q.x, numer.x, denom.x)
        sig_off()
        return q

    def __mod__(self, right):
        """
        EXAMPLES::

            sage: R.<x> = Integers(81)[]
            sage: f = x^7 + x + 1; g = x^3
            sage: r = f % g; r
            x + 1
            sage: g * x^4 + r
            x^7 + x + 1
        """
        if PY_TYPE(self) != PY_TYPE(right) or (<Element>self)._parent is not (<Element>right)._parent:
            self, right = canonical_coercion(self, right)
            return self % right
        cdef Polynomial_dense_modn_ntl_zz numer = <Polynomial_dense_modn_ntl_zz>self
        cdef Polynomial_dense_modn_ntl_zz denom = <Polynomial_dense_modn_ntl_zz>right
        cdef Polynomial_dense_modn_ntl_zz r = numer._new()
        sig_on()
        numer.c.restore_c()
        zz_pX_mod(r.x, numer.x, denom.x)
        sig_off()
        return r

    def shift(self, n):
        """
        Shift self to left by `n`, which is multiplication by `x^n`,
        truncating if `n` is negative.

        EXAMPLES::

            sage: R.<x> = Integers(77)[]
            sage: f = x^7 + x + 1
            sage: f.shift(1)
            x^8 + x^2 + x
            sage: f.shift(-1)
            x^6 + 1
            sage: f.shift(10).shift(-10) == f
            True

        TESTS::

            sage: p = R(0)
            sage: p.shift(3).is_zero()
            True
            sage: p.shift(-3).is_zero()
            True

        """
        return self << n

    def __lshift__(Polynomial_dense_modn_ntl_zz self, long n):
        """
        TEST::

            sage: R.<x> = Integers(77)[]
            sage: f = x^5 + 2*x + 1
            sage: f << 3
            x^8 + 2*x^4 + x^3
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        zz_pX_LeftShift(r.x, self.x, n)
        return r

    def __rshift__(Polynomial_dense_modn_ntl_zz self, long n):
        """
        TEST::

            sage: R.<x> = Integers(77)[]
            sage: f = x^5 + 2*x + 1
            sage: f >> 3
            x^2
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        zz_pX_RightShift(r.x, self.x, n)
        return r

    def _derivative(self, var=None):
        r"""
        Returns the formal derivative of self with respect to var.

        var must be either the generator of the polynomial ring to which
        this polynomial belongs, or None (either way the behaviour is the
        same).

        .. seealso:: :meth:`.derivative`

        EXAMPLES::

            sage: R.<x> = Integers(77)[]
            sage: f = x^4 - x - 1
            sage: f._derivative()
            4*x^3 + 76
            sage: f._derivative(None)
            4*x^3 + 76

            sage: f._derivative(2*x)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to 2*x

            sage: y = var("y")
            sage: f._derivative(y)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to y
        """
        if var is not None and var is not self._parent.gen():
            raise ValueError("cannot differentiate with respect to %s" % var)

        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        zz_pX_diff(r.x, self.x)
        return r

    def reverse(self):
        """
        Reverses the coefficients of self. The reverse of `f(x)` is `x^n f(1/x)`.

        The degree will go down if the constant term is zero.

        EXAMPLES::

            sage: R.<x> = Integers(77)[]
            sage: f = x^4 - x - 1
            sage: f.reverse()
            76*x^4 + 76*x^3 + 1
            sage: f = x^3 - x
            sage: f.reverse()
            76*x^2 + 1
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        zz_pX_reverse(r.x, self.x)
        return r

    def is_gen(self):
        return zz_pX_IsX(self.x)

    def __nonzero__(self):
        """
        TESTS::

            sage: R.<x> = Integers(77)[]
            sage: f = x^4 - x - 1
            sage: not f
            False
            sage: not (x-x)
            True
        """
        return not zz_pX_IsZero(self.x)

    def valuation(self):
        """
        Returns the valuation of self, that is, the power of the
        lowest non-zero monomial of self.

        EXAMPLES::

            sage: R.<x> = Integers(10)[]
            sage: x.valuation()
            1
            sage: f = x-3; f.valuation()
            0
            sage: f = x^99; f.valuation()
            99
            sage: f = x-x; f.valuation()
            +Infinity
        """
        cdef long n
        for n from 0 <= n <= zz_pX_deg(self.x):
            if zz_p_rep(zz_pX_GetCoeff(self.x, n)):
                return n
        return infinity

    def degree(self):
        """
        EXAMPLES::

            sage: R.<x> = Integers(77)[]
            sage: f = x^4 - x - 1
            sage: f.degree()
            4
            sage: f = 77*x + 1
            sage: f.degree()
            0
        """
        return zz_pX_deg(self.x)

    cpdef Polynomial truncate(self, long n):
        """
        Returns this polynomial mod `x^n`.

        EXAMPLES::

            sage: R.<x> = Integers(77)[]
            sage: f = sum(x^n for n in range(10)); f
            x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
            sage: f.truncate(6)
            x^5 + x^4 + x^3 + x^2 + x + 1
        """
        cdef Polynomial_dense_modn_ntl_zz r = self._new()
        zz_pX_trunc(r.x, self.x, n)
        return r

    def __call__(self, *args, **kwds):
        """
        Evaluate self at x. If x is a single argument coercible into
        the base ring of self, this is done directly in NTL, otherwise
        the generic Polynomial call code is used.

        EXAMPLES::

            sage: R.<x> = Integers(100)[]
            sage: f = x^3+7
            sage: f(5)
            32
            sage: f(5r)
            32
            sage: f(mod(5, 1000))
            32
            sage: f(x)
            x^3 + 7
            sage: S.<y> = Integers(5)[]
            sage: f(y)
            y^3 + 2
        """
        if len(args) != 1 or len(kwds) != 0:
            return Polynomial.__call__(self, *args, **kwds)
        arg = args[0]
        cdef ntl_zz_p fx = ntl_zz_p(0, self.c), x = None
        if PY_TYPE_CHECK(arg, int):
            x = ntl_zz_p(arg, self.c)
        elif PY_TYPE_CHECK(arg, Integer):
            x = ntl_zz_p(arg, self.c)
        elif PY_TYPE_CHECK(arg, Element):
            if <void *>self._parent._base == <void *>(<Element>arg)._parent: # c++ pointer hack
                x = ntl_zz_p(arg, self.c)
            else:
                map = self._parent._base.coerce_map_from((<Element>arg)._parent)
                if map is not None:
                    x = ntl_zz_p(map(arg), self.c)
        if <PyObject *>x == <PyObject *>None: # c++ pointer compare error
            return Polynomial.__call__(self, *args, **kwds)
        else:
            zz_pX_eval(fx.x, self.x, x.x)
            return self._parent(int(fx))



cdef class Polynomial_dense_modn_ntl_ZZ(Polynomial_dense_mod_n):

    def __init__(self, parent, v=None, check=True, is_gen=False, construct=False):
        if isinstance(v, Polynomial):
            if (<Element>v)._parent == parent:
                Polynomial.__init__(self, parent, is_gen=is_gen)
                self.x = (<Polynomial_dense_modn_ntl_ZZ>v).x
                self.c = (<Polynomial_dense_modn_ntl_ZZ>v).c
                return

        Polynomial_dense_mod_n.__init__(self, parent, v, check=check, is_gen=is_gen, construct=construct)
        cdef ntl_ZZ_pX ntl = self.__poly
        self.__poly = None # this will eventually go away
        self.x = ntl.x
        self.c = ntl.c

    def __dealloc__(self):
        if <object>self.c is not None:
            self.c.restore_c()
        ZZ_pX_destruct(&self.x)


    cdef Polynomial_dense_modn_ntl_ZZ _new(self):
        cdef Polynomial_dense_modn_ntl_ZZ y = <Polynomial_dense_modn_ntl_ZZ>PY_NEW(Polynomial_dense_modn_ntl_ZZ)
        y.c = self.c
        y._parent = self._parent
        return y


    def list(self):
        return [self._parent._base(self[n]) for n from 0 <= n <= self.degree()]

    def __getitem__(self, n):
        """
        Returns coefficient of the monomial of degree `n` if `n` is an integer,
        returns the monomials of self of degree in slice `n` if `n` is a slice.

        EXAMPLES::

            sage: R.<x> = Integers(10^30)[]
            sage: from sage.rings.polynomial.polynomial_modn_dense_ntl import Polynomial_dense_modn_ntl_ZZ
            sage: f = Polynomial_dense_modn_ntl_ZZ(R,[2,1])^7
            sage: f[3]
            560
            sage: f[3:6]
            84*x^5 + 280*x^4 + 560*x^3
            sage: f[-5:50] == f
            True
            sage: f[6:]
            x^7 + 14*x^6
        """
        if isinstance(n, slice):
            start, stop = n.start, n.stop
            R = self.base_ring()
            if start < 0:
                start = 0
            if stop > ZZ_pX_deg(self.x)+1 or stop is None:
                stop = ZZ_pX_deg(self.x)+1
            v = [ self[t] for t from start <= t < stop ]
            return Polynomial_dense_modn_ntl_ZZ(self._parent, v, check=False) << start
        else:
            R = self._parent._base
            if n < 0 or n > ZZ_pX_deg(self.x):
                return R(0)

        self.c.restore_c()
        cdef Integer z

        # TODO, make this faster
        cdef ntl_ZZ_p ntl = ntl_ZZ_p(0, self.c)
        ntl.x = ZZ_pX_coeff(self.x, n)
        return R(ntl._integer_())

    def _unsafe_mutate(self, n, value):
        self.c.restore_c()
        cdef Integer a
        if PY_TYPE_CHECK(value, Integer):
            a = <Integer>value
        else:
            a = ZZ(value)
        cdef ntl_ZZ_p val = ntl_ZZ_p(a, self.c)
        ZZ_pX_SetCoeff(self.x, n, val.x)

    cpdef ModuleElement _add_(self, ModuleElement _right):
        """
        TESTS::

            sage: R.<x> = Integers(10^30)[]
            sage: (x+5) + (x^2 - 6)
            x^2 + x + 999999999999999999999999999999
        """
        cdef Polynomial_dense_modn_ntl_ZZ right = <Polynomial_dense_modn_ntl_ZZ>_right
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef bint do_sig = (ZZ_pX_deg(self.x) + ZZ_pX_deg(right.x)) * self.c.p_bits > 1e7
        if do_sig: sig_on()
        self.c.restore_c()
        ZZ_pX_add(r.x, self.x, right.x)
        if do_sig: sig_off()
        return r

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        """
        TESTS::

            sage: R.<x> = Integers(10^30)[]
            sage: (x+5) - (x^2 - 6)
            999999999999999999999999999999*x^2 + x + 11
        """
        cdef Polynomial_dense_modn_ntl_ZZ right = <Polynomial_dense_modn_ntl_ZZ>_right
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef bint do_sig = (ZZ_pX_deg(self.x) + ZZ_pX_deg(right.x)) * self.c.p_bits > 1e7
        if do_sig: sig_on()
        self.c.restore_c()
        ZZ_pX_sub(r.x, self.x, right.x)
        if do_sig: sig_off()
        return r

    cpdef RingElement _mul_(self, RingElement _right):
        """
        TESTS::

            sage: R.<x> = Integers(10^30)[]
            sage: (x+5) * (x^2 - 1)
            x^3 + 5*x^2 + 999999999999999999999999999999*x + 999999999999999999999999999995
        """
        cdef Polynomial_dense_modn_ntl_ZZ right = <Polynomial_dense_modn_ntl_ZZ>_right
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef bint do_sig = (ZZ_pX_deg(self.x) + ZZ_pX_deg(right.x)) * self.c.p_bits > 1e5
        if do_sig: sig_on()
        self.c.restore_c()
        if self is right:
            ZZ_pX_sqr(r.x, self.x)
        else:
            ZZ_pX_mul(r.x, self.x, right.x)
        if do_sig: sig_off()
        return r

    cpdef Polynomial_dense_modn_ntl_ZZ _mul_trunc(self, Polynomial_dense_modn_ntl_ZZ right, long n):
        """
        Return the product of ``self`` and ``right`` truncated to the
        given length `n`, only return terms of degree less than `n`.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation="NTL")
            sage: f = x - 2
            sage: g = x^2 - 8*x + 16
            sage: f*g
            x^3 + 999999999999999999999999999990*x^2 + 32*x + 999999999999999999999999999968
            sage: f._mul_trunc(g, 42)
            x^3 + 999999999999999999999999999990*x^2 + 32*x + 999999999999999999999999999968
            sage: f._mul_trunc(g, 3)
            999999999999999999999999999990*x^2 + 32*x + 999999999999999999999999999968
            sage: f._mul_trunc(g, 2)
            32*x + 999999999999999999999999999968
            sage: f._mul_trunc(g, 1)
            999999999999999999999999999968
            sage: f._mul_trunc(g, 0)
            0
            sage: f = x^2 - 8*x + 16
            sage: f._mul_trunc(f, 2)
            999999999999999999999999999744*x + 256
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef bint do_sig = (ZZ_pX_deg(self.x) + ZZ_pX_deg(right.x)) * self.c.p_bits > 1e5
        if do_sig: sig_on()
        self.c.restore_c()
        if self is right:
            ZZ_pX_SqrTrunc(r.x, self.x, n)
        else:
            ZZ_pX_MulTrunc(r.x, self.x, right.x, n)
        if do_sig: sig_off()
        return r

    cpdef ModuleElement _rmul_(self, RingElement c):
        """
        TESTS::

            sage: R.<x> = Integers(10^30)[]
            sage: (x+5) * 3
            3*x + 15
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef bint do_sig = ZZ_pX_deg(self.x) * self.c.p_bits > 1e7
        if do_sig: sig_on()
        self.c.restore_c()
        cdef ntl_ZZ_p value = ntl_ZZ_p(c, self.c)
        ZZ_pX_rmul(r.x, self.x, value.x)
        if do_sig: sig_off()
        return r

    cpdef ModuleElement _lmul_(self, RingElement c):
        """
        TESTS::

            sage: R.<x> = Integers(10^30)[]
            sage: 3 * (x+5)
            3*x + 15
        """
        return self._rmul_(c)

    def __pow__(Polynomial_dense_modn_ntl_ZZ self, ee, modulus):
        """
        TESTS::

            sage: R.<x> = Integers(10^30)[]
            sage: (x+1)^5
            x^5 + 5*x^4 + 10*x^3 + 10*x^2 + 5*x + 1

        We define ``0^0`` to be unity, :trac:`13895`::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: R(0)^0
            1

        The value returned from ``0^0`` should belong to our ring::

            sage: R.<x> = PolynomialRing(Integers(10^30), implementation='NTL')
            sage: type(R(0)^0) == type(R(0))
            True

        """
        cdef bint recip = 0, do_sig
        cdef long e = ee
        if e != ee:
            raise TypeError("Only integral powers defined.")
        elif e < 0:
            recip = 1 # delay because powering frac field elements is slow
            e = -e
        if self == 0 and e == 0:
            return self.parent(1)
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef ZZ_pX_Modulus_c *mod

        self.c.restore_c()

        if modulus is None:
            if ZZ_pX_IsX(self.x):
                ZZ_pX_LeftShift(r.x, self.x, e - 1)
            else:
                do_sig = ZZ_pX_deg(self.x) * e * self.c.p_bits > 1e5
                if do_sig: sig_on()
                ZZ_pX_power(r.x, self.x, e)
                if do_sig: sig_off()
        else:
            if not PY_TYPE_CHECK(modulus, Polynomial_dense_modn_ntl_ZZ):
                modulus = self.parent()._coerce_(modulus)
            ZZ_pX_Modulus_construct(mod)
            ZZ_pX_Modulus_build(mod[0], (<Polynomial_dense_modn_ntl_ZZ>modulus).x)

            do_sig = ZZ_pX_deg(self.x) * e * self.c.p_bits > 1e5
            if do_sig: sig_on()
            ZZ_pX_PowerMod_long_pre(r.x, self.x, e, mod[0])
            if do_sig: sig_off()
            ZZ_pX_Modulus_destruct(mod)
        if recip:
            return ~r
        else:
            return r

    @coerce_binop
    def quo_rem(self, right):
        """
        Returns `q` and `r`, with the degree of `r` less than the degree of `right`,
        such that `q * right + r = self`.

        EXAMPLES::

            sage: R.<x> = Integers(10^30)[]
            sage: f = x^5+1; g = (x+1)^2
            sage: q, r = f.quo_rem(g)
            sage: q
            x^3 + 999999999999999999999999999998*x^2 + 3*x + 999999999999999999999999999996
            sage: r
            5*x + 5
            sage: q*g + r
            x^5 + 1
        """
        cdef Polynomial_dense_modn_ntl_ZZ q = self._new()
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        cdef Polynomial_dense_modn_ntl_ZZ denom = <Polynomial_dense_modn_ntl_ZZ>right
        sig_on()
        self.c.restore_c()
        ZZ_pX_DivRem(q.x, r.x, self.x, denom.x)
        sig_off()
        return q, r

    def __floordiv__(self, right):
        """
        Returns the whole part of self/right, without remainder.

        For q = n // d, we have deg(n - q*d) < deg(d)

        EXAMPLES::

            sage: R.<x> = Integers(10^30)[]
            sage: f = x^7 + 1; g = x^2 - 1
            sage: q = f // g; q
            x^5 + x^3 + x
            sage: f - q*g
            x + 1
        """
        if PY_TYPE(self) != PY_TYPE(right) or (<Element>self)._parent is not (<Element>right)._parent:
            self, right = canonical_coercion(self, right)
            return self // right
        cdef Polynomial_dense_modn_ntl_ZZ numer = <Polynomial_dense_modn_ntl_ZZ>self
        cdef Polynomial_dense_modn_ntl_ZZ denom = <Polynomial_dense_modn_ntl_ZZ>right
        cdef Polynomial_dense_modn_ntl_ZZ q = numer._new()
        sig_on()
        numer.c.restore_c()
        ZZ_pX_div(q.x, numer.x, denom.x)
        sig_off()
        return q

    def __mod__(self, right):
        """
        EXAMPLES::

            sage: R.<x> = Integers(9^30)[]
            sage: f = x^7 + x + 1; g = x^3 - 1
            sage: r = f % g; r
            2*x + 1
            sage: g * (x^4 + x) + r
            x^7 + x + 1
        """
        if PY_TYPE(self) != PY_TYPE(right) or (<Element>self)._parent is not (<Element>right)._parent:
            self, right = canonical_coercion(self, right)
            return self % right
        cdef Polynomial_dense_modn_ntl_ZZ numer = <Polynomial_dense_modn_ntl_ZZ>self
        cdef Polynomial_dense_modn_ntl_ZZ denom = <Polynomial_dense_modn_ntl_ZZ>right
        cdef Polynomial_dense_modn_ntl_ZZ r = numer._new()
        sig_on()
        numer.c.restore_c()
        ZZ_pX_rem(r.x, numer.x, denom.x)
        sig_off()
        return r

    def shift(self, n):
        """
        Shift self to left by `n`, which is multiplication by `x^n`,
        truncating if `n` is negative.

        EXAMPLES::

            sage: R.<x> = Integers(12^30)[]
            sage: f = x^7 + x + 1
            sage: f.shift(1)
            x^8 + x^2 + x
            sage: f.shift(-1)
            x^6 + 1
            sage: f.shift(10).shift(-10) == f
            True

        TESTS::

            sage: p = R(0)
            sage: p.shift(3).is_zero()
            True
            sage: p.shift(-3).is_zero()
            True

        """
        return self << n

    def __lshift__(Polynomial_dense_modn_ntl_ZZ self, long n):
        """
        TEST::

            sage: R.<x> = Integers(14^30)[]
            sage: f = x^5 + 2*x + 1
            sage: f << 3
            x^8 + 2*x^4 + x^3
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        ZZ_pX_LeftShift(r.x, self.x, n)
        return r

    def __rshift__(Polynomial_dense_modn_ntl_ZZ self, long n):
        """
        TEST::

            sage: R.<x> = Integers(15^30)[]
            sage: f = x^5 + 2*x + 1
            sage: f >> 3
            x^2
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        ZZ_pX_RightShift(r.x, self.x, n)
        return r


    def _derivative(self, var=None):
        r"""
        Returns the formal derivative of self with respect to var.

        var must be either the generator of the polynomial ring to which
        this polynomial belongs, or None (either way the behaviour is the
        same).

        .. seealso:: :meth:`.derivative`

        EXAMPLES::

            sage: R.<x> = Integers(12^29)[]
            sage: f = x^4 + x + 5
            sage: f._derivative()
            4*x^3 + 1
            sage: f._derivative(None)
            4*x^3 + 1

            sage: f._derivative(2*x)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to 2*x

            sage: y = var("y")
            sage: f._derivative(y)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to y
        """
        if var is not None and var is not self._parent.gen():
            raise ValueError("cannot differentiate with respect to %s" % var)

        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        ZZ_pX_diff(r.x, self.x)
        return r


    def reverse(self):
        """
        Reverses the coefficients of self. The reverse of `f(x)` is `x^n
        f(1/x)`.

        The degree will go down if the constant term is zero.

        EXAMPLES::

            sage: R.<x> = Integers(12^29)[]
            sage: f = x^4 + 2*x + 5
            sage: f.reverse()
            5*x^4 + 2*x^3 + 1
            sage: f = x^3 + x
            sage: f.reverse()
            x^2 + 1
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        ZZ_pX_reverse(r.x, self.x)
        return r

    def is_gen(self):
        return ZZ_pX_IsX(self.x)

    def valuation(self):
        """
        Returns the valuation of self, that is, the power of the
        lowest non-zero monomial of self.

        EXAMPLES::

            sage: R.<x> = Integers(10^50)[]
            sage: x.valuation()
            1
            sage: f = x-3; f.valuation()
            0
            sage: f = x^99; f.valuation()
            99
            sage: f = x-x; f.valuation()
            +Infinity
        """
        cdef long n
        cdef ZZ_p_c coeff
        for n from 0 <= n <= ZZ_pX_deg(self.x):
            coeff = ZZ_pX_coeff(self.x, n)
            if not ZZ_p_IsZero(coeff):
                return n
        return infinity

    def __nonzero__(self):
        """
        TESTS::

            sage: R.<x> = Integers(12^29)[]
            sage: f = x^4 + 1
            sage: not f
            False
            sage: not (x-x)
            True
        """
        return not ZZ_pX_IsZero(self.x)

    def degree(self):
        """
        EXAMPLES::

            sage: R.<x> = Integers(14^34)[]
            sage: f = x^4 - x - 1
            sage: f.degree()
            4
            sage: f = 14^43*x + 1
            sage: f.degree()
            0
        """
        return ZZ_pX_deg(self.x)

    cpdef Polynomial truncate(self, long n):
        """
        Returns this polynomial mod `x^n`.

        EXAMPLES::

            sage: R.<x> = Integers(15^30)[]
            sage: f = sum(x^n for n in range(10)); f
            x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
            sage: f.truncate(6)
            x^5 + x^4 + x^3 + x^2 + x + 1
        """
        cdef Polynomial_dense_modn_ntl_ZZ r = self._new()
        ZZ_pX_trunc(r.x, self.x, n)
        return r

    def __call__(self, *args, **kwds):
        """
        Evaluate self at x. If x is a single argument coercible into
        the base ring of self, this is done directly in NTL, otherwise
        the generic Polynomial call code is used.

        EXAMPLES::

            sage: R.<x> = Integers(10^30)[]
            sage: f = x^3+7
            sage: f(5)
            132
            sage: f(5r)
            132
            sage: f(mod(5, 10^50))
            132
            sage: f(x)
            x^3 + 7
            sage: S.<y> = Integers(5)[]
            sage: f(y)
            y^3 + 2
        """
        if len(args) != 1 or len(kwds) != 0:
            return Polynomial.__call__(self, *args, **kwds)
        arg = args[0]
        cdef ntl_ZZ_p fx = ntl_ZZ_p(0, self.c), x = None
        if PY_TYPE_CHECK(arg, int) or PY_TYPE_CHECK(arg, Integer):
            x = ntl_ZZ_p(arg, self.c)
        elif PY_TYPE_CHECK(arg, Element):
            if <void *>self._parent._base == <void *>(<Element>arg)._parent: # c++ pointer hack
                x = ntl_ZZ_p(arg, self.c)
            else:
                map = self._parent._base.coerce_map_from((<Element>arg)._parent)
                if map is not None:
                    x = ntl_ZZ_p(map(arg), self.c)
        if <PyObject *>x == <PyObject *>None: # c++ pointer compare error
            return Polynomial.__call__(self, *args, **kwds)
        else:
            ZZ_pX_eval(fx.x, self.x, x.x)
            return self._parent(fx._integer_())


cdef class Polynomial_dense_mod_p(Polynomial_dense_mod_n):
    """
    A dense polynomial over the integers modulo p, where p is prime.
    """

    @coerce_binop
    def gcd(self, right):
        """
        Return the greatest common divisor of this polynomial and ``other``, as
        a monic polynomial.

        INPUT:

        - ``other`` -- a polynomial defined over the same ring as ``self``

        EXAMPLES::

            sage: R.<x> = PolynomialRing(GF(3),implementation="NTL")
            sage: f,g = x + 2, x^2 - 1
            sage: f.gcd(g)
            x + 2

        """
        g = self.ntl_ZZ_pX().gcd(right.ntl_ZZ_pX())
        return self.parent()(g, construct=True)

    @coerce_binop
    def xgcd(self, other):
        r"""
        Compute the extended gcd of this element and ``other``.

        INPUT:

        - ``other`` -- an element in the same polynomial ring

        OUTPUT:

        A tuple ``(r,s,t)`` of elements in the polynomial ring such
        that ``r = s*self + t*other``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(GF(3),implementation='NTL')
            sage: x.xgcd(x)
            (x, 0, 1)
            sage: (x^2 - 1).xgcd(x - 1)
            (x + 2, 0, 1)
            sage: R.zero().xgcd(R.one())
            (1, 0, 1)
            sage: (x^3 - 1).xgcd((x - 1)^2)
            (x^2 + x + 1, 0, 1)
            sage: ((x - 1)*(x + 1)).xgcd(x*(x - 1))
            (x + 2, 1, 2)

        """
        r, s, t = self.ntl_ZZ_pX().xgcd(other.ntl_ZZ_pX())
        return self.parent()(r, construct=True), self.parent()(s, construct=True), self.parent()(t, construct=True)

    @coerce_binop
    def resultant(self, other):
        """
        Returns the resultant of self and other, which must lie in the same
        polynomial ring.

        INPUT:

        - ``other`` -- a polynomial

        OUTPUT: an element of the base ring of the polynomial ring

        EXAMPLES::

            sage: R.<x> = GF(19)['x']
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: r = f.resultant(g); r
            11
            sage: r.parent() is GF(19)
            True
        """
        other = self.parent()._coerce_(other)
        return self.base_ring()(str(self.ntl_ZZ_pX().resultant(other.ntl_ZZ_pX())))

    def discriminant(self):
        """
        EXAMPLES::

            sage: _.<x> = PolynomialRing(GF(19))
            sage: f = x^3 + 3*x - 17
            sage: f.discriminant()
            12
        """
        return self.base_ring()(str(self.ntl_ZZ_pX().discriminant()))
