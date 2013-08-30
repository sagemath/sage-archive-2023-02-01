"""
PARI C-library interface

AUTHORS:

- William Stein (2006-03-01): updated to work with PARI 2.2.12-beta

- William Stein (2006-03-06): added newtonpoly

- Justin Walker: contributed some of the function definitions

- Gonzalo Tornaria: improvements to conversions; much better error
  handling.

- Robert Bradshaw, Jeroen Demeyer, William Stein (2010-08-15):
  Upgrade to PARI 2.4.3 (#9343)

- Jeroen Demeyer (2011-11-12): rewrite various conversion routines
  (#11611, #11854, #11952)


EXAMPLES::

    sage: pari('5! + 10/x')
    (120*x + 10)/x
    sage: pari('intnum(x=0,13,sin(x)+sin(x^2) + x)')
    85.1885681951527
    sage: f = pari('x^3-1')
    sage: v = f.factor(); v
    [x - 1, 1; x^2 + x + 1, 1]
    sage: v[0]   # indexing is 0-based unlike in GP.
    [x - 1, x^2 + x + 1]~
    sage: v[1]
    [1, 1]~

Arithmetic obeys the usual coercion rules::

    sage: type(pari(1) + 1)
    <type 'sage.libs.pari.gen.gen'>
    sage: type(1 + pari(1))
    <type 'sage.libs.pari.gen.gen'>

GUIDE TO REAL PRECISION AND THE PARI LIBRARY

The default real precision in communicating with the Pari library
is the same as the default Sage real precision, which is 53 bits.
Inexact Pari objects are therefore printed by default to 15 decimal
digits (even if they are actually more precise).

Default precision example (53 bits, 15 significant decimals)::

    sage: a = pari(1.23); a
    1.23000000000000
    sage: a.sin()
    0.942488801931698

Example with custom precision of 200 bits (60 significant
decimals)::

    sage: R = RealField(200)
    sage: a = pari(R(1.23)); a   # only 15 significant digits printed
    1.23000000000000
    sage: R(a)         # but the number is known to precision of 200 bits
    1.2300000000000000000000000000000000000000000000000000000000
    sage: a.sin()      # only 15 significant digits printed
    0.942488801931698
    sage: R(a.sin())   # but the number is known to precision of 200 bits
    0.94248880193169751002382356538924454146128740562765030213504

It is possible to change the number of printed decimals::

    sage: R = RealField(200)    # 200 bits of precision in computations
    sage: old_prec = pari.set_real_precision(60)  # 60 decimals printed
    sage: a = pari(R(1.23)); a
    1.23000000000000000000000000000000000000000000000000000000000
    sage: a.sin()
    0.942488801931697510023823565389244541461287405627650302135038
    sage: pari.set_real_precision(old_prec)  # restore the default printing behavior
    60

Unless otherwise indicated in the docstring, most Pari functions
that return inexact objects use the precision of their arguments to
decide the precision of the computation. However, if some of these
arguments happen to be exact numbers (integers, rationals, etc.),
an optional parameter indicates the precision (in bits) to which
these arguments should be converted before the computation. If this
precision parameter is missing, the default precision of 53 bits is
used. The following first converts 2 into a real with 53-bit
precision::

    sage: R = RealField()
    sage: R(pari(2).sin())
    0.909297426825682

We can ask for a better precision using the optional parameter::

    sage: R = RealField(150)
    sage: R(pari(2).sin(precision=150))
    0.90929742682568169539601986591174484270225497

Warning regarding conversions Sage - Pari - Sage: Some care must be
taken when juggling inexact types back and forth between Sage and
Pari. In theory, calling p=pari(s) creates a Pari object p with the
same precision as s; in practice, the Pari library's precision is
word-based, so it will go up to the next word. For example, a
default 53-bit Sage real s will be bumped up to 64 bits by adding
bogus 11 bits. The function p.python() returns a Sage object with
exactly the same precision as the Pari object p. So
pari(s).python() is definitely not equal to s, since it has 64 bits
of precision, including the bogus 11 bits. The correct way of
avoiding this is to coerce pari(s).python() back into a domain with
the right precision. This has to be done by the user (or by Sage
functions that use Pari library functions in gen.pyx). For
instance, if we want to use the Pari library to compute sqrt(pi)
with a precision of 100 bits::

    sage: R = RealField(100)
    sage: s = R(pi); s
    3.1415926535897932384626433833
    sage: p = pari(s).sqrt()
    sage: x = p.python(); x  # wow, more digits than I expected!
    1.7724538509055160272981674833410973484
    sage: x.prec()           # has precision 'improved' from 100 to 128?
    128
    sage: x == RealField(128)(pi).sqrt()  # sadly, no!
    False
    sage: R(x)               # x should be brought back to precision 100
    1.7724538509055160272981674833
    sage: R(x) == s.sqrt()
    True

Elliptic curves and precision: If you are working with elliptic
curves and want to compute with a precision other than the default
53 bits, you should use the precision parameter of ellinit()::

    sage: R = RealField(150)
    sage: e = pari([0,0,0,-82,0]).ellinit(precision=150)
    sage: eta1 = e.elleta()[0]
    sage: R(eta1)
    3.6054636014326520859158205642077267748102690

Number fields and precision: TODO

TESTS:

Check that output from PARI's print command is actually seen by
Sage (ticket #9636)::

    sage: pari('print("test")')
    test

"""

#*****************************************************************************
#       Copyright (C) 2006,2010 William Stein <wstein@gmail.com>
#       Copyright (C) ???? Justin Walker
#       Copyright (C) ???? Gonzalo Tornaria
#       Copyright (C) 2010 Robert Bradshaw <robertwb@math.washington.edu>
#       Copyright (C) 2010,2011 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import sys
import math
import types
import operator
import sage.structure.element
from sage.structure.element cimport ModuleElement, RingElement, Element
from sage.structure.parent cimport Parent
from sage.misc.randstate cimport randstate, current_randstate

from sage.misc.misc_c import is_64_bit

include 'pari_err.pxi'
include 'sage/ext/stdsage.pxi'
include 'sage/ext/python.pxi'

cdef extern from "mpz_pylong.h":
    cdef int mpz_set_pylong(mpz_t dst, src) except -1

# Make sure we don't use mpz_t_offset before initializing it by
# putting in a value that's likely to provoke a segmentation fault,
# rather than silently corrupting memory.
cdef long mpz_t_offset = 1000000000

# so Galois groups are represented in a sane way
# See the polgalois section of the PARI users manual.
new_galois_format = 1

cdef pari_sp mytop

# real precision in decimal digits: see documentation for
# get_real_precision() and set_real_precision().  This variable is used
# in gp to set the precision of input quantities (e.g. sqrt(2)), and for
# determining the number of digits to be printed.  It is *not* used as
# a "default precision" for internal computations, which always use
# the actual precision of arguments together (where relevant) with a
# "prec" parameter.  In ALL cases (for real computations) the prec
# parameter is a WORD precision and NOT decimal precision.  Pari reals
# with word precision w have bit precision (of the mantissa) equal to
# 32*(w-2) or 64*(w-2).
#
# Hence the only relevance of this parameter in Sage is (1) for the
# output format of components of objects of type
# 'sage.libs.pari.gen.gen'; (2) for setting the precision of pari
# variables created from strings (e.g. via sage: pari('1.2')).
#
# WARNING: Many pari library functions take a last parameter "prec"
# which should be a words precision.  In many cases this is redundant
# and is simply ignored.  In our wrapping of these functions we use
# the variable prec here for convenience only.
cdef unsigned long prec

#################################################################
# conversions between various real precision models
#################################################################

def prec_bits_to_dec(int prec_in_bits):
    r"""
    Convert from precision expressed in bits to precision expressed in
    decimal.

    EXAMPLES::

        sage: import sage.libs.pari.gen as gen
        sage: gen.prec_bits_to_dec(53)
        15
        sage: [(32*n,gen.prec_bits_to_dec(32*n)) for n in range(1,9)]
        [(32, 9),
        (64, 19),
        (96, 28),
        (128, 38),
        (160, 48),
        (192, 57),
        (224, 67),
        (256, 77)]
    """
    log_2 = 0.301029995663981
    return int(prec_in_bits*log_2)

def prec_dec_to_bits(int prec_in_dec):
    r"""
    Convert from precision expressed in decimal to precision expressed
    in bits.

    EXAMPLES::

        sage: import sage.libs.pari.gen as gen
        sage: gen.prec_dec_to_bits(15)
        49
        sage: [(n,gen.prec_dec_to_bits(n)) for n in range(10,100,10)]
        [(10, 33),
        (20, 66),
        (30, 99),
        (40, 132),
        (50, 166),
        (60, 199),
        (70, 232),
        (80, 265),
        (90, 298)]
    """
    log_10 = 3.32192809488736
    return int(prec_in_dec*log_10)

def prec_bits_to_words(int prec_in_bits=0):
    r"""
    Convert from precision expressed in bits to pari real precision
    expressed in words. Note: this rounds up to the nearest word,
    adjusts for the two codewords of a pari real, and is
    architecture-dependent.

    EXAMPLES::

        sage: import sage.libs.pari.gen as gen
        sage: gen.prec_bits_to_words(70)
        5   # 32-bit
        4   # 64-bit

    ::

        sage: [(32*n,gen.prec_bits_to_words(32*n)) for n in range(1,9)]
        [(32, 3), (64, 4), (96, 5), (128, 6), (160, 7), (192, 8), (224, 9), (256, 10)] # 32-bit
        [(32, 3), (64, 3), (96, 4), (128, 4), (160, 5), (192, 5), (224, 6), (256, 6)] # 64-bit
    """
    if not prec_in_bits:
        return prec
    cdef int wordsize
    wordsize = BITS_IN_LONG

    # increase prec_in_bits to the nearest multiple of wordsize
    cdef int padded_bits
    padded_bits = (prec_in_bits + wordsize - 1) & ~(wordsize - 1)
    return int(padded_bits/wordsize + 2)

pbw = prec_bits_to_words

def prec_words_to_bits(int prec_in_words):
    r"""
    Convert from pari real precision expressed in words to precision
    expressed in bits. Note: this adjusts for the two codewords of a
    pari real, and is architecture-dependent.

    EXAMPLES::

        sage: import sage.libs.pari.gen as gen
        sage: gen.prec_words_to_bits(10)
        256   # 32-bit
        512   # 64-bit
        sage: [(n,gen.prec_words_to_bits(n)) for n in range(3,10)]
        [(3, 32), (4, 64), (5, 96), (6, 128), (7, 160), (8, 192), (9, 224)]  # 32-bit
        [(3, 64), (4, 128), (5, 192), (6, 256), (7, 320), (8, 384), (9, 448)] # 64-bit
    """
    # see user's guide to the pari library, page 10
    return int((prec_in_words - 2) * BITS_IN_LONG)

def prec_dec_to_words(int prec_in_dec):
    r"""
    Convert from precision expressed in decimal to precision expressed
    in words. Note: this rounds up to the nearest word, adjusts for the
    two codewords of a pari real, and is architecture-dependent.

    EXAMPLES::

        sage: import sage.libs.pari.gen as gen
        sage: gen.prec_dec_to_words(38)
        6   # 32-bit
        4   # 64-bit
        sage: [(n,gen.prec_dec_to_words(n)) for n in range(10,90,10)]
        [(10, 4), (20, 5), (30, 6), (40, 7), (50, 8), (60, 9), (70, 10), (80, 11)] # 32-bit
        [(10, 3), (20, 4), (30, 4), (40, 5), (50, 5), (60, 6), (70, 6), (80, 7)] # 64-bit
    """
    return prec_bits_to_words(prec_dec_to_bits(prec_in_dec))

def prec_words_to_dec(int prec_in_words):
    r"""
    Convert from precision expressed in words to precision expressed in
    decimal. Note: this adjusts for the two codewords of a pari real,
    and is architecture-dependent.

    EXAMPLES::

        sage: import sage.libs.pari.gen as gen
        sage: gen.prec_words_to_dec(5)
        28   # 32-bit
        57   # 64-bit
        sage: [(n,gen.prec_words_to_dec(n)) for n in range(3,10)]
        [(3, 9), (4, 19), (5, 28), (6, 38), (7, 48), (8, 57), (9, 67)] # 32-bit
        [(3, 19), (4, 38), (5, 57), (6, 77), (7, 96), (8, 115), (9, 134)] # 64-bit
    """
    return prec_bits_to_dec(prec_words_to_bits(prec_in_words))


# The unique running Pari instance.
cdef PariInstance pari_instance, P
pari_instance = PariInstance(16000000, 500000)
P = pari_instance   # shorthand notation

# PariInstance.__init__ must not create gen objects because their parent is not constructed yet
pari_catch_sig_on()
pari_instance.PARI_ZERO = pari_instance.new_gen_noclear(gen_0)
pari_instance.PARI_ONE  = pari_instance.new_gen_noclear(gen_1)
pari_instance.PARI_TWO  = pari_instance.new_gen_noclear(gen_2)
pari_catch_sig_off()

# Also a copy of PARI accessible from external pure python code.
pari = pari_instance

# temp variables
cdef GEN t0,t1,t2,t3,t4
t0heap = [0]*5

cdef t0GEN(x):
    global t0
    t0 = P.toGEN(x, 0)

cdef t1GEN(x):
    global t1
    t1 = P.toGEN(x, 1)

cdef t2GEN(x):
    global t2
    t2 = P.toGEN(x, 2)

cdef t3GEN(x):
    global t3
    t3 = P.toGEN(x, 3)

cdef t4GEN(x):
    global t4
    t4 = P.toGEN(x, 4)

cdef object Integer

cdef void late_import():
    global Integer

    if Integer is not None:
        return

    import sage.rings.integer
    Integer = sage.rings.integer.Integer

    global mpz_t_offset
    mpz_t_offset = sage.rings.integer.mpz_t_offset_python

cdef class gen(sage.structure.element.RingElement):
    """
    Python extension class that models the PARI GEN type.
    """
    def __init__(self):
        self.b = 0
        self._parent = P
        self._refers_to = {}

    def parent(self):
        return pari_instance

    cdef void init(self, GEN g, pari_sp b):
        """
        g - PARI GEN b - pointer to memory chunk where PARI gen lives (if
        nonzero then this memory is freed when the object goes out of
        scope)
        """
        self.g = g
        self.b = b
        self._parent = P
        self._refers_to = {}

    def __dealloc__(self):
        if self.b:
            sage_free(<void*> self.b)

    def __repr__(self):
        pari_catch_sig_on()
        return P.new_gen_to_string(self.g)

    def __hash__(self):
        """
        Return the hash of self, computed using PARI's hash_GEN().

        TESTS::

            sage: type(pari('1 + 2.0*I').__hash__())
            <type 'int'>
        """
        cdef long h
        pari_catch_sig_on()
        h = hash_GEN(self.g)
        pari_catch_sig_off()
        return h

    def _testclass(self):
        import test
        T = test.testclass()
        T._init(self)
        return T

    cdef GEN _gen(self):
        return self.g

    def list(self):
        """
        Convert self to a list of PARI gens.

        EXAMPLES:

        A PARI vector becomes a Sage list::

            sage: L = pari("vector(10,i,i^2)").list()
            sage: L
            [1, 4, 9, 16, 25, 36, 49, 64, 81, 100]
            sage: type(L)
            <type 'list'>
            sage: type(L[0])
            <type 'sage.libs.pari.gen.gen'>

        For polynomials, list() behaves as for ordinary Sage polynomials::

            sage: pol = pari("x^3 + 5/3*x"); pol.list()
            [0, 5/3, 0, 1]

        For power series or laurent series, we get all coefficients starting
        from the lowest degree term.  This includes trailing zeros::

            sage: R.<x> = LaurentSeriesRing(QQ)
            sage: s = x^2 + O(x^8)
            sage: s.list()
            [1]
            sage: pari(s).list()
            [1, 0, 0, 0, 0, 0]
            sage: s = x^-2 + O(x^0)
            sage: s.list()
            [1]
            sage: pari(s).list()
            [1, 0]

        For matrices, we get a list of columns::

            sage: M = matrix(ZZ,3,2,[1,4,2,5,3,6]); M
            [1 4]
            [2 5]
            [3 6]
            sage: pari(M).list()
            [[1, 2, 3]~, [4, 5, 6]~]

        For "scalar" types, we get a 1-element list containing ``self``::

            sage: pari("42").list()
            [42]
        """
        if typ(self.g) == t_POL:
            return list(self.Vecrev())
        return list(self.Vec())

    def __reduce__(self):
        """
        EXAMPLES::

            sage: f = pari('x^3 - 3')
            sage: loads(dumps(f)) == f
            True
        """
        s = str(self)
        import sage.libs.pari.gen_py
        return sage.libs.pari.gen_py.pari, (s,)

    cpdef ModuleElement _add_(self, ModuleElement right):
        pari_catch_sig_on()
        return P.new_gen(gadd(self.g, (<gen>right).g))

    def _add_unsafe(gen self, gen right):
        """
        VERY FAST addition of self and right on stack (and leave on stack)
        without any type checking.

        Basically, this is often about 10 times faster than just typing
        "self + right". The drawback is that (1) if self + right would give
        an error in PARI, it will totally crash Sage, and (2) the memory
        used by self + right is *never* returned - it gets allocated on
        the PARI stack and will never be freed.

        EXAMPLES::

            sage: pari(2)._add_unsafe(pari(3))
            5
        """
        global mytop
        cdef GEN z
        cdef gen w
        z = gadd(self.g, right.g)
        w = PY_NEW(gen)
        w.init(z,0)
        mytop = avma
        return w

    cpdef ModuleElement _sub_(self, ModuleElement right):
        pari_catch_sig_on()
        return P.new_gen(gsub(self.g, (<gen> right).g))

    def _sub_unsafe(gen self, gen right):
        """
        VERY FAST subtraction of self and right on stack (and leave on
        stack) without any type checking.

        Basically, this is often about 10 times faster than just typing
        "self - right". The drawback is that (1) if self - right would give
        an error in PARI, it will totally crash Sage, and (2) the memory
        used by self - right is *never* returned - it gets allocated on
        the PARI stack and will never be freed.

        EXAMPLES::

            sage: pari(2)._sub_unsafe(pari(3))
            -1
        """
        global mytop
        cdef GEN z
        cdef gen w
        z = gsub(self.g, right.g)
        w = PY_NEW(gen)
        w.init(z, 0)
        mytop = avma
        return w

    cpdef RingElement _mul_(self, RingElement right):
        pari_catch_sig_on()
        return P.new_gen(gmul(self.g, (<gen>right).g))

    def _mul_unsafe(gen self, gen right):
        """
        VERY FAST multiplication of self and right on stack (and leave on
        stack) without any type checking.

        Basically, this is often about 10 times faster than just typing
        "self \* right". The drawback is that (1) if self \* right would
        give an error in PARI, it will totally crash Sage, and (2) the
        memory used by self \* right is *never* returned - it gets
        allocated on the PARI stack and will never be freed.

        EXAMPLES::

            sage: pari(2)._mul_unsafe(pari(3))
            6
        """
        global mytop
        cdef GEN z
        cdef gen w
        z = gmul(self.g, right.g)
        w = PY_NEW(gen)
        w.init(z, 0)
        mytop = avma
        return w

    cpdef RingElement _div_(self, RingElement right):
        pari_catch_sig_on()
        return P.new_gen(gdiv(self.g, (<gen>right).g))

    def _div_unsafe(gen self, gen right):
        """
        VERY FAST division of self and right on stack (and leave on stack)
        without any type checking.

        Basically, this is often about 10 times faster than just typing
        "self / right". The drawback is that (1) if self / right would give
        an error in PARI, it will totally crash Sage, and (2) the memory
        used by self / right is *never* returned - it gets allocated on
        the PARI stack and will never be freed.

        EXAMPLES::

            sage: pari(2)._div_unsafe(pari(3))
            2/3
        """
        global mytop
        cdef GEN z
        cdef gen w
        z = gdiv(self.g, right.g)
        w = PY_NEW(gen)
        w.init(z, 0)
        mytop = avma
        return w

    #################################################################

    def _add_one(gen self):
        """
        Return self + 1.

        OUTPUT: pari gen

        EXAMPLES::

            sage: n = pari(5)
            sage: n._add_one()
            6
            sage: n = pari('x^3')
            sage: n._add_one()
            x^3 + 1
        """
        pari_catch_sig_on()
        return P.new_gen(gaddsg(1, self.g))

    def __mod__(self, other):
        cdef gen selfgen
        cdef gen othergen
        if isinstance(other, gen) and isinstance(self, gen):
            selfgen = self
            othergen = other
            pari_catch_sig_on()
            return P.new_gen(gmod(selfgen.g, othergen.g))
        return sage.structure.element.bin_op(self, other, operator.mod)

    def __pow__(self, n, m):
        t0GEN(self)
        t1GEN(n)
        pari_catch_sig_on()
        # Note: the prec parameter here has no effect when t0,t1 are
        # real; the precision of the result is the minimum of the
        # precisions of t0 and t1.  In any case the 3rd parameter to
        # gpow should be a word-precision, not a decimal precision.
        return P.new_gen(gpow(t0, t1, prec))

    def __neg__(gen self):
        pari_catch_sig_on()
        return P.new_gen(gneg(self.g))

    def __xor__(gen self, n):
        raise RuntimeError, "Use ** for exponentiation, not '^', which means xor\n"+\
              "in Python, and has the wrong precedence."

    def __rshift__(gen self, long n):
        pari_catch_sig_on()
        return P.new_gen(gshift(self.g, -n))

    def __lshift__(gen self, long n):
        pari_catch_sig_on()
        return P.new_gen(gshift(self.g, n))

    def __invert__(gen self):
        pari_catch_sig_on()
        return P.new_gen(ginv(self.g))

    ###########################################
    # ACCESS
    ###########################################
    def getattr(self, attr):
        t0GEN(str(self) + '.' + str(attr))
        pari_catch_sig_on()
        return self.new_gen(t0)

    def mod(self):
        """
        Given an INTMOD or POLMOD ``Mod(a,m)``, return the modulus `m`.

        EXAMPLES::

            sage: pari(4).Mod(5).mod()
            5
            sage: pari("Mod(x, x*y)").mod()
            y*x
            sage: pari("[Mod(4,5)]").mod()
            Traceback (most recent call last):
            ...
            TypeError: Not an INTMOD or POLMOD in mod()
        """
        if typ(self.g) != t_INTMOD and typ(self.g) != t_POLMOD:
            raise TypeError("Not an INTMOD or POLMOD in mod()")
        pari_catch_sig_on()
        # The hardcoded 1 below refers to the position in the internal
        # representation of a INTMOD or POLDMOD where the modulus is
        # stored.
        return self.new_gen(gel(self.g, 1))

    cdef GEN get_nf(self) except NULL:
        """
        Given a PARI object `self`, convert it to a proper PARI number
        field (nf) structure.

        INPUT:

        - ``self`` -- A PARI number field being the output of ``nfinit()``,
                      ``bnfinit()`` or ``bnrinit()``.

        TESTS:

        We test this indirectly through `nf_get_pol()`::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^4 - 4*x^2 + 1)
            sage: K.pari_nf().nf_get_pol()
            y^4 - 4*y^2 + 1
            sage: K.pari_bnf().nf_get_pol()
            y^4 - 4*y^2 + 1
            sage: bnr = pari("K = bnfinit(x^4 - 4*x^2 + 1); bnrinit(K, 2*x)")
            sage: bnr.nf_get_pol()
            x^4 - 4*x^2 + 1

        It does not work with ``rnfinit()`` or garbage input::

            sage: K.extension(x^2 - 5, 'b').pari_rnf().nf_get_pol()
            Traceback (most recent call last):
            ...
            TypeError: Not a PARI number field
            sage: pari("[0]").nf_get_pol()
            Traceback (most recent call last):
            ...
            TypeError: Not a PARI number field
        """
        cdef GEN nf
        cdef long nftyp
        pari_catch_sig_on()
        nf = get_nf(self.g, &nftyp)
        pari_catch_sig_off()
        if not nf:
            raise TypeError("Not a PARI number field")
        return nf

    def nf_get_pol(self):
        """
        Returns the defining polynomial of this number field.

        INPUT:

        - ``self`` -- A PARI number field being the output of ``nfinit()``,
                      ``bnfinit()`` or ``bnrinit()``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 4*x^2 + 1)
            sage: pari(K).nf_get_pol()
            y^4 - 4*y^2 + 1
            sage: bnr = pari("K = bnfinit(x^4 - 4*x^2 + 1); bnrinit(K, 2*x)")
            sage: bnr.nf_get_pol()
            x^4 - 4*x^2 + 1

        For relative extensions, we can only get the absolute polynomial,
        not the relative one::

            sage: L.<b> = K.extension(x^2 - 5)
            sage: pari(L).nf_get_pol()   # Absolute polynomial
            y^8 - 28*y^6 + 208*y^4 - 408*y^2 + 36
            sage: L.pari_rnf().nf_get_pol()
            Traceback (most recent call last):
            ...
            TypeError: Not a PARI number field
        """
        cdef GEN nf = self.get_nf()
        pari_catch_sig_on()
        return self.new_gen(nf_get_pol(nf))

    def nf_get_diff(self):
        """
        Returns the different of this number field as a PARI ideal.

        INPUT:

        - ``self`` -- A PARI number field being the output of ``nfinit()``,
                      ``bnfinit()`` or ``bnrinit()``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 4*x^2 + 1)
            sage: pari(K).nf_get_diff()
            [12, 0, 0, 0; 0, 12, 8, 0; 0, 0, 4, 0; 0, 0, 0, 4]
        """
        cdef GEN nf = self.get_nf()
        pari_catch_sig_on()
        # Very bad code, but there doesn't seem to be a better way
        return self.new_gen(gel(gel(nf, 5), 5))

    def nf_get_sign(self):
        """
        Returns a Python list ``[r1, r2]``, where ``r1`` and ``r2`` are
        Python ints representing the number of real embeddings and pairs
        of complex embeddings of this number field, respectively.

        INPUT:

        - ``self`` -- A PARI number field being the output of ``nfinit()``,
                      ``bnfinit()`` or ``bnrinit()``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 4*x^2 + 1)
            sage: s = K.pari_nf().nf_get_sign(); s
            [4, 0]
            sage: type(s); type(s[0])
            <type 'list'>
            <type 'int'>
            sage: CyclotomicField(15).pari_nf().nf_get_sign()
            [0, 4]
        """
        cdef long r1
        cdef long r2
        cdef GEN nf = self.get_nf()
        nf_get_sign(nf, &r1, &r2)
        return [r1, r2]

    def nf_get_zk(self):
        """
        Returns a vector with a `\ZZ`-basis for the ring of integers of
        this number field. The first element is always `1`.

        INPUT:

        - ``self`` -- A PARI number field being the output of ``nfinit()``,
                      ``bnfinit()`` or ``bnrinit()``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 4*x^2 + 1)
            sage: pari(K).nf_get_zk()
            [1, y, y^3 - 4*y, y^2 - 2]
        """
        cdef GEN nf = self.get_nf()
        pari_catch_sig_on()
        return self.new_gen(nf_get_zk(nf))

    def bnf_get_no(self):
        """
        Returns the class number of ``self``, a "big number field" (``bnf``).

        EXAMPLES::

            sage: K.<a> = QuadraticField(-65)
            sage: K.pari_bnf().bnf_get_no()
            8
        """
        pari_catch_sig_on()
        return self.new_gen(bnf_get_no(self.g))

    def bnf_get_cyc(self):
        """
        Returns the structure of the class group of this number field as
        a vector of SNF invariants.

        NOTE: ``self`` must be a "big number field" (``bnf``).

        EXAMPLES::

            sage: K.<a> = QuadraticField(-65)
            sage: K.pari_bnf().bnf_get_cyc()
            [4, 2]
        """
        pari_catch_sig_on()
        return self.new_gen(bnf_get_cyc(self.g))

    def bnf_get_gen(self):
        """
        Returns a vector of generators of the class group of this
        number field.

        NOTE: ``self`` must be a "big number field" (``bnf``).

        EXAMPLES::

            sage: K.<a> = QuadraticField(-65)
            sage: G = K.pari_bnf().bnf_get_gen(); G
            [[3, 2; 0, 1], [2, 1; 0, 1]]
            sage: map(lambda J: K.ideal(J), G)
            [Fractional ideal (3, a + 2), Fractional ideal (2, a + 1)]
        """
        pari_catch_sig_on()
        return self.new_gen(bnf_get_gen(self.g))

    def bnf_get_reg(self):
        """
        Returns the regulator of this number field.

        NOTE: ``self`` must be a "big number field" (``bnf``).

        EXAMPLES::

            sage: K.<a> = NumberField(x^4 - 4*x^2 + 1)
            sage: K.pari_bnf().bnf_get_reg()
            2.66089858019037...
        """
        pari_catch_sig_on()
        return self.new_gen(bnf_get_reg(self.g))

    def pr_get_p(self):
        """
        Returns the prime of `\ZZ` lying below this prime ideal.

        NOTE: ``self`` must be a PARI prime ideal (as returned by
        ``idealfactor`` for example).

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: F = pari(K).idealfactor(K.ideal(5)); F
            [[5, [-2, 1]~, 1, 1, [2, 1]~], 1; [5, [2, 1]~, 1, 1, [-2, 1]~], 1]
            sage: F[0,0].pr_get_p()
            5
        """
        pari_catch_sig_on()
        return self.new_gen(pr_get_p(self.g))

    def pr_get_e(self):
        """
        Returns the ramification index (over `\QQ`) of this prime ideal.

        NOTE: ``self`` must be a PARI prime ideal (as returned by
        ``idealfactor`` for example).

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: pari(K).idealfactor(K.ideal(2))[0,0].pr_get_e()
            2
            sage: pari(K).idealfactor(K.ideal(3))[0,0].pr_get_e()
            1
            sage: pari(K).idealfactor(K.ideal(5))[0,0].pr_get_e()
            1
        """
        cdef long e
        pari_catch_sig_on()
        e = pr_get_e(self.g)
        pari_catch_sig_off()
        return e

    def pr_get_f(self):
        """
        Returns the residue class degree (over `\QQ`) of this prime ideal.

        NOTE: ``self`` must be a PARI prime ideal (as returned by
        ``idealfactor`` for example).

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: pari(K).idealfactor(K.ideal(2))[0,0].pr_get_f()
            1
            sage: pari(K).idealfactor(K.ideal(3))[0,0].pr_get_f()
            2
            sage: pari(K).idealfactor(K.ideal(5))[0,0].pr_get_f()
            1
        """
        cdef long f
        pari_catch_sig_on()
        f = pr_get_f(self.g)
        pari_catch_sig_off()
        return f

    def pr_get_gen(self):
        """
        Returns the second generator of this PARI prime ideal, where the
        first generator is ``self.pr_get_p()``.

        NOTE: ``self`` must be a PARI prime ideal (as returned by
        ``idealfactor`` for example).

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: g = pari(K).idealfactor(K.ideal(2))[0,0].pr_get_gen(); g; K(g)
            [1, 1]~
            i + 1
            sage: g = pari(K).idealfactor(K.ideal(3))[0,0].pr_get_gen(); g; K(g)
            [3, 0]~
            3
            sage: g = pari(K).idealfactor(K.ideal(5))[0,0].pr_get_gen(); g; K(g)
            [-2, 1]~
            i - 2
        """
        pari_catch_sig_on()
        return self.new_gen(pr_get_gen(self.g))

    def bid_get_cyc(self):
        """
        Returns the structure of the group `(O_K/I)^*`, where `I` is the
        ideal represented by ``self``.

        NOTE: ``self`` must be a "big ideal" (``bid``) as returned by
        ``idealstar`` for example.

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: J = pari(K).idealstar(K.ideal(4*i + 2))
            sage: J.bid_get_cyc()
            [4, 2]
        """
        pari_catch_sig_on()
        return self.new_gen(bid_get_cyc(self.g))

    def bid_get_gen(self):
        """
        Returns a vector of generators of the group `(O_K/I)^*`, where
        `I` is the ideal represented by ``self``.

        NOTE: ``self`` must be a "big ideal" (``bid``) with generators,
        as returned by ``idealstar`` with ``flag`` = 2.

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: J = pari(K).idealstar(K.ideal(4*i + 2), 2)
            sage: J.bid_get_gen()
            [7, [-2, -1]~]

        We get an exception if we do not supply ``flag = 2`` to
        ``idealstar``::

            sage: J = pari(K).idealstar(K.ideal(3))
            sage: J.bid_get_gen()
            Traceback (most recent call last):
            ...
            PariError:  (5)
        """
        pari_catch_sig_on()
        return self.new_gen(bid_get_gen(self.g))

    def __getitem__(gen self, n):
        """
        Return the nth entry of self. The indexing is 0-based, like in
        Python. Note that this is *different* than the default behavior
        of the PARI/GP interpreter.

        EXAMPLES::

            sage: p = pari('1 + 2*x + 3*x^2')
            sage: p[0]
            1
            sage: p[2]
            3
            sage: p[100]
            0
            sage: p[-1]
            0
            sage: q = pari('x^2 + 3*x^3 + O(x^6)')
            sage: q[3]
            3
            sage: q[5]
            0
            sage: q[6]
            Traceback (most recent call last):
            ...
            IndexError: index out of bounds
            sage: m = pari('[1,2;3,4]')
            sage: m[0]
            [1, 3]~
            sage: m[1,0]
            3
            sage: l = pari('List([1,2,3])')
            sage: l[1]
            2
            sage: s = pari('"hello, world!"')
            sage: s[0]
            'h'
            sage: s[4]
            'o'
            sage: s[12]
            '!'
            sage: s[13]
            Traceback (most recent call last):
            ...
            IndexError: index out of bounds
            sage: v = pari('[1,2,3]')
            sage: v[0]
            1
            sage: c = pari('Col([1,2,3])')
            sage: c[1]
            2
            sage: sv = pari('Vecsmall([1,2,3])')
            sage: sv[2]
            3
            sage: type(sv[2])
            <type 'int'>
            sage: tuple(pari(3/5))
            (3, 5)
            sage: tuple(pari('1 + 5*I'))
            (1, 5)
            sage: tuple(pari('Qfb(1, 2, 3)'))
            (1, 2, 3)
            sage: pari(57)[0]
            Traceback (most recent call last):
            ...
            TypeError: unindexable object
            sage: m = pari("[[1,2;3,4],5]") ; m[0][1,0]
            3
            sage: v = pari(xrange(20))
            sage: v[2:5]
            [2, 3, 4]
            sage: v[:]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: v[10:]
            [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: v[:5]
            [0, 1, 2, 3, 4]
            sage: v
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: v[-1]
            Traceback (most recent call last):
            ...
            IndexError: index out of bounds
            sage: v[:-3]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
            sage: v[5:]
            [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: pari([])[::]
            []
        """
        cdef int pari_type

        pari_type = typ(self.g)

        if isinstance(n, tuple):
            if pari_type != t_MAT:
                raise TypeError, "self must be of pari type t_MAT"
            if len(n) != 2:
                raise IndexError, "index must be an integer or a 2-tuple (i,j)"
            i = int(n[0])
            j = int(n[1])

            if i < 0 or i >= glength(<GEN>(self.g[1])):
                raise IndexError, "row index out of bounds"
            if j < 0 or j >= glength(self.g):
                raise IndexError, "column index out of bounds"

            ind = (i,j)

            if PyDict_Contains(self._refers_to, ind):
                return self._refers_to[ind]
            else:
                ## In this case, we're being asked to return
                ## a GEN that has no gen pointing to it, so
                ## we need to create such a gen, add it to
                ## self._refers_to, and return it.
                val = P.new_ref(gmael(self.g, j+1, i+1), self)
                self._refers_to[ind] = val
                return val

        elif isinstance(n, slice):
            l = glength(self.g)
            start,stop,step = n.indices(l)
            inds = xrange(start,stop,step)
            k = len(inds)
            # fast exit
            if k==0:
                return P.vector(0)
            # fast call, beware pari is one based
            if pari_type == t_VEC:
                if step==1:
                    return self.vecextract('"'+str(start+1)+".."+str(stop)+'"')
                if step==-1:
                    return self.vecextract('"'+str(start+1)+".."+str(stop+2)+'"')
            # slow call
            v = P.vector(k)
            for i,j in enumerate(inds):
                v[i] = self[j]
            return v

        ## there are no "out of bounds" problems
        ## for a polynomial or power series, so these go before
        ## bounds testing
        if pari_type == t_POL:
            return self.polcoeff(n)

        elif pari_type == t_SER:
            bound = valp(self.g) + lg(self.g) - 2
            if n >= bound:
                raise IndexError, "index out of bounds"
            return self.polcoeff(n)

        elif pari_type in (t_INT, t_REAL, t_PADIC, t_QUAD):
            # these are definitely scalar!
            raise TypeError, "unindexable object"

        elif n < 0 or n >= glength(self.g):
            raise IndexError, "index out of bounds"

        elif pari_type == t_VEC or pari_type == t_MAT:
            #t_VEC    : row vector        [ code ] [  x_1  ] ... [  x_k  ]
            #t_MAT    : matrix            [ code ] [ col_1 ] ... [ col_k ]
            if PyDict_Contains(self._refers_to, n):
                return self._refers_to[n]
            else:
                ## In this case, we're being asked to return
                ## a GEN that has no gen pointing to it, so
                ## we need to create such a gen, add it to
                ## self._refers_to, and return it.
                val = P.new_ref(gel(self.g, n+1), self)
                self._refers_to[n] = val
                return val

        elif pari_type == t_VECSMALL:
            #t_VECSMALL: vec. small ints  [ code ] [ x_1 ] ... [ x_k ]
            return self.g[n+1]

        elif pari_type == t_STR:
            #t_STR    : string            [ code ] [ man_1 ] ... [ man_k ]
            return chr( (<char *>(self.g+1))[n] )

        elif pari_type == t_LIST:
            #t_LIST   : list              [ code ] [ n ] [ nmax ][ vec ]
            return self.component(n+1)
            #code from previous version, now segfaults:
            #return P.new_ref(gel(self.g,n+2), self)

        elif pari_type in (t_INTMOD, t_POLMOD):
            #t_INTMOD : integermods       [ code ] [ mod  ] [ integer ]
            #t_POLMOD : poly mod          [ code ] [ mod  ] [ polynomial ]

            # if we keep going we would get:
            #   [0] = modulus
            #   [1] = lift to t_INT or t_POL
            # do we want this? maybe the other way around?
            raise TypeError, "unindexable object"

        #elif pari_type in (t_FRAC, t_RFRAC):
            # generic code gives us:
            #   [0] = numerator
            #   [1] = denominator

        #elif pari_type == t_COMPLEX:
            # generic code gives us
            #   [0] = real part
            #   [1] = imag part

        #elif type(self.g) in (t_QFR, t_QFI):
            # generic code works ok

        else:
            ## generic code, which currently handles cases
            ## as mentioned above
            return P.new_ref(gel(self.g,n+1), self)

    def _gen_length(gen self):
        """
        Return the length of self as a Pari object, *including*
        codewords.

        EXAMPLES::

            sage: n = pari(30)
            sage: n.length()
            1
            sage: n._gen_length()
            3
        """
        return lg(self.g)

    def __setitem__(gen self, n, y):
        r"""
        Set the nth entry to a reference to y.


            -  The indexing is 0-based, like everywhere else in Python, but
               *unlike* in Pari/GP.

            -  Assignment sets the nth entry to a reference to y, assuming y is
               an object of type gen. This is the same as in Python, but
               *different* than what happens in the gp interpreter, where
               assignment makes a copy of y.

            -  Because setting creates references it is *possible* to make
               circular references, unlike in GP. Do *not* do this (see the
               example below). If you need circular references, work at the Python
               level (where they work well), not the PARI object level.



        EXAMPLES::

            sage: v = pari(range(10))
            sage: v
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: v[0] = 10
            sage: w = pari([5,8,-20])
            sage: v
            [10, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: v[1] = w
            sage: v
            [10, [5, 8, -20], 2, 3, 4, 5, 6, 7, 8, 9]
            sage: w[0] = -30
            sage: v
            [10, [-30, 8, -20], 2, 3, 4, 5, 6, 7, 8, 9]
            sage: t = v[1]; t[1] = 10 ; v
            [10, [-30, 10, -20], 2, 3, 4, 5, 6, 7, 8, 9]
            sage: v[1][0] = 54321 ; v
            [10, [54321, 10, -20], 2, 3, 4, 5, 6, 7, 8, 9]
            sage: w
            [54321, 10, -20]
            sage: v = pari([[[[0,1],2],3],4]) ; v[0][0][0][1] = 12 ; v
            [[[[0, 12], 2], 3], 4]
            sage: m = pari(matrix(2,2,range(4))) ; l = pari([5,6]) ; n = pari(matrix(2,2,[7,8,9,0])) ; m[1,0] = l ; l[1] = n ; m[1,0][1][1,1] = 1111 ; m
            [0, 1; [5, [7, 8; 9, 1111]], 3]
            sage: m = pari("[[1,2;3,4],5,6]") ; m[0][1,1] = 11 ; m
            [[1, 2; 3, 11], 5, 6]

        Finally, we create a circular reference::

            sage: v = pari([0])
            sage: w = pari([v])
            sage: v
            [0]
            sage: w
            [[0]]
            sage: v[0] = w

        Now there is a circular reference. Accessing v[0] will crash Sage.

        ::

            sage: s=pari.vector(2,[0,0])
            sage: s[:1]
            [0]
            sage: s[:1]=[1]
            sage: s
            [1, 0]
            sage: type(s[0])
            <type 'sage.libs.pari.gen.gen'>
            sage: s = pari(range(20)) ; s
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: s[0:10:2] = range(50,55) ; s
            [50, 1, 51, 3, 52, 5, 53, 7, 54, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: s[10:20:3] = range(100,150) ; s
            [50, 1, 51, 3, 52, 5, 53, 7, 54, 9, 100, 11, 12, 101, 14, 15, 102, 17, 18, 103]

        TESTS::

            sage: v = pari(xrange(10)) ; v
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: v[:] = [20..29]
            sage: v
            [20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
            sage: type(v[0])
            <type 'sage.libs.pari.gen.gen'>
        """
        cdef int i, j
        cdef gen x
        cdef long l
        cdef Py_ssize_t ii, jj, step

        pari_catch_sig_on()
        try:
            if isinstance(y, gen):
                x = y
            else:
                x = pari(y)

            if isinstance(n, tuple):
                if typ(self.g) != t_MAT:
                    raise TypeError, "cannot index Pari type %s by tuple"%typ(self.g)

                if len(n) != 2:
                    raise ValueError, "matrix index must be of the form [row, column]"

                i = int(n[0])
                j = int(n[1])
                ind = (i,j)

                if i < 0 or i >= glength(<GEN>(self.g[1])):
                    raise IndexError, "row i(=%s) must be between 0 and %s"%(i,self.nrows()-1)
                if j < 0 or j >= glength(self.g):
                    raise IndexError, "column j(=%s) must be between 0 and %s"%(j,self.ncols()-1)
                self._refers_to[ind] = x

                (<GEN>(self.g)[j+1])[i+1] = <long>(x.g)
                return

            elif isinstance(n, slice):
                l = glength(self.g)
                inds = xrange(*n.indices(l))
                k = len(inds)
                if k > len(y):
                    raise ValueError, "attempt to assign sequence of size %s to slice of size %s"%(len(y), k)

                # actually set the values
                for i,j in enumerate(inds):
                    self[j] = y[i]
                return

            i = int(n)

            if i < 0 or i >= glength(self.g):
                raise IndexError, "index (%s) must be between 0 and %s"%(i,glength(self.g)-1)

            # so python memory manager will work correctly
            # and not free x if PARI part of self is the
            # only thing pointing to it.
            self._refers_to[i] = x

            ## correct indexing for t_POLs
            if typ(self.g) == t_POL:
                i = i + 1

            ## actually set the value
            (self.g)[i+1] = <long>(x.g)
            return
        finally:
            pari_catch_sig_off()

    def __len__(gen self):
        return glength(self.g)



    ###########################################
    # comparisons
    # I had to rewrite PARI's compare, since
    # otherwise trapping signals and other horrible,
    # memory-leaking and slow stuff occurs.
    ###########################################

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Comparisons

        First uses PARI's cmp routine; if it decides the objects are not
        comparable, it then compares the underlying strings (since in
        Python all objects are supposed to be comparable).

        EXAMPLES::

            sage: a = pari(5)
            sage: b = 10
            sage: a < b
            True
            sage: a <= b
            True
            sage: a <= 5
            True
            sage: a > b
            False
            sage: a >= b
            False
            sage: a >= pari(10)
            False
            sage: a == 5
            True
            sage: a is 5
            False

        ::

            sage: pari(2.5) > None
            True
            sage: pari(3) == pari(3)
            True
            sage: pari('x^2 + 1') == pari('I-1')
            False
            sage: pari(I) == pari(I)
            True
        """
        return gcmp_sage(left.g, (<gen>right).g)

    def __copy__(gen self):
        pari_catch_sig_on()
        return P.new_gen(gcopy(self.g))

    ###########################################
    # Conversion --> Python
    # Try to convert to a meaningful python object
    # in various ways
    ###########################################

    def list_str(gen self):
        """
        Return str that might correctly evaluate to a Python-list.
        """
        s = str(self)
        if s[:4] == "Mat(":
            s = "[" + s[4:-1] + "]"
        s = s.replace("~","")
        if s.find(";") != -1:
            s = s.replace(";","], [")
            s = "[" + s + "]"
            return eval(s)
        else:
            return eval(s)

    def __hex__(gen self):
        """
        Return the hexadecimal digits of self in lower case.

        EXAMPLES::

            sage: print hex(pari(0))
            0
            sage: print hex(pari(15))
            f
            sage: print hex(pari(16))
            10
            sage: print hex(pari(16938402384092843092843098243))
            36bb1e3929d1a8fe2802f083
            sage: print hex(long(16938402384092843092843098243))
            0x36bb1e3929d1a8fe2802f083L
            sage: print hex(pari(-16938402384092843092843098243))
            -36bb1e3929d1a8fe2802f083
        """
        cdef GEN x
        cdef long lx, *xp
        cdef long w
        cdef char *s, *sp
        cdef char *hexdigits
        hexdigits = "0123456789abcdef"
        cdef int i, j
        cdef int size
        x = self.g
        if typ(x) != t_INT:
            raise TypeError, "gen must be of PARI type t_INT"
        if not signe(x):
            return "0"
        lx = lgefint(x)-2  # number of words
        size = lx*2*sizeof(long)
        s = <char *>sage_malloc(size+2) # 1 char for sign, 1 char for '\0'
        sp = s + size+1
        sp[0] = 0
        xp = int_LSW(x)
        for i from 0 <= i < lx:
            w = xp[0]
            for j from 0 <= j < 2*sizeof(long):
                sp = sp-1
                sp[0] = hexdigits[w & 15]
                w = w>>4
            xp = int_nextW(xp)
        # remove leading zeros!
        while sp[0] == c'0':
            sp = sp+1
        if signe(x) < 0:
            sp = sp-1
            sp[0] = c'-'
        k = <object>sp
        sage_free(s)
        return k

    def __int__(gen self):
        """
        Return Python int. Very fast, and if the number is too large to fit
        into a C int, a Python long is returned instead.

        EXAMPLES::

            sage: int(pari(0))
            0
            sage: int(pari(10))
            10
            sage: int(pari(-10))
            -10
            sage: int(pari(123456789012345678901234567890))
            123456789012345678901234567890L
            sage: int(pari(-123456789012345678901234567890))
            -123456789012345678901234567890L
            sage: int(pari(2^31-1))
            2147483647
            sage: int(pari(-2^31))
            -2147483648
            sage: int(pari("Pol(10)"))
            10
        """
        cdef GEN x
        cdef long lx, *xp
        if  typ(self.g)==t_POL and self.poldegree()<=0:
            # Change a constant polynomial to its constant term
            x = constant_term(self.g)
        else:
            x = self.g
        if typ(x) != t_INT:
            raise TypeError, "gen must be of PARI type t_INT or t_POL of degree 0"
        if not signe(x):
            return 0
        lx = lgefint(x)-3   # take 1 to account for the MSW
        xp = int_MSW(x)
        # special case 1 word so we return int if it fits
        if not lx:
            if   signe(x) |  xp[0] > 0:     # both positive
                return xp[0]
            elif signe(x) & -xp[0] < 0:     # both negative
                return -xp[0]
        i = <ulong>xp[0]
        while lx:
            xp = int_precW(xp)
            i = i << BITS_IN_LONG | <ulong>xp[0]
            lx = lx-1
        if signe(x) > 0:
            return i
        else:
            return -i
        # NOTE: Could use int_unsafe below, which would be much faster, but
        # the default PARI prints annoying stuff to the screen when
        # the number is large.

    def int_unsafe(gen self):
        """
        Returns int form of self, but raises an exception if int does not
        fit into a long integer.

        This is about 5 times faster than the usual int conversion.
        """
        return gtolong(self.g)

    def intvec_unsafe(self):
        """
        Returns Python int list form of entries of self, but raises an
        exception if int does not fit into a long integer. Here self must
        be a vector.

        EXAMPLES::

            sage: pari('[3,4,5]').type()
            't_VEC'
            sage: pari('[3,4,5]').intvec_unsafe()
            [3, 4, 5]
            sage: type(pari('[3,4,5]').intvec_unsafe()[0])
            <type 'int'>

        TESTS::

            sage: pari(3).intvec_unsafe()
            Traceback (most recent call last):
            ...
            TypeError: gen must be of PARI type t_VEC
            sage: pari('[2^150,1]').intvec_unsafe()
            Traceback (most recent call last):
            ...
            PariError:  (15)
        """
        cdef int n, L
        cdef object v
        cdef GEN g
        g = self.g
        if typ(g) != t_VEC:
            raise TypeError, "gen must be of PARI type t_VEC"

        pari_catch_sig_on()
        L = glength(g)
        v = []
        for n from 0 <= n < L:
            v.append(gtolong(<GEN> (g[n+1])))
        pari_catch_sig_off()
        return v

    def python_list_small(gen self):
        """
        Return a Python list of the PARI gens. This object must be of type
        t_VECSMALL, and the resulting list contains python 'int's

        EXAMPLES::

            sage: v=pari([1,2,3,10,102,10]).Vecsmall()
            sage: w = v.python_list_small()
            sage: w
            [1, 2, 3, 10, 102, 10]
            sage: type(w[0])
            <type 'int'>
        """
        cdef long n
        if typ(self.g) != t_VECSMALL:
            raise TypeError, "Object (=%s) must be of type t_VECSMALL."%self
        V = []
        m = glength(self.g)
        for n from 0 <= n < m:
            V.append(self.g[n+1])
        return V

    def python_list(gen self):
        """
        Return a Python list of the PARI gens. This object must be of type
        t_VEC

        INPUT: NoneOUTPUT:


        -  ``list`` - Python list whose elements are the
           elements of the input gen.


        EXAMPLES::

            sage: v=pari([1,2,3,10,102,10])
            sage: w = v.python_list()
            sage: w
            [1, 2, 3, 10, 102, 10]
            sage: type(w[0])
            <type 'sage.libs.pari.gen.gen'>
            sage: pari("[1,2,3]").python_list()
            [1, 2, 3]
        """
        cdef long n, m
        cdef gen t

        if typ(self.g) != t_VEC:
            raise TypeError, "Object (=%s) must be of type t_VEC."%self
        m = glength(self.g)
        V = []
        for n from 0 <= n < m:
##            t = P.new_ref(<GEN> (self.g[n+1]), V)
##            V.append(t)
            V.append(self.__getitem__(n))
        return V

    def python(self, locals=None):
        """
        Return Python eval of self.

        Note: is self is a real (type t_REAL) the result will be a
        RealField element of the equivalent precision; if self is a complex
        (type t_COMPLEX) the result will be a ComplexField element of
        precision the minimum precision of the real and imaginary parts.

        EXAMPLES::

            sage: pari('389/17').python()
            389/17
            sage: f = pari('(2/3)*x^3 + x - 5/7 + y'); f
            2/3*x^3 + x + (y - 5/7)
            sage: var('x,y')
            (x, y)
            sage: f.python({'x':x, 'y':y})
            2/3*x^3 + x + y - 5/7

        You can also use .sage, which is a psynonym::

            sage: f.sage({'x':x, 'y':y})
            2/3*x^3 + x + y - 5/7
        """
        import sage.libs.pari.gen_py
        return sage.libs.pari.gen_py.python(self, locals=locals)

    sage = python

#  This older version illustrates how irrelevant the real_precision
#  global variable is in the pari library.  The output of this
#  function when self is a t_REAL is completely independent of the
#  precision parameter.  The new version above has no precision
#  parameter at all: the caller can coerce into a lower-precision
#  RealField if desired (or even in to a higher precision one, but
#  that will just pad with 0 bits).

#     def python(self, precision=0, bits_prec=None):
#         """
#         Return Python eval of self.
#         """
#         if not bits_prec is None:
#             precision = int(bits_prec * 3.4) + 1
#         import sage.libs.pari.gen_py
#         cdef long orig_prec
#         if precision:
#             orig_prec = P.get_real_precision()
#             P.set_real_precision(precision)
#             print "self.precision() = ",self.precision()
#             x = sage.libs.pari.gen_py.python(self)
#             print "x.prec() = ",x.prec()
#             P.set_real_precision(orig_prec)
#             return x
#         else:
#             return sage.libs.pari.gen_py.python(self)

    _sage_ = _eval_ = python

    def __long__(gen self):
        """
        Return Python long.
        """
        return long(int(self))

    def __float__(gen self):
        """
        Return Python float.
        """
        cdef double d
        pari_catch_sig_on()
        d = gtodouble(self.g)
        pari_catch_sig_off()
        return d

    def __complex__(self):
        r"""
        Return ``self`` as a Python ``complex``
        value.

        EXAMPLES::

            sage: g = pari(-1.0)^(1/5); g
            0.809016994374947 + 0.587785252292473*I
            sage: g.__complex__()
            (0.8090169943749475+0.5877852522924731j)
            sage: complex(g)
            (0.8090169943749475+0.5877852522924731j)

        ::

            sage: g = pari(Integers(5)(3)); g
            Mod(3, 5)
            sage: complex(g)
            Traceback (most recent call last):
            ...
            PariError: incorrect type (11)
        """
        cdef double re, im
        pari_catch_sig_on()
        re = gtodouble(greal(self.g))
        im = gtodouble(gimag(self.g))
        pari_catch_sig_off()
        return complex(re, im)

    def __nonzero__(self):
        """
        EXAMPLES::

            sage: pari('1').__nonzero__()
            True
            sage: pari('x').__nonzero__()
            True
            sage: bool(pari(0))
            False
            sage: a = pari('Mod(0,3)')
            sage: a.__nonzero__()
            False
        """
        return not gequal0(self.g)


    ###########################################
    # Comparisons (from PARI)
    ###########################################

    def gequal(gen a, b):
        r"""
        Check whether `a` and `b` are equal using PARI's ``gequal``.

        EXAMPLES::

            sage: a = pari(1); b = pari(1.0); c = pari('"some_string"')
            sage: a.gequal(a)
            True
            sage: b.gequal(b)
            True
            sage: c.gequal(c)
            True
            sage: a.gequal(b)
            True
            sage: a.gequal(c)
            False

        WARNING: this relation is not transitive::

            sage: a = pari('[0]'); b = pari(0); c = pari('[0,0]')
            sage: a.gequal(b)
            True
            sage: b.gequal(c)
            True
            sage: a.gequal(c)
            False
        """
        t0GEN(b)
        pari_catch_sig_on()
        cdef int ret = gequal(a.g, t0)
        pari_catch_sig_off()
        return ret != 0

    def gequal0(gen a):
        r"""
        Check whether `a` is equal to zero.

        EXAMPLES::

            sage: pari(0).gequal0()
            True
            sage: pari(1).gequal0()
            False
            sage: pari(1e-100).gequal0()
            False
            sage: pari("0.0 + 0.0*I").gequal0()
            True
            sage: pari(GF(3^20,'t')(0)).gequal0()
            True
        """
        pari_catch_sig_on()
        cdef int ret = gequal0(a.g)
        pari_catch_sig_off()
        return ret != 0

    def gequal_long(gen a, long b):
        r"""
        Check whether `a` is equal to the ``long int`` `b` using PARI's ``gequalsg``.

        EXAMPLES::

            sage: a = pari(1); b = pari(2.0); c = pari('3*matid(3)')
            sage: a.gequal_long(1)
            True
            sage: a.gequal_long(-1)
            False
            sage: a.gequal_long(0)
            False
            sage: b.gequal_long(2)
            True
            sage: b.gequal_long(-2)
            False
            sage: c.gequal_long(3)
            True
            sage: c.gequal_long(-3)
            False
        """
        pari_catch_sig_on()
        cdef int ret = gequalsg(b, a.g)
        pari_catch_sig_off()
        return ret != 0


    ###########################################
    # arith1.c
    ###########################################
    def isprime(gen self, flag=0):
        """
        isprime(x, flag=0): Returns True if x is a PROVEN prime number, and
        False otherwise.

        INPUT:


        -  ``flag`` - int 0 (default): use a combination of
           algorithms. 1: certify primality using the Pocklington-Lehmer Test.
           2: certify primality using the APRCL test.


        OUTPUT:


        -  ``bool`` - True or False


        EXAMPLES::

            sage: pari(9).isprime()
            False
            sage: pari(17).isprime()
            True
            sage: n = pari(561)    # smallest Carmichael number
            sage: n.isprime()      # not just a pseudo-primality test!
            False
            sage: n.isprime(1)
            False
            sage: n.isprime(2)
            False
        """
        cdef bint t
        pari_catch_sig_on()
        t = (signe(gisprime(self.g, flag)) != 0)
        pari_catch_sig_off()
        return t

    def qfbhclassno(gen n):
        r"""
        Computes the Hurwitz-Kronecker class number of `n`.

        INPUT:

        - `n` (gen) -- a non-negative integer

        .. note::

           If `n` is large (more than `5*10^5`), the result is
           conditional upon GRH.

        EXAMPLES:

        The Hurwitx class number is 0 is n is congruent to 1 or 2 modulo 4::
            sage: pari(-10007).qfbhclassno()
            0
            sage: pari(-2).qfbhclassno()
            0

        It is -1/12 for n=0::

            sage: pari(0).qfbhclassno()
            -1/12

        Otherwise it is the number of classes of positive definite
        binary quadratic forms with discriminant `-n`, weighted by
        `1/m` where `m` is the number of automorphisms of the form::

            sage: pari(4).qfbhclassno()
            1/2
            sage: pari(3).qfbhclassno()
            1/3
            sage: pari(23).qfbhclassno()
            3

        """
        pari_catch_sig_on()
        return P.new_gen(hclassno(n.g))

    def ispseudoprime(gen self, flag=0):
        """
        ispseudoprime(x, flag=0): Returns True if x is a pseudo-prime
        number, and False otherwise.

        INPUT:


        -  ``flag`` - int 0 (default): checks whether x is a
           Baillie-Pomerance-Selfridge-Wagstaff pseudo prime (strong
           Rabin-Miller pseudo prime for base 2, followed by strong Lucas test
           for the sequence (P,-1), P smallest positive integer such that
           `P^2 - 4` is not a square mod x). 0: checks whether x is a
           strong Miller-Rabin pseudo prime for flag randomly chosen bases
           (with end-matching to catch square roots of -1).


        OUTPUT:


        -  ``bool`` - True or False


        EXAMPLES::

            sage: pari(9).ispseudoprime()
            False
            sage: pari(17).ispseudoprime()
            True
            sage: n = pari(561)     # smallest Carmichael number
            sage: n.ispseudoprime(2)
            False
        """
        cdef long z
        pari_catch_sig_on()
        z = ispseudoprime(self.g, flag)
        pari_catch_sig_off()
        return (z != 0)

    def ispower(gen self, k=None):
        r"""
        Determine whether or not self is a perfect k-th power. If k is not
        specified, find the largest k so that self is a k-th power.

        .. note:::

           There is a BUG in the PARI C-library function (at least in
           PARI 2.2.12-beta) that is used to implement this function! This is
           in GP::

              ? p=nextprime(10^100); n=p^100; m=p^2; m^50==n; ispower(n,50)


        INPUT:


        -  ``k`` - int (optional)


        OUTPUT:


        -  ``power`` - int, what power it is

        -  ``g`` - what it is a power of


        EXAMPLES::

            sage: pari(9).ispower()
            (2, 3)
            sage: pari(17).ispower()
            (1, 17)
            sage: pari(17).ispower(2)
            (False, None)
            sage: pari(17).ispower(1)
            (1, 17)
            sage: pari(2).ispower()
            (1, 2)
        """
        cdef int n
        cdef GEN x

        if k is None:
            pari_catch_sig_on()
            n = gisanypower(self.g, &x)
            if n == 0:
                pari_catch_sig_off()
                return 1, self
            else:
                return n, P.new_gen(x)
        else:
            k = int(k)
            t0GEN(k)
            pari_catch_sig_on()
            n = ispower(self.g, t0, &x)
            if n == 0:
                pari_catch_sig_off()
                return False, None
            else:
                return k, P.new_gen(x)

    ###########################################
    # 1: Standard monadic or dyadic OPERATORS
    ###########################################
    def divrem(gen x, y, var=-1):
        """
        divrem(x, y, v): Euclidean division of x by y giving as a
        2-dimensional column vector the quotient and the remainder, with
        respect to v (to main variable if v is omitted).
        """
        t0GEN(y)
        pari_catch_sig_on()
        return P.new_gen(divrem(x.g, t0, P.get_var(var)))

    def lex(gen x, y):
        """
        lex(x,y): Compare x and y lexicographically (1 if xy, 0 if x==y, -1
        if xy)
        """
        t0GEN(y)
        pari_catch_sig_on()
        return lexcmp(x.g, t0)

    def max(gen x, y):
        """
        max(x,y): Return the maximum of x and y.
        """
        t0GEN(y)
        pari_catch_sig_on()
        return P.new_gen(gmax(x.g, t0))

    def min(gen x, y):
        """
        min(x,y): Return the minimum of x and y.
        """
        t0GEN(y)
        pari_catch_sig_on()
        return P.new_gen(gmin(x.g, t0))

    def shift(gen x, long n):
        """
        shift(x,n): shift x left n bits if n=0, right -n bits if n0.
        """
        pari_catch_sig_on()
        return P.new_gen(gshift(x.g, n))

    def shiftmul(gen x, long n):
        """
        shiftmul(x,n): Return the product of x by `2^n`.
        """
        pari_catch_sig_on()
        return P.new_gen(gmul2n(x.g, n))

    def moebius(gen x):
        """
        moebius(x): Moebius function of x.
        """
        pari_catch_sig_on()
        return P.new_gen(gmoebius(x.g))

    def sign(gen x):
        """
        sign(x): Return the sign of x, where x is of type integer, real or
        fraction.
        """
        # Pari throws an error if you attempt to take the sign of
        # a complex number.
        pari_catch_sig_on()
        return gsigne(x.g)

    def vecmax(gen x):
        """
        vecmax(x): Return the maximum of the elements of the vector/matrix
        x.
        """
        pari_catch_sig_on()
        return P.new_gen(vecmax(x.g))


    def vecmin(gen x):
        """
        vecmin(x): Return the maximum of the elements of the vector/matrix
        x.
        """
        pari_catch_sig_on()
        return P.new_gen(vecmin(x.g))



    ###########################################
    # 2: CONVERSIONS and similar elementary functions
    ###########################################

    def Col(gen x, long n = 0):
        """
        Transform the object `x` into a column vector with minimal size `|n|`.

        INPUT:

        - ``x`` -- gen

        - ``n`` -- Make the column vector of minimal length `|n|`. If `n > 0`,
          append zeros; if `n < 0`, prepend zeros.

        OUTPUT:

        A PARI column vector (type ``t_COL``)

        EXAMPLES::

            sage: pari(1.5).Col()
            [1.50000000000000]~
            sage: pari([1,2,3,4]).Col()
            [1, 2, 3, 4]~
            sage: pari('[1,2; 3,4]').Col()
            [[1, 2], [3, 4]]~
            sage: pari('"Sage"').Col()
            ["S", "a", "g", "e"]~
            sage: pari('x + 3*x^3').Col()
            [3, 0, 1, 0]~
            sage: pari('x + 3*x^3 + O(x^5)').Col()
            [1, 0, 3, 0]~

        We demonstate the `n` argument::

            sage: pari([1,2,3,4]).Col(2)
            [1, 2, 3, 4]~
            sage: pari([1,2,3,4]).Col(-2)
            [1, 2, 3, 4]~
            sage: pari([1,2,3,4]).Col(6)
            [1, 2, 3, 4, 0, 0]~
            sage: pari([1,2,3,4]).Col(-6)
            [0, 0, 1, 2, 3, 4]~

        See also :meth:`Vec` (create a row vector) for more examples
        and :meth:`Colrev` (create a column in reversed order).
        """
        pari_catch_sig_on()
        return P.new_gen(_Vec_append(gtocol(x.g), gen_0, n))

    def Colrev(gen x, long n = 0):
        """
        Transform the object `x` into a column vector with minimal size `|n|`.
        The order of the resulting vector is reversed compared to :meth:`Col`.

        INPUT:

        - ``x`` -- gen

        - ``n`` -- Make the vector of minimal length `|n|`. If `n > 0`,
          prepend zeros; if `n < 0`, append zeros.

        OUTPUT:

        A PARI column vector (type ``t_COL``)

        EXAMPLES::

            sage: pari(1.5).Colrev()
            [1.50000000000000]~
            sage: pari([1,2,3,4]).Colrev()
            [4, 3, 2, 1]~
            sage: pari('[1,2; 3,4]').Colrev()
            [[3, 4], [1, 2]]~
            sage: pari('x + 3*x^3').Colrev()
            [0, 1, 0, 3]~

        We demonstate the `n` argument::

            sage: pari([1,2,3,4]).Colrev(2)
            [4, 3, 2, 1]~
            sage: pari([1,2,3,4]).Colrev(-2)
            [4, 3, 2, 1]~
            sage: pari([1,2,3,4]).Colrev(6)
            [0, 0, 4, 3, 2, 1]~
            sage: pari([1,2,3,4]).Colrev(-6)
            [4, 3, 2, 1, 0, 0]~
        """
        pari_catch_sig_on()
        # Create a non-reversed column vector
        cdef GEN v = _Vec_append(gtocol(x.g), gen_0, n)
        # Reverse it in-place
        cdef GEN L = v + 1
        cdef GEN R = v + (lg(v)-1)
        cdef long t
        while (L < R):
            t = L[0]
            L[0] = R[0]
            R[0] = t
            L += 1
            R -= 1
        return P.new_gen(v)

    def List(gen x):
        """
        List(x): transforms the PARI vector or list x into a list.

        EXAMPLES::

            sage: v = pari([1,2,3])
            sage: v
            [1, 2, 3]
            sage: v.type()
            't_VEC'
            sage: w = v.List()
            sage: w
            List([1, 2, 3])
            sage: w.type()
            't_LIST'
        """
        pari_catch_sig_on()
        return P.new_gen(gtolist(x.g))

    def Mat(gen x):
        """
        Mat(x): Returns the matrix defined by x.

        - If x is already a matrix, a copy of x is created and returned.

        - If x is not a vector or a matrix, this function returns a 1x1
          matrix.

        - If x is a row (resp. column) vector, this functions returns
          a 1-row (resp. 1-column) matrix, *unless* all elements are
          column (resp. row) vectors of the same length, in which case
          the vectors are concatenated sideways and the associated big
          matrix is returned.

        INPUT:


        -  ``x`` - gen


        OUTPUT:


        -  ``gen`` - a PARI matrix


        EXAMPLES::

            sage: x = pari(5)
            sage: x.type()
            't_INT'
            sage: y = x.Mat()
            sage: y
            Mat(5)
            sage: y.type()
            't_MAT'
            sage: x = pari('[1,2;3,4]')
            sage: x.type()
            't_MAT'
            sage: x = pari('[1,2,3,4]')
            sage: x.type()
            't_VEC'
            sage: y = x.Mat()
            sage: y
            Mat([1, 2, 3, 4])
            sage: y.type()
            't_MAT'

        ::

            sage: v = pari('[1,2;3,4]').Vec(); v
            [[1, 3]~, [2, 4]~]
            sage: v.Mat()
            [1, 2; 3, 4]
            sage: v = pari('[1,2;3,4]').Col(); v
            [[1, 2], [3, 4]]~
            sage: v.Mat()
            [1, 2; 3, 4]
        """
        pari_catch_sig_on()
        return P.new_gen(gtomat(x.g))

    def Mod(gen x, y):
        """
        Mod(x, y): Returns the object x modulo y, denoted Mod(x, y).

        The input y must be a an integer or a polynomial:

        - If y is an INTEGER, x must also be an integer, a rational
          number, or a p-adic number compatible with the modulus y.

        - If y is a POLYNOMIAL, x must be a scalar (which is not a
          polmod), a polynomial, a rational function, or a power
          series.

        .. warning::

           This function is not the same as ``x % y`` which is an
           integer or a polynomial.

        INPUT:


        -  ``x`` - gen

        -  ``y`` - integer or polynomial


        OUTPUT:


        -  ``gen`` - intmod or polmod


        EXAMPLES::

            sage: z = pari(3)
            sage: x = z.Mod(pari(7))
            sage: x
            Mod(3, 7)
            sage: x^2
            Mod(2, 7)
            sage: x^100
            Mod(4, 7)
            sage: x.type()
            't_INTMOD'

        ::

            sage: f = pari("x^2 + x + 1")
            sage: g = pari("x")
            sage: a = g.Mod(f)
            sage: a
            Mod(x, x^2 + x + 1)
            sage: a*a
            Mod(-x - 1, x^2 + x + 1)
            sage: a.type()
            't_POLMOD'
        """
        t0GEN(y)
        pari_catch_sig_on()
        return P.new_gen(gmodulo(x.g,t0))

    def Pol(self, v=-1):
        """
        Pol(x, v): convert x into a polynomial with main variable v and
        return the result.

        - If x is a scalar, returns a constant polynomial.

        - If x is a power series, the effect is identical to
          ``truncate``, i.e. it chops off the `O(X^k)`.

        - If x is a vector, this function creates the polynomial whose
          coefficients are given in x, with x[0] being the leading
          coefficient (which can be zero).

        .. warning::

           This is *not* a substitution function. It will not
           transform an object containing variables of higher priority
           than v::

               sage: pari('x+y').Pol('y')
               Traceback (most recent call last):
               ...
               PariError:  (5)

        INPUT:


        -  ``x`` - gen

        -  ``v`` - (optional) which variable, defaults to 'x'


        OUTPUT:


        -  ``gen`` - a polynomial


        EXAMPLES::

            sage: v = pari("[1,2,3,4]")
            sage: f = v.Pol()
            sage: f
            x^3 + 2*x^2 + 3*x + 4
            sage: f*f
            x^6 + 4*x^5 + 10*x^4 + 20*x^3 + 25*x^2 + 24*x + 16

        ::

            sage: v = pari("[1,2;3,4]")
            sage: v.Pol()
            [1, 3]~*x + [2, 4]~
        """
        pari_catch_sig_on()
        return P.new_gen(gtopoly(self.g, P.get_var(v)))

    def Polrev(self, v=-1):
        """
        Polrev(x, v): Convert x into a polynomial with main variable v and
        return the result. This is the reverse of Pol if x is a vector,
        otherwise it is identical to Pol. By "reverse" we mean that the
        coefficients are reversed.

        INPUT:

        -  ``x`` - gen

        OUTPUT:

        -  ``gen`` - a polynomial

        EXAMPLES::

            sage: v = pari("[1,2,3,4]")
            sage: f = v.Polrev()
            sage: f
            4*x^3 + 3*x^2 + 2*x + 1
            sage: v.Pol()
            x^3 + 2*x^2 + 3*x + 4
            sage: v.Polrev('y')
            4*y^3 + 3*y^2 + 2*y + 1

        Note that Polrev does *not* reverse the coefficients of a
        polynomial! ::

            sage: f
            4*x^3 + 3*x^2 + 2*x + 1
            sage: f.Polrev()
            4*x^3 + 3*x^2 + 2*x + 1
            sage: v = pari("[1,2;3,4]")
            sage: v.Polrev()
            [2, 4]~*x + [1, 3]~
        """
        pari_catch_sig_on()
        return P.new_gen(gtopolyrev(self.g, P.get_var(v)))

    def Qfb(gen a, b, c, D=0):
        """
        Qfb(a,b,c,D=0.): Returns the binary quadratic form

        .. math::

                                ax^2 + bxy + cy^2.


        The optional D is 0 by default and initializes Shank's distance if
        `b^2 - 4ac > 0`.  The discriminant of the quadratic form must not
        be a perfect square.

        .. note::

           Negative definite forms are not implemented, so use their
           positive definite counterparts instead. (I.e., if f is a
           negative definite quadratic form, then -f is positive
           definite.)

        INPUT:


        -  ``a`` - gen

        -  ``b`` - gen

        -  ``c`` - gen

        -  ``D`` - gen (optional, defaults to 0)


        OUTPUT:


        -  ``gen`` - binary quadratic form


        EXAMPLES::

            sage: pari(3).Qfb(7, 1)
            Qfb(3, 7, 1, 0.E-19)
            sage: pari(3).Qfb(7, 2)  # discriminant is 25
            Traceback (most recent call last):
            ...
            PariError:  (5)
        """
        t0GEN(b); t1GEN(c); t2GEN(D)
        pari_catch_sig_on()
        return P.new_gen(Qfb0(a.g, t0, t1, t2, prec))


    def Ser(gen x, v=-1, long seriesprecision = 16):
        """
        Ser(x,v=x): Create a power series from x with main variable v and
        return the result.

        - If x is a scalar, this gives a constant power series with
          precision given by the default series precision, as returned
          by get_series_precision().

        - If x is a polynomial, the precision is the greatest of
          get_series_precision() and the degree of the polynomial.

        - If x is a vector, the precision is similarly given, and the
          coefficients of the vector are understood to be the
          coefficients of the power series starting from the constant
          term (i.e. the reverse of the function Pol).

        .. warning::

           This is *not* a substitution function. It will not
           transform an object containing variables of higher priority
           than v.

        INPUT:


        -  ``x`` - gen

        -  ``v`` - PARI variable (default: x)


        OUTPUT:


        -  ``gen`` - PARI object of PARI type t_SER


        EXAMPLES::

            sage: pari(2).Ser()
            2 + O(x^16)
            sage: x = pari([1,2,3,4,5])
            sage: x.Ser()
            1 + 2*x + 3*x^2 + 4*x^3 + 5*x^4 + O(x^5)
            sage: f = x.Ser('v'); print f
            1 + 2*v + 3*v^2 + 4*v^3 + 5*v^4 + O(v^5)
            sage: pari(1)/f
            1 - 2*v + v^2 + O(v^5)
            sage: pari('x^5').Ser(seriesprecision = 20)
            x^5 + O(x^25)
            sage: pari('1/x').Ser(seriesprecision = 1)
            x^-1 + O(x^0)
        """
        pari_catch_sig_on()
        return P.new_gen(gtoser(x.g, P.get_var(v), seriesprecision))


    def Set(gen x):
        """
        Set(x): convert x into a set, i.e. a row vector of strings in
        increasing lexicographic order.

        INPUT:


        -  ``x`` - gen


        OUTPUT:


        -  ``gen`` - a vector of strings in increasing
           lexicographic order.


        EXAMPLES::

            sage: pari([1,5,2]).Set()
            ["1", "2", "5"]
            sage: pari([]).Set()     # the empty set
            []
            sage: pari([1,1,-1,-1,3,3]).Set()
            ["-1", "1", "3"]
            sage: pari(1).Set()
            ["1"]
            sage: pari('1/(x*y)').Set()
            ["1/(y*x)"]
            sage: pari('["bc","ab","bc"]').Set()
            ["\"ab\"", "\"bc\""]
        """
        pari_catch_sig_on()
        return P.new_gen(gtoset(x.g))


    def Str(self):
        """
        Str(self): Return the print representation of self as a PARI
        object.

        INPUT:


        -  ``self`` - gen


        OUTPUT:


        -  ``gen`` - a PARI gen of type t_STR, i.e., a PARI
           string


        EXAMPLES::

            sage: pari([1,2,['abc',1]]).Str()
            "[1, 2, [abc, 1]]"
            sage: pari([1,1, 1.54]).Str()
            "[1, 1, 1.54000000000000]"
            sage: pari(1).Str()       # 1 is automatically converted to string rep
            "1"
            sage: x = pari('x')       # PARI variable "x"
            sage: x.Str()             # is converted to string rep.
            "x"
            sage: x.Str().type()
            't_STR'
        """
        cdef char* c
        pari_catch_sig_on()
        c = GENtostr(self.g)
        v = self.new_gen(strtoGENstr(c))
        pari_free(c)
        return v


    def Strchr(gen x):
        """
        Strchr(x): converts x to a string, translating each integer into a
        character (in ASCII).

        .. note::

           :meth:`.Vecsmall` is (essentially) the inverse to :meth:`.Strchr`.

        INPUT:


        -  ``x`` - PARI vector of integers


        OUTPUT:


        -  ``gen`` - a PARI string


        EXAMPLES::

            sage: pari([65,66,123]).Strchr()
            "AB{"
            sage: pari('"Sage"').Vecsmall()   # pari('"Sage"') --> PARI t_STR
            Vecsmall([83, 97, 103, 101])
            sage: _.Strchr()
            "Sage"
            sage: pari([83, 97, 103, 101]).Strchr()
            "Sage"
        """
        pari_catch_sig_on()
        return P.new_gen(Strchr(x.g))

    def Strexpand(gen x):
        """
        Strexpand(x): Concatenate the entries of the vector x into a single
        string, performing tilde expansion.

        .. note::

           I have no clue what the point of this function is. - William
        """
        if typ(x.g) != t_VEC:
            raise TypeError, "x must be of type t_VEC."
        pari_catch_sig_on()
        return P.new_gen(Strexpand(x.g))


    def Strtex(gen x):
        r"""
        Strtex(x): Translates the vector x of PARI gens to TeX format and
        returns the resulting concatenated strings as a PARI t_STR.

        INPUT:


        -  ``x`` - gen


        OUTPUT:


        -  ``gen`` - PARI t_STR (string)


        EXAMPLES::

            sage: v=pari('x^2')
            sage: v.Strtex()
            "x^2"
            sage: v=pari(['1/x^2','x'])
            sage: v.Strtex()
            "\\frac{1}{x^2}x"
            sage: v=pari(['1 + 1/x + 1/(y+1)','x-1'])
            sage: v.Strtex()
            "\\frac{ \\left(y\n + 2\\right)  x\n + \\left(y\n + 1\\right) }{ \\left(y\n + 1\\right)  x}x\n - 1"
        """
        if typ(x.g) != t_VEC:
            x = P.vector(1, [x])
        pari_catch_sig_on()
        return P.new_gen(Strtex(x.g))

    def printtex(gen x):
        return x.Strtex()

    def Vec(gen x, long n = 0):
        """
        Transform the object `x` into a vector with minimal size `|n|`.

        INPUT:

        - ``x`` -- gen

        - ``n`` -- Make the vector of minimal length `|n|`. If `n > 0`,
          append zeros; if `n < 0`, prepend zeros.

        OUTPUT:

        A PARI vector (type ``t_VEC``)

        EXAMPLES::

            sage: pari(1).Vec()
            [1]
            sage: pari('x^3').Vec()
            [1, 0, 0, 0]
            sage: pari('x^3 + 3*x - 2').Vec()
            [1, 0, 3, -2]
            sage: pari([1,2,3]).Vec()
            [1, 2, 3]
            sage: pari('[1, 2; 3, 4]').Vec()
            [[1, 3]~, [2, 4]~]
            sage: pari('"Sage"').Vec()
            ["S", "a", "g", "e"]
            sage: pari('2*x^2 + 3*x^3 + O(x^5)').Vec()
            [2, 3, 0]
            sage: pari('2*x^-2 + 3*x^3 + O(x^5)').Vec()
            [2, 0, 0, 0, 0, 3, 0]

        Note the different term ordering for polynomials and series::

            sage: pari('1 + x + 3*x^3 + O(x^5)').Vec()
            [1, 1, 0, 3, 0]
            sage: pari('1 + x + 3*x^3').Vec()
            [3, 0, 1, 1]

        We demonstate the `n` argument::

            sage: pari([1,2,3,4]).Vec(2)
            [1, 2, 3, 4]
            sage: pari([1,2,3,4]).Vec(-2)
            [1, 2, 3, 4]
            sage: pari([1,2,3,4]).Vec(6)
            [1, 2, 3, 4, 0, 0]
            sage: pari([1,2,3,4]).Vec(-6)
            [0, 0, 1, 2, 3, 4]

        See also :meth:`Col` (create a column vector) and :meth:`Vecrev`
        (create a vector in reversed order).
        """
        pari_catch_sig_on()
        return P.new_gen(_Vec_append(gtovec(x.g), gen_0, n))

    def Vecrev(gen x, long n = 0):
        """
        Transform the object `x` into a vector with minimal size `|n|`.
        The order of the resulting vector is reversed compared to :meth:`Vec`.

        INPUT:

        - ``x`` -- gen

        - ``n`` -- Make the vector of minimal length `|n|`. If `n > 0`,
          prepend zeros; if `n < 0`, append zeros.

        OUTPUT:

        A PARI vector (type ``t_VEC``)

        EXAMPLES::

            sage: pari(1).Vecrev()
            [1]
            sage: pari('x^3').Vecrev()
            [0, 0, 0, 1]
            sage: pari('x^3 + 3*x - 2').Vecrev()
            [-2, 3, 0, 1]
            sage: pari([1, 2, 3]).Vecrev()
            [3, 2, 1]
            sage: pari('Col([1, 2, 3])').Vecrev()
            [3, 2, 1]
            sage: pari('[1, 2; 3, 4]').Vecrev()
            [[2, 4]~, [1, 3]~]
            sage: pari('"Sage"').Vecrev()
            ["e", "g", "a", "S"]

        We demonstate the `n` argument::

            sage: pari([1,2,3,4]).Vecrev(2)
            [4, 3, 2, 1]
            sage: pari([1,2,3,4]).Vecrev(-2)
            [4, 3, 2, 1]
            sage: pari([1,2,3,4]).Vecrev(6)
            [0, 0, 4, 3, 2, 1]
            sage: pari([1,2,3,4]).Vecrev(-6)
            [4, 3, 2, 1, 0, 0]
        """
        pari_catch_sig_on()
        return P.new_gen(_Vec_append(gtovecrev(x.g), gen_0, -n))

    def Vecsmall(gen x, long n = 0):
        """
        Transform the object `x` into a ``t_VECSMALL`` with minimal size `|n|`.

        INPUT:

        - ``x`` -- gen

        - ``n`` -- Make the vector of minimal length `|n|`. If `n > 0`,
          append zeros; if `n < 0`, prepend zeros.

        OUTPUT:

        A PARI vector of small integers (type ``t_VECSMALL``)

        EXAMPLES::

            sage: pari([1,2,3]).Vecsmall()
            Vecsmall([1, 2, 3])
            sage: pari('"Sage"').Vecsmall()
            Vecsmall([83, 97, 103, 101])
            sage: pari(1234).Vecsmall()
            Vecsmall([1234])
            sage: pari('x^2 + 2*x + 3').Vecsmall()
            Traceback (most recent call last):
            ...
            PariError: incorrect type (11)

        We demonstate the `n` argument::

            sage: pari([1,2,3]).Vecsmall(2)
            Vecsmall([1, 2, 3])
            sage: pari([1,2,3]).Vecsmall(-2)
            Vecsmall([1, 2, 3])
            sage: pari([1,2,3]).Vecsmall(6)
            Vecsmall([1, 2, 3, 0, 0, 0])
            sage: pari([1,2,3]).Vecsmall(-6)
            Vecsmall([0, 0, 0, 1, 2, 3])
        """
        pari_catch_sig_on()
        return P.new_gen(_Vec_append(gtovecsmall(x.g), <GEN>0, n))

    def binary(gen x):
        """
        binary(x): gives the vector formed by the binary digits of abs(x),
        where x is of type t_INT.

        INPUT:


        -  ``x`` - gen of type t_INT


        OUTPUT:


        -  ``gen`` - of type t_VEC


        EXAMPLES::

            sage: pari(0).binary()
            [0]
            sage: pari(-5).binary()
            [1, 0, 1]
            sage: pari(5).binary()
            [1, 0, 1]
            sage: pari(2005).binary()
            [1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1]

        ::

            sage: pari('"2"').binary()
            Traceback (most recent call last):
            ...
            TypeError: x (="2") must be of type t_INT, but is of type t_STR.
        """
        if typ(x.g) != t_INT:
            raise TypeError, "x (=%s) must be of type t_INT, but is of type %s."%(x,x.type())
        pari_catch_sig_on()
        return P.new_gen(binaire(x.g))

    def bitand(gen x, y):
        """
        bitand(x,y): Bitwise and of two integers x and y. Negative numbers
        behave as if modulo some large power of 2.

        INPUT:


        -  ``x`` - gen (of type t_INT)

        -  ``y`` - coercible to gen (of type t_INT)


        OUTPUT:


        -  ``gen`` - of type type t_INT


        EXAMPLES::

            sage: pari(8).bitand(4)
            0
            sage: pari(8).bitand(8)
            8
            sage: pari(6).binary()
            [1, 1, 0]
            sage: pari(7).binary()
            [1, 1, 1]
            sage: pari(6).bitand(7)
            6
            sage: pari(19).bitand(-1)
            19
            sage: pari(-1).bitand(-1)
            -1
        """
        t0GEN(y)
        pari_catch_sig_on()
        return P.new_gen(gbitand(x.g, t0))


    def bitneg(gen x, long n=-1):
        r"""
        bitneg(x,n=-1): Bitwise negation of the integer x truncated to n
        bits. n=-1 (the default) represents an infinite sequence of the bit
        1. Negative numbers behave as if modulo some large power of 2.

        With n=-1, this function returns -n-1. With n = 0, it returns a
        number a such that `a\cong -n-1 \pmod{2^n}`.

        INPUT:


        -  ``x`` - gen (t_INT)

        -  ``n`` - long, default = -1


        OUTPUT:


        -  ``gen`` - t_INT


        EXAMPLES::

            sage: pari(10).bitneg()
            -11
            sage: pari(1).bitneg()
            -2
            sage: pari(-2).bitneg()
            1
            sage: pari(-1).bitneg()
            0
            sage: pari(569).bitneg()
            -570
            sage: pari(569).bitneg(10)
            454
            sage: 454 % 2^10
            454
            sage: -570 % 2^10
            454
        """
        pari_catch_sig_on()
        return P.new_gen(gbitneg(x.g,n))


    def bitnegimply(gen x, y):
        """
        bitnegimply(x,y): Bitwise negated imply of two integers x and y, in
        other words, x BITAND BITNEG(y). Negative numbers behave as if
        modulo big power of 2.

        INPUT:


        -  ``x`` - gen (of type t_INT)

        -  ``y`` - coercible to gen (of type t_INT)


        OUTPUT:


        -  ``gen`` - of type type t_INT


        EXAMPLES::

            sage: pari(14).bitnegimply(0)
            14
            sage: pari(8).bitnegimply(8)
            0
            sage: pari(8+4).bitnegimply(8)
            4
        """
        t0GEN(y)
        pari_catch_sig_on()
        return P.new_gen(gbitnegimply(x.g, t0))


    def bitor(gen x, y):
        """
        bitor(x,y): Bitwise or of two integers x and y. Negative numbers
        behave as if modulo big power of 2.

        INPUT:


        -  ``x`` - gen (of type t_INT)

        -  ``y`` - coercible to gen (of type t_INT)


        OUTPUT:


        -  ``gen`` - of type type t_INT


        EXAMPLES::

            sage: pari(14).bitor(0)
            14
            sage: pari(8).bitor(4)
            12
            sage: pari(12).bitor(1)
            13
            sage: pari(13).bitor(1)
            13
        """
        t0GEN(y)
        pari_catch_sig_on()
        return P.new_gen(gbitor(x.g, t0))


    def bittest(gen x, long n):
        """
        bittest(x, long n): Returns bit number n (coefficient of
        `2^n` in binary) of the integer x. Negative numbers behave
        as if modulo a big power of 2.

        INPUT:


        -  ``x`` - gen (pari integer)


        OUTPUT:


        -  ``bool`` - a Python bool


        EXAMPLES::

            sage: x = pari(6)
            sage: x.bittest(0)
            False
            sage: x.bittest(1)
            True
            sage: x.bittest(2)
            True
            sage: x.bittest(3)
            False
            sage: pari(-3).bittest(0)
            True
            sage: pari(-3).bittest(1)
            False
            sage: [pari(-3).bittest(n) for n in range(10)]
            [True, False, True, True, True, True, True, True, True, True]
        """
        pari_catch_sig_on()
        b = bool(bittest(x.g, n))
        pari_catch_sig_off()
        return b

    def bitxor(gen x, y):
        """
        bitxor(x,y): Bitwise exclusive or of two integers x and y. Negative
        numbers behave as if modulo big power of 2.

        INPUT:


        -  ``x`` - gen (of type t_INT)

        -  ``y`` - coercible to gen (of type t_INT)


        OUTPUT:


        -  ``gen`` - of type type t_INT


        EXAMPLES::

            sage: pari(6).bitxor(4)
            2
            sage: pari(0).bitxor(4)
            4
            sage: pari(6).bitxor(0)
            6
        """
        t0GEN(y)
        pari_catch_sig_on()
        return P.new_gen(gbitxor(x.g, t0))


    def ceil(gen x):
        """
        For real x: return the smallest integer = x. For rational
        functions: the quotient of numerator by denominator. For lists:
        apply componentwise.

        INPUT:


        -  ``x`` - gen


        OUTPUT:


        -  ``gen`` - depends on type of x


        EXAMPLES::

            sage: pari(1.4).ceil()
            2
            sage: pari(-1.4).ceil()
            -1
            sage: pari(3/4).ceil()
            1
            sage: pari(x).ceil()
            x
            sage: pari((x^2+x+1)/x).ceil()
            x + 1

        This may be unexpected: but it is correct, treating the argument as
        a rational function in RR(x).

        ::

            sage: pari(x^2+5*x+2.5).ceil()
            x^2 + 5*x + 2.50000000000000
        """
        pari_catch_sig_on()
        return P.new_gen(gceil(x.g))

    def centerlift(gen x, v=-1):
        """
        centerlift(x,v): Centered lift of x. This function returns exactly
        the same thing as lift, except if x is an integer mod.

        INPUT:


        -  ``x`` - gen

        -  ``v`` - var (default: x)


        OUTPUT: gen

        EXAMPLES::

            sage: x = pari(-2).Mod(5)
            sage: x.centerlift()
            -2
            sage: x.lift()
            3
            sage: f = pari('x-1').Mod('x^2 + 1')
            sage: f.centerlift()
            x - 1
            sage: f.lift()
            x - 1
            sage: f = pari('x-y').Mod('x^2+1')
            sage: f
            Mod(x - y, x^2 + 1)
            sage: f.centerlift('x')
            x - y
            sage: f.centerlift('y')
            Mod(x - y, x^2 + 1)
        """
        pari_catch_sig_on()
        return P.new_gen(centerlift0(x.g, P.get_var(v)))


    def component(gen x, long n):
        """
        component(x, long n): Return n'th component of the internal
        representation of x. This function is 1-based instead of 0-based.

        .. note::

           For vectors or matrices, it is simpler to use x[n-1]. For
           list objects such as is output by nfinit, it is easier to
           use member functions.

        INPUT:


        -  ``x`` - gen

        -  ``n`` - C long (coercible to)


        OUTPUT: gen

        EXAMPLES::

            sage: pari([0,1,2,3,4]).component(1)
            0
            sage: pari([0,1,2,3,4]).component(2)
            1
            sage: pari([0,1,2,3,4]).component(4)
            3
            sage: pari('x^3 + 2').component(1)
            2
            sage: pari('x^3 + 2').component(2)
            0
            sage: pari('x^3 + 2').component(4)
            1

        ::

            sage: pari('x').component(0)
            Traceback (most recent call last):
            ...
            PariError:  (5)
        """
        pari_catch_sig_on()
        return P.new_gen(compo(x.g, n))

    def conj(gen x):
        """
        conj(x): Return the algebraic conjugate of x.

        INPUT:


        -  ``x`` - gen


        OUTPUT: gen

        EXAMPLES::

            sage: pari('x+1').conj()
            x + 1
            sage: pari('x+I').conj()
            x - I
            sage: pari('1/(2*x+3*I)').conj()
            1/(2*x - 3*I)
            sage: pari([1,2,'2-I','Mod(x,x^2+1)', 'Mod(x,x^2-2)']).conj()
            [1, 2, 2 + I, Mod(-x, x^2 + 1), Mod(-x, x^2 - 2)]
            sage: pari('Mod(x,x^2-2)').conj()
            Mod(-x, x^2 - 2)
            sage: pari('Mod(x,x^3-3)').conj()
            Traceback (most recent call last):
            ...
            PariError: incorrect type (11)
        """
        pari_catch_sig_on()
        return P.new_gen(gconj(x.g))

    def conjvec(gen x):
        """
        conjvec(x): Returns the vector of all conjugates of the algebraic
        number x. An algebraic number is a polynomial over Q modulo an
        irreducible polynomial.

        INPUT:


        -  ``x`` - gen


        OUTPUT: gen

        EXAMPLES::

            sage: pari('Mod(1+x,x^2-2)').conjvec()
            [-0.414213562373095, 2.41421356237310]~
            sage: pari('Mod(x,x^3-3)').conjvec()
            [1.44224957030741, -0.721124785153704 + 1.24902476648341*I, -0.721124785153704 - 1.24902476648341*I]~
        """
        pari_catch_sig_on()
        return P.new_gen(conjvec(x.g, prec))

    def denominator(gen x):
        """
        denominator(x): Return the denominator of x. When x is a vector,
        this is the least common multiple of the denominators of the
        components of x.

        what about poly? INPUT:


        -  ``x`` - gen


        OUTPUT: gen

        EXAMPLES::

            sage: pari('5/9').denominator()
            9
            sage: pari('(x+1)/(x-2)').denominator()
            x - 2
            sage: pari('2/3 + 5/8*x + 7/3*x^2 + 1/5*y').denominator()
            1
            sage: pari('2/3*x').denominator()
            1
            sage: pari('[2/3, 5/8, 7/3, 1/5]').denominator()
            120
        """
        pari_catch_sig_on()
        return P.new_gen(denom(x.g))

    def floor(gen x):
        """
        For real x: return the largest integer = x. For rational functions:
        the quotient of numerator by denominator. For lists: apply
        componentwise.

        INPUT:


        -  ``x`` - gen


        OUTPUT: gen

        EXAMPLES::

            sage: pari(5/9).floor()
            0
            sage: pari(11/9).floor()
            1
            sage: pari(1.17).floor()
            1
            sage: pari([1.5,2.3,4.99]).floor()
            [1, 2, 4]
            sage: pari([[1.1,2.2],[3.3,4.4]]).floor()
            [[1, 2], [3, 4]]
            sage: pari(x).floor()
            x
            sage: pari((x^2+x+1)/x).floor()
            x + 1
            sage: pari(x^2+5*x+2.5).floor()
            x^2 + 5*x + 2.50000000000000

        ::

            sage: pari('"hello world"').floor()
            Traceback (most recent call last):
            ...
            PariError: incorrect type (11)
        """
        pari_catch_sig_on()
        return P.new_gen(gfloor(x.g))

    def frac(gen x):
        """
        frac(x): Return the fractional part of x, which is x - floor(x).

        INPUT:


        -  ``x`` - gen


        OUTPUT: gen

        EXAMPLES::

            sage: pari(1.75).frac()
            0.750000000000000
            sage: pari(sqrt(2)).frac()
            0.414213562373095
            sage: pari('sqrt(-2)').frac()
            Traceback (most recent call last):
            ...
            PariError: incorrect type (11)
        """
        pari_catch_sig_on()
        return P.new_gen(gfrac(x.g))

    def imag(gen x):
        """
        imag(x): Return the imaginary part of x. This function also works
        component-wise.

        INPUT:


        -  ``x`` - gen


        OUTPUT: gen

        EXAMPLES::

            sage: pari('1+2*I').imag()
            2
            sage: pari(sqrt(-2)).imag()
            1.41421356237310
            sage: pari('x+I').imag()
            1
            sage: pari('x+2*I').imag()
            2
            sage: pari('(1+I)*x^2+2*I').imag()
            x^2 + 2
            sage: pari('[1,2,3] + [4*I,5,6]').imag()
            [4, 0, 0]
        """
        pari_catch_sig_on()
        return P.new_gen(gimag(x.g))

    def length(self):
        """

        """
        return glength(self.g)

    def lift(gen x, v=-1):
        """
        lift(x,v): Returns the lift of an element of Z/nZ to Z or R[x]/(P)
        to R[x] for a type R if v is omitted. If v is given, lift only
        polymods with main variable v. If v does not occur in x, lift only
        intmods.

        INPUT:


        -  ``x`` - gen

        -  ``v`` - (optional) variable


        OUTPUT: gen

        EXAMPLES::

            sage: x = pari("x")
            sage: a = x.Mod('x^3 + 17*x + 3')
            sage: a
            Mod(x, x^3 + 17*x + 3)
            sage: b = a^4; b
            Mod(-17*x^2 - 3*x, x^3 + 17*x + 3)
            sage: b.lift()
            -17*x^2 - 3*x

        ??? more examples
        """
        pari_catch_sig_on()
        if v == -1:
            return P.new_gen(lift(x.g))
        return P.new_gen(lift0(x.g, P.get_var(v)))

    def numbpart(gen x):
        """
        numbpart(x): returns the number of partitions of x.

        EXAMPLES::

            sage: pari(20).numbpart()
            627
            sage: pari(100).numbpart()
            190569292
        """
        pari_catch_sig_on()
        return P.new_gen(numbpart(x.g))

    def numerator(gen x):
        """
        numerator(x): Returns the numerator of x.

        INPUT:


        -  ``x`` - gen


        OUTPUT: gen

        EXAMPLES:
        """
        pari_catch_sig_on()
        return P.new_gen(numer(x.g))


    def numtoperm(gen k, long n):
        """
        numtoperm(k, n): Return the permutation number k (mod n!) of n
        letters, where n is an integer.

        INPUT:


        -  ``k`` - gen, integer

        -  ``n`` - int


        OUTPUT:


        -  ``gen`` - vector (permutation of 1,...,n)


        EXAMPLES:
        """
        pari_catch_sig_on()
        return P.new_gen(numtoperm(n, k.g))


    def padicprec(gen x, p):
        """
        padicprec(x,p): Return the absolute p-adic precision of the object
        x.

        INPUT:


        -  ``x`` - gen


        OUTPUT: int

        EXAMPLES::

            sage: K = Qp(11,5)
            sage: x = K(11^-10 + 5*11^-7 + 11^-6)
            sage: y = pari(x)
            sage: y.padicprec(11)
            -5
            sage: y.padicprec(17)
            Traceback (most recent call last):
            ...
            ValueError: not the same prime in padicprec
        """
        cdef gen _p
        _p = pari(p)
        if typ(_p.g) != t_INT:
            raise TypeError("p (=%s) must be of type t_INT, but is of type %s."%(
                _p, _p.type()))
        if not gequal(gel(x.g, 2), _p.g):
            raise ValueError("not the same prime in padicprec")
        return padicprec(x.g, _p.g)

    def padicprime(gen x):
        """
        The uniformizer of the p-adic ring this element lies in, as a t_INT.

        INPUT:

        - ``x`` - gen, of type t_PADIC

        OUTPUT:

        - ``p`` - gen, of type t_INT

        EXAMPLES::

            sage: K = Qp(11,5)
            sage: x = K(11^-10 + 5*11^-7 + 11^-6)
            sage: y = pari(x)
            sage: y.padicprime()
            11
            sage: y.padicprime().type()
            't_INT'
        """
        pari_catch_sig_on()
        return P.new_gen(gel(x.g, 2))

    def permtonum(gen x):
        """
        permtonum(x): Return the ordinal (between 1 and n!) of permutation
        vector x. ??? Huh ??? say more. what is a perm vector. 0 to n-1 or
        1-n.

        INPUT:


        -  ``x`` - gen (vector of integers)


        OUTPUT:


        -  ``gen`` - integer


        EXAMPLES:
        """
        if typ(x.g) != t_VEC:
            raise TypeError, "x (=%s) must be of type t_VEC, but is of type %s."%(x,x.type())
        pari_catch_sig_on()
        return P.new_gen(permtonum(x.g))

    def precision(gen x, long n=-1):
        """
        precision(x,n): Change the precision of x to be n, where n is a
        C-integer). If n is omitted, output the real precision of x.

        INPUT:


        -  ``x`` - gen

        -  ``n`` - (optional) int


        OUTPUT: nothing or gen if n is omitted

        EXAMPLES:
        """
        if n <= -1:
            return precision(x.g)
        pari_catch_sig_on()
        return P.new_gen(precision0(x.g, n))

    def random(gen N):
        r"""
        ``random(N=2^31)``: Return a pseudo-random integer
        between 0 and `N-1`.

        INPUT:


        -``N`` - gen, integer


        OUTPUT:


        -  ``gen`` - integer


        EXAMPLES:
        """
        if typ(N.g) != t_INT:
            raise TypeError, "x (=%s) must be of type t_INT, but is of type %s."%(N,N.type())
        pari_catch_sig_on()
        return P.new_gen(genrand(N.g))

    def real(gen x):
        """
        real(x): Return the real part of x.

        INPUT:


        -  ``x`` - gen


        OUTPUT: gen

        EXAMPLES:
        """
        pari_catch_sig_on()
        return P.new_gen(greal(x.g))

    def round(gen x, estimate=False):
        """
        round(x,estimate=False): If x is a real number, returns x rounded
        to the nearest integer (rounding up). If the optional argument
        estimate is True, also returns the binary exponent e of the
        difference between the original and the rounded value (the
        "fractional part") (this is the integer ceiling of log_2(error)).

        When x is a general PARI object, this function returns the result
        of rounding every coefficient at every level of PARI object. Note
        that this is different than what the truncate function does (see
        the example below).

        One use of round is to get exact results after a long approximate
        computation, when theory tells you that the coefficients must be
        integers.

        INPUT:


        -  ``x`` - gen

        -  ``estimate`` - (optional) bool, False by default


        OUTPUT:

        - if estimate is False, return a single gen.

        - if estimate is True, return rounded version of x and error
          estimate in bits, both as gens.

        EXAMPLES::

            sage: pari('1.5').round()
            2
            sage: pari('1.5').round(True)
            (2, -1)
            sage: pari('1.5 + 2.1*I').round()
            2 + 2*I
            sage: pari('1.0001').round(True)
            (1, -14)
            sage: pari('(2.4*x^2 - 1.7)/x').round()
            (2*x^2 - 2)/x
            sage: pari('(2.4*x^2 - 1.7)/x').truncate()
            2.40000000000000*x
        """
        cdef int n
        cdef long e
        cdef gen y
        pari_catch_sig_on()
        if not estimate:
            return P.new_gen(ground(x.g))
        y = P.new_gen(grndtoi(x.g, &e))
        return y, e

    def simplify(gen x):
        """
        simplify(x): Simplify the object x as much as possible, and return
        the result.

        A complex or quadratic number whose imaginary part is an exact 0
        (i.e., not an approximate one such as O(3) or 0.E-28) is converted
        to its real part, and a a polynomial of degree 0 is converted to
        its constant term. Simplification occurs recursively.

        This function is useful before using arithmetic functions, which
        expect integer arguments:

        EXAMPLES::

            sage: y = pari('y')
            sage: x = pari('9') + y - y
            sage: x
            9
            sage: x.type()
            't_POL'
            sage: x.factor()
            matrix(0,2)
            sage: pari('9').factor()
            Mat([3, 2])
            sage: x.simplify()
            9
            sage: x.simplify().factor()
            Mat([3, 2])
            sage: x = pari('1.5 + 0*I')
            sage: x.type()
            't_REAL'
            sage: x.simplify()
            1.50000000000000
            sage: y = x.simplify()
            sage: y.type()
            't_REAL'
        """
        pari_catch_sig_on()
        return P.new_gen(simplify(x.g))

    def sizeword(gen x):
        """
        Return the total number of machine words occupied by the
        complete tree of the object x.  A machine word is 32 or
        64 bits, depending on the computer.

        INPUT:

        -  ``x`` - gen

        OUTPUT: int (a Python int)

        EXAMPLES::

            sage: pari('0').sizeword()
            2
            sage: pari('1').sizeword()
            3
            sage: pari('1000000').sizeword()
            3
            sage: pari('10^100').sizeword()
            13      # 32-bit
            8       # 64-bit
            sage: pari(RDF(1.0)).sizeword()
            4       # 32-bit
            3       # 64-bit
            sage: pari('x').sizeword()
            9
            sage: pari('x^20').sizeword()
            66
            sage: pari('[x, I]').sizeword()
            20
        """
        return gsizeword(x.g)

    def sizebyte(gen x):
        """
        Return the total number of bytes occupied by the complete tree
        of the object x. Note that this number depends on whether the
        computer is 32-bit or 64-bit.

        INPUT:

        -  ``x`` - gen

        OUTPUT: int (a Python int)

        EXAMPLE::

            sage: pari('1').sizebyte()
            12           # 32-bit
            24           # 64-bit
        """
        return gsizebyte(x.g)

    def sizedigit(gen x):
        """
        sizedigit(x): Return a quick estimate for the maximal number of
        decimal digits before the decimal point of any component of x.

        INPUT:


        -  ``x`` - gen


        OUTPUT:


        -  ``int`` - Python integer


        EXAMPLES::

            sage: x = pari('10^100')
            sage: x.Str().length()
            101
            sage: x.sizedigit()
            101

        Note that digits after the decimal point are ignored.

        ::

            sage: x = pari('1.234')
            sage: x
            1.23400000000000
            sage: x.sizedigit()
            1

        The estimate can be one too big::

            sage: pari('7234.1').sizedigit()
            4
            sage: pari('9234.1').sizedigit()
            5
        """
        return sizedigit(x.g)

    def truncate(gen x, estimate=False):
        """
        truncate(x,estimate=False): Return the truncation of x. If estimate
        is True, also return the number of error bits.

        When x is in the real numbers, this means that the part after the
        decimal point is chopped away, e is the binary exponent of the
        difference between the original and truncated value (the
        "fractional part"). If x is a rational function, the result is the
        integer part (Euclidean quotient of numerator by denominator) and
        if requested the error estimate is 0.

        When truncate is applied to a power series (in X), it transforms it
        into a polynomial or a rational function with denominator a power
        of X, by chopping away the `O(X^k)`. Similarly, when
        applied to a p-adic number, it transforms it into an integer or a
        rational number by chopping away the `O(p^k)`.

        INPUT:


        -  ``x`` - gen

        -  ``estimate`` - (optional) bool, which is False by
           default


        OUTPUT:

        - if estimate is False, return a single gen.

        - if estimate is True, return rounded version of x and error
          estimate in bits, both as gens.

        EXAMPLES::

            sage: pari('(x^2+1)/x').round()
            (x^2 + 1)/x
            sage: pari('(x^2+1)/x').truncate()
            x
            sage: pari('1.043').truncate()
            1
            sage: pari('1.043').truncate(True)
            (1, -5)
            sage: pari('1.6').truncate()
            1
            sage: pari('1.6').round()
            2
            sage: pari('1/3 + 2 + 3^2 + O(3^3)').truncate()
            34/3
            sage: pari('sin(x+O(x^10))').truncate()
            1/362880*x^9 - 1/5040*x^7 + 1/120*x^5 - 1/6*x^3 + x
            sage: pari('sin(x+O(x^10))').round()   # each coefficient has abs < 1
            x + O(x^10)
        """
        cdef long e
        cdef gen y
        pari_catch_sig_on()
        if not estimate:
            return P.new_gen(gtrunc(x.g))
        y = P.new_gen(gcvtoi(x.g, &e))
        return y, e

    def valuation(gen x, p):
        """
        valuation(x,p): Return the valuation of x with respect to p.

        The valuation is the highest exponent of p dividing x.

        - If p is an integer, x must be an integer, an intmod whose
          modulus is divisible by p, a rational number, a p-adic
          number, or a polynomial or power series in which case the
          valuation is the minimum of the valuations of the
          coefficients.

        - If p is a polynomial, x must be a polynomial or a rational
          function. If p is a monomial then x may also be a power
          series.

        - If x is a vector, complex or quadratic number, then the
          valuation is the minimum of the component valuations.

        - If x = 0, the result is `2^31-1` on 32-bit machines or
          `2^63-1` on 64-bit machines if x is an exact
          object. If x is a p-adic number or power series, the result
          is the exponent of the zero.

        INPUT:


        -  ``x`` - gen

        -  ``p`` - coercible to gen


        OUTPUT:


        -  ``gen`` - integer


        EXAMPLES::

            sage: pari(9).valuation(3)
            2
            sage: pari(9).valuation(9)
            1
            sage: x = pari(9).Mod(27); x.valuation(3)
            2
            sage: pari('5/3').valuation(3)
            -1
            sage: pari('9 + 3*x + 15*x^2').valuation(3)
            1
            sage: pari([9,3,15]).valuation(3)
            1
            sage: pari('9 + 3*x + 15*x^2 + O(x^5)').valuation(3)
            1

        ::

            sage: pari('x^2*(x+1)^3').valuation(pari('x+1'))
            3
            sage: pari('x + O(x^5)').valuation('x')
            1
            sage: pari('2*x^2 + O(x^5)').valuation('x')
            2

        ::

            sage: pari(0).valuation(3)
            2147483647            # 32-bit
            9223372036854775807   # 64-bit
        """
        cdef long v
        t0GEN(p)
        pari_catch_sig_on()
        v = ggval(x.g, t0)
        pari_catch_sig_off()
        return v

    def _valp(gen x):
        """
        Return the valuation of x where x is a p-adic number (t_PADIC)
        or a laurent series (t_SER).  If x is a different type, this
        will give a bogus number.

        EXAMPLES::

            sage: pari('1/x^2 + O(x^10)')._valp()
            -2
            sage: pari('O(x^10)')._valp()
            10
            sage: pari('(1145234796 + O(3^10))/771966234')._valp()
            -2
            sage: pari('O(2^10)')._valp()
            10
            sage: pari('x')._valp()   # random
            -35184372088832
        """
        # This is a simple macro, so we don't need pari_catch_sig_on()
        return valp(x.g)

    def variable(gen x):
        """
        variable(x): Return the main variable of the object x, or p if x is
        a p-adic number.

        This function raises a TypeError exception on scalars, i.e., on
        objects with no variable associated to them.

        INPUT:


        -  ``x`` - gen


        OUTPUT: gen

        EXAMPLES::

            sage: pari('x^2 + x -2').variable()
            x
            sage: pari('1+2^3 + O(2^5)').variable()
            2
            sage: pari('x+y0').variable()
            x
            sage: pari('y0+z0').variable()
            y0
        """
        pari_catch_sig_on()
        return P.new_gen(gpolvar(x.g))


    ###########################################
    # 3: TRANSCENDENTAL functions
    # AUTHORS: Pyrex Code, docs -- Justin Walker (justin@mac.com)
    #          Examples, docs   -- William Stein
    ###########################################

    def abs(gen x):
        """
        Returns the absolute value of x (its modulus, if x is complex).
        Rational functions are not allowed. Contrary to most transcendental
        functions, an exact argument is not converted to a real number
        before applying abs and an exact result is returned if possible.

        EXAMPLES::

            sage: x = pari("-27.1")
            sage: x.abs()
            27.1000000000000

        If x is a polynomial, returns -x if the leading coefficient is real
        and negative else returns x. For a power series, the constant
        coefficient is considered instead.

        EXAMPLES::

            sage: pari('x-1.2*x^2').abs()
            1.20000000000000*x^2 - x
        """
        pari_catch_sig_on()
        # the prec parameter here has no effect
        return P.new_gen(gabs(x.g, prec))

    def acos(gen x, long precision=0):
        r"""
        The principal branch of `\cos^{-1}(x)`, so that
        `\RR e(\mathrm{acos}(x))` belongs to `[0,Pi]`. If `x`
        is real and `|x| > 1`, then `\mathrm{acos}(x)` is complex.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(0.5).acos()
            1.04719755119660
            sage: pari(1/2).acos()
            1.04719755119660
            sage: pari(1.1).acos()
            0.443568254385115*I
            sage: C.<i> = ComplexField()
            sage: pari(1.1+i).acos()
            0.849343054245252 - 1.09770986682533*I
        """
        pari_catch_sig_on()
        return P.new_gen(gacos(x.g, pbw(precision)))

    def acosh(gen x, precision=0):
        r"""
        The principal branch of `\cosh^{-1}(x)`, so that
        `\Im(\mathrm{acosh}(x))` belongs to `[0,Pi]`. If
        `x` is real and `x < 1`, then
        `\mathrm{acosh}(x)` is complex.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(2).acosh()
            1.31695789692482
            sage: pari(0).acosh()
            1.57079632679490*I
            sage: C.<i> = ComplexField()
            sage: pari(i).acosh()
            0.881373587019543 + 1.57079632679490*I
        """
        pari_catch_sig_on()
        return P.new_gen(gach(x.g, pbw(precision)))

    def agm(gen x, y, precision=0):
        r"""
        The arithmetic-geometric mean of x and y. In the case of complex or
        negative numbers, the principal square root is always chosen.
        p-adic or power series arguments are also allowed. Note that a
        p-adic AGM exists only if x/y is congruent to 1 modulo p (modulo 16
        for p=2). x and y cannot both be vectors or matrices.

        If any of `x` or `y` is an exact argument, it is
        first converted to a real or complex number using the optional
        parameter precision (in bits). If the arguments are inexact (e.g.
        real), the smallest of their two precisions is used in the
        computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(2).agm(2)
            2.00000000000000
            sage: pari(0).agm(1)
            0
            sage: pari(1).agm(2)
            1.45679103104691
            sage: C.<i> = ComplexField()
            sage: pari(1+i).agm(-3)
            -0.964731722290876 + 1.15700282952632*I
        """
        t0GEN(y)
        pari_catch_sig_on()
        return P.new_gen(agm(x.g, t0, pbw(precision)))

    def arg(gen x, precision=0):
        r"""
        arg(x): argument of x,such that `-\pi < \arg(x) \leq \pi`.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: C.<i> = ComplexField()
            sage: pari(2+i).arg()
            0.463647609000806
        """
        pari_catch_sig_on()
        return P.new_gen(garg(x.g, pbw(precision)))

    def asin(gen x, precision=0):
        r"""
        The principal branch of `\sin^{-1}(x)`, so that
        `\RR e(\mathrm{asin}(x))` belongs to `[-\pi/2,\pi/2]`. If
        `x` is real and `|x| > 1` then `\mathrm{asin}(x)`
        is complex.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(pari(0.5).sin()).asin()
            0.500000000000000
            sage: pari(2).asin()
            1.57079632679490 - 1.31695789692482*I
        """
        pari_catch_sig_on()
        return P.new_gen(gasin(x.g, pbw(precision)))

    def asinh(gen x, precision=0):
        r"""
        The principal branch of `\sinh^{-1}(x)`, so that
        `\Im(\mathrm{asinh}(x))` belongs to `[-\pi/2,\pi/2]`.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(2).asinh()
            1.44363547517881
            sage: C.<i> = ComplexField()
            sage: pari(2+i).asinh()
            1.52857091948100 + 0.427078586392476*I
        """
        pari_catch_sig_on()
        return P.new_gen(gash(x.g, pbw(precision)))

    def atan(gen x, precision=0):
        r"""
        The principal branch of `\tan^{-1}(x)`, so that
        `\RR e(\mathrm{atan}(x))` belongs to `]-\pi/2, \pi/2[`.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(1).atan()
            0.785398163397448
            sage: C.<i> = ComplexField()
            sage: pari(1.5+i).atan()
            1.10714871779409 + 0.255412811882995*I
        """
        pari_catch_sig_on()
        return P.new_gen(gatan(x.g, pbw(precision)))

    def atanh(gen x, precision=0):
        r"""
        The principal branch of `\tanh^{-1}(x)`, so that
        `\Im(\mathrm{atanh}(x))` belongs to `]-\pi/2,\pi/2]`. If
        `x` is real and `|x| > 1` then `\mathrm{atanh}(x)`
        is complex.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(0).atanh()
            0.E-19
            sage: pari(2).atanh()
            0.549306144334055 - 1.57079632679490*I
        """
        pari_catch_sig_on()
        return P.new_gen(gath(x.g, pbw(precision)))

    def bernfrac(gen x):
        r"""
        The Bernoulli number `B_x`, where `B_0 = 1`,
        `B_1 = -1/2`, `B_2 = 1/6,\ldots,` expressed as a
        rational number. The argument `x` should be of type
        integer.

        EXAMPLES::

            sage: pari(18).bernfrac()
            43867/798
            sage: [pari(n).bernfrac() for n in range(10)]
            [1, -1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0]
        """
        pari_catch_sig_on()
        return P.new_gen(bernfrac(x))

    def bernreal(gen x):
        r"""
        The Bernoulli number `B_x`, as for the function bernfrac,
        but `B_x` is returned as a real number (with the current
        precision).

        EXAMPLES::

            sage: pari(18).bernreal()
            54.9711779448622
        """
        pari_catch_sig_on()
        # the argument prec has no effect
        return P.new_gen(bernreal(x, prec))

    def bernvec(gen x):
        r"""
        Creates a vector containing, as rational numbers, the Bernoulli
        numbers `B_0, B_2,\ldots, B_{2x}`. This routine is
        obsolete. Use bernfrac instead each time you need a Bernoulli
        number in exact form.

        Note: this routine is implemented using repeated independent calls
        to bernfrac, which is faster than the standard recursion in exact
        arithmetic.

        EXAMPLES::

            sage: pari(8).bernvec()
            [1, 1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510]
            sage: [pari(2*n).bernfrac() for n in range(9)]
            [1, 1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510]
        """
        pari_catch_sig_on()
        return P.new_gen(bernvec(x))

    def besselh1(gen nu, x, precision=0):
        r"""
        The `H^1`-Bessel function of index `\nu` and
        argument `x`.

        If `nu` or `x` is an exact argument, it is first
        converted to a real or complex number using the optional parameter
        precision (in bits). If the arguments are inexact (e.g. real), the
        smallest of their precisions is used in the computation, and the
        parameter precision is ignored.

        EXAMPLES::

            sage: pari(2).besselh1(3)
            0.486091260585891 - 0.160400393484924*I
        """
        t0GEN(x)
        pari_catch_sig_on()
        return P.new_gen(hbessel1(nu.g, t0, pbw(precision)))

    def besselh2(gen nu, x, precision=0):
        r"""
        The `H^2`-Bessel function of index `\nu` and
        argument `x`.

        If `nu` or `x` is an exact argument, it is first
        converted to a real or complex number using the optional parameter
        precision (in bits). If the arguments are inexact (e.g. real), the
        smallest of their precisions is used in the computation, and the
        parameter precision is ignored.

        EXAMPLES::

            sage: pari(2).besselh2(3)
            0.486091260585891 + 0.160400393484924*I
        """
        t0GEN(x)
        pari_catch_sig_on()
        return P.new_gen(hbessel2(nu.g, t0, pbw(precision)))

    def besselj(gen nu, x, precision=0):
        r"""
        Bessel J function (Bessel function of the first kind), with index
        `\nu` and argument `x`. If `x` converts to
        a power series, the initial factor
        `(x/2)^{\nu}/\Gamma(\nu+1)` is omitted (since it cannot be
        represented in PARI when `\nu` is not integral).

        If `nu` or `x` is an exact argument, it is first
        converted to a real or complex number using the optional parameter
        precision (in bits). If the arguments are inexact (e.g. real), the
        smallest of their precisions is used in the computation, and the
        parameter precision is ignored.

        EXAMPLES::

            sage: pari(2).besselj(3)
            0.486091260585891
        """
        t0GEN(x)
        pari_catch_sig_on()
        return P.new_gen(jbessel(nu.g, t0, pbw(precision)))

    def besseljh(gen nu, x, precision=0):
        """
        J-Bessel function of half integral index (Spherical Bessel
        function of the first kind). More precisely, besseljh(n,x) computes
        `J_{n+1/2}(x)` where n must an integer, and x is any
        complex value. In the current implementation (PARI, version
        2.2.11), this function is not very accurate when `x` is
        small.

        If `nu` or `x` is an exact argument, it is first
        converted to a real or complex number using the optional parameter
        precision (in bits). If the arguments are inexact (e.g. real), the
        smallest of their precisions is used in the computation, and the
        parameter precision is ignored.

        EXAMPLES::

            sage: pari(2).besseljh(3)
            0.4127100324          # 32-bit
            0.412710032209716     # 64-bit
        """
        t0GEN(x)
        pari_catch_sig_on()
        return P.new_gen(jbesselh(nu.g, t0, pbw(precision)))

    def besseli(gen nu, x, precision=0):
        r"""
        Bessel I function (Bessel function of the second kind), with index
        `\nu` and argument `x`. If `x` converts to
        a power series, the initial factor
        `(x/2)^{\nu}/\Gamma(\nu+1)` is omitted (since it cannot be
        represented in PARI when `\nu` is not integral).

        If `nu` or `x` is an exact argument, it is first
        converted to a real or complex number using the optional parameter
        precision (in bits). If the arguments are inexact (e.g. real), the
        smallest of their precisions is used in the computation, and the
        parameter precision is ignored.

        EXAMPLES::

            sage: pari(2).besseli(3)
            2.24521244092995
            sage: C.<i> = ComplexField()
            sage: pari(2).besseli(3+i)
            1.12539407613913 + 2.08313822670661*I
        """
        t0GEN(x)
        pari_catch_sig_on()
        return P.new_gen(ibessel(nu.g, t0, pbw(precision)))

    def besselk(gen nu, x, long flag=0, precision=0):
        """
        nu.besselk(x, flag=0): K-Bessel function (modified Bessel function
        of the second kind) of index nu, which can be complex, and argument
        x.

        If `nu` or `x` is an exact argument, it is first
        converted to a real or complex number using the optional parameter
        precision (in bits). If the arguments are inexact (e.g. real), the
        smallest of their precisions is used in the computation, and the
        parameter precision is ignored.

        INPUT:


        -  ``nu`` - a complex number

        -  ``x`` - real number (positive or negative)

        -  ``flag`` - default: 0 or 1: use hyperu (hyperu is
           much slower for small x, and doesn't work for negative x).


        EXAMPLES::

            sage: C.<i> = ComplexField()
            sage: pari(2+i).besselk(3)
            0.0455907718407551 + 0.0289192946582081*I

        ::

            sage: pari(2+i).besselk(-3)
            -4.34870874986752 - 5.38744882697109*I

        ::

            sage: pari(2+i).besselk(300, flag=1)
            3.74224603319728 E-132 + 2.49071062641525 E-134*I
        """
        t0GEN(x)
        pari_catch_sig_on()
        return P.new_gen(kbessel(nu.g, t0, pbw(precision)))

    def besseln(gen nu, x, precision=0):
        """
        nu.besseln(x): Bessel N function (Spherical Bessel function of the
        second kind) of index nu and argument x.

        If `nu` or `x` is an exact argument, it is first
        converted to a real or complex number using the optional parameter
        precision (in bits). If the arguments are inexact (e.g. real), the
        smallest of their precisions is used in the computation, and the
        parameter precision is ignored.

        EXAMPLES::

            sage: C.<i> = ComplexField()
            sage: pari(2+i).besseln(3)
            -0.280775566958244 - 0.486708533223726*I
        """
        t0GEN(x)
        pari_catch_sig_on()
        return P.new_gen(nbessel(nu.g, t0, pbw(precision)))

    def cos(gen x, precision=0):
        """
        The cosine function.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(1.5).cos()
            0.0707372016677029
            sage: C.<i> = ComplexField()
            sage: pari(1+i).cos()
            0.833730025131149 - 0.988897705762865*I
            sage: pari('x+O(x^8)').cos()
            1 - 1/2*x^2 + 1/24*x^4 - 1/720*x^6 + 1/40320*x^8 + O(x^9)
        """
        pari_catch_sig_on()
        return P.new_gen(gcos(x.g, pbw(precision)))

    def cosh(gen x, precision=0):
        """
        The hyperbolic cosine function.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(1.5).cosh()
            2.35240961524325
            sage: C.<i> = ComplexField()
            sage: pari(1+i).cosh()
            0.833730025131149 + 0.988897705762865*I
            sage: pari('x+O(x^8)').cosh()
            1 + 1/2*x^2 + 1/24*x^4 + 1/720*x^6 + O(x^8)
        """
        pari_catch_sig_on()
        return P.new_gen(gch(x.g, pbw(precision)))

    def cotan(gen x, precision=0):
        """
        The cotangent of x.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(5).cotan()
            -0.295812915532746

        Computing the cotangent of `\pi` doesn't raise an error,
        but instead just returns a very large (positive or negative)
        number.

        ::

            sage: x = RR(pi)
            sage: pari(x).cotan()         # random
            -8.17674825 E15
        """
        pari_catch_sig_on()
        return P.new_gen(gcotan(x.g, pbw(precision)))

    def dilog(gen x, precision=0):
        r"""
        The principal branch of the dilogarithm of `x`, i.e. the
        analytic continuation of the power series
        `\log_2(x) = \sum_{n>=1} x^n/n^2`.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(1).dilog()
            1.64493406684823
            sage: C.<i> = ComplexField()
            sage: pari(1+i).dilog()
            0.616850275068085 + 1.46036211675312*I
        """
        pari_catch_sig_on()
        return P.new_gen(dilog(x.g, pbw(precision)))

    def eint1(gen x, long n=0, precision=0):
        r"""
        x.eint1(n): exponential integral E1(x):

        .. math::

                         \int_{x}^{\infty} \frac{e^{-t}}{t} dt


        If n is present, output the vector [eint1(x), eint1(2\*x), ...,
        eint1(n\*x)]. This is faster than repeatedly calling eint1(i\*x).

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        REFERENCE:

        - See page 262, Prop 5.6.12, of Cohen's book "A Course in
          Computational Algebraic Number Theory".

        EXAMPLES:
        """
        pari_catch_sig_on()
        if n <= 0:
            return P.new_gen(eint1(x.g, pbw(precision)))
        else:
            return P.new_gen(veceint1(x.g, stoi(n), pbw(precision)))

    def erfc(gen x, precision=0):
        r"""
        Return the complementary error function:

        .. math::

            (2/\sqrt{\pi}) \int_{x}^{\infty} e^{-t^2} dt.



        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(1).erfc()
            0.157299207050285
        """
        pari_catch_sig_on()
        return P.new_gen(gerfc(x.g, pbw(precision)))

    def eta(gen x, flag=0, precision=0):
        r"""
        x.eta(flag=0): if flag=0, `\eta` function without the
        `q^{1/24}`; otherwise `\eta` of the complex number
        `x` in the upper half plane intelligently computed using
        `\mathrm{SL}(2,\ZZ)` transformations.

        DETAILS: This functions computes the following. If the input
        `x` is a complex number with positive imaginary part, the
        result is `\prod_{n=1}^{\infty} (q-1^n)`, where
        `q=e^{2 i \pi x}`. If `x` is a power series
        (or can be converted to a power series) with positive valuation,
        the result is `\prod_{n=1}^{\infty} (1-x^n)`.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: C.<i> = ComplexField()
            sage: pari(i).eta()
            0.998129069925959
        """
        pari_catch_sig_on()
        if flag == 1:
            return P.new_gen(trueeta(x.g, pbw(precision)))
        return P.new_gen(eta(x.g, pbw(precision)))

    def exp(gen self, precision=0):
        """
        x.exp(): exponential of x.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(0).exp()
            1.00000000000000
            sage: pari(1).exp()
            2.71828182845905
            sage: pari('x+O(x^8)').exp()
            1 + x + 1/2*x^2 + 1/6*x^3 + 1/24*x^4 + 1/120*x^5 + 1/720*x^6 + 1/5040*x^7 + O(x^8)
        """
        pari_catch_sig_on()
        return P.new_gen(gexp(self.g, pbw(precision)))

    def gamma(gen s, precision=0):
        """
        s.gamma(precision): Gamma function at s.

        If `s` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `s` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(2).gamma()
            1.00000000000000
            sage: pari(5).gamma()
            24.0000000000000
            sage: C.<i> = ComplexField()
            sage: pari(1+i).gamma()
            0.498015668118356 - 0.154949828301811*I

        TESTS::

            sage: pari(-1).gamma()
            Traceback (most recent call last):
            ...
            PariError:  (5)
        """
        pari_catch_sig_on()
        return P.new_gen(ggamma(s.g, pbw(precision)))

    def gammah(gen s, precision=0):
        """
        s.gammah(): Gamma function evaluated at the argument x+1/2.

        If `s` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `s` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(2).gammah()
            1.32934038817914
            sage: pari(5).gammah()
            52.3427777845535
            sage: C.<i> = ComplexField()
            sage: pari(1+i).gammah()
            0.575315188063452 + 0.0882106775440939*I
        """
        pari_catch_sig_on()
        return P.new_gen(ggamd(s.g, pbw(precision)))

    def hyperu(gen a, b, x, precision=0):
        r"""
        a.hyperu(b,x): U-confluent hypergeometric function.

        If `a`, `b`, or `x` is an exact argument,
        it is first converted to a real or complex number using the
        optional parameter precision (in bits). If the arguments are
        inexact (e.g. real), the smallest of their precisions is used in
        the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(1).hyperu(2,3)
            0.333333333333333
        """
        t0GEN(b)
        t1GEN(x)
        pari_catch_sig_on()
        return P.new_gen(hyperu(a.g, t0, t1, pbw(precision)))

    def incgam(gen s, x, y=None, precision=0):
        r"""
        s.incgam(x, y, precision): incomplete gamma function. y is optional
        and is the precomputed value of gamma(s).

        If `s` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `s` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: C.<i> = ComplexField()
            sage: pari(1+i).incgam(3-i)
            -0.0458297859919946 + 0.0433696818726677*I
        """
        t0GEN(x)
        pari_catch_sig_on()
        if y is None:
            return P.new_gen(incgam(s.g, t0, pbw(precision)))
        else:
            t1GEN(y)
            return P.new_gen(incgam0(s.g, t0, t1, pbw(precision)))

    def incgamc(gen s, x, precision=0):
        r"""
        s.incgamc(x): complementary incomplete gamma function.

        The arguments `x` and `s` are complex numbers such
        that `s` is not a pole of `\Gamma` and
        `|x|/(|s|+1)` is not much larger than `1`
        (otherwise, the convergence is very slow). The function returns the
        value of the integral
        `\int_{0}^{x} e^{-t} t^{s-1} dt.`

        If `s` or `x` is an exact argument, it is first
        converted to a real or complex number using the optional parameter
        precision (in bits). If the arguments are inexact (e.g. real), the
        smallest of their precisions is used in the computation, and the
        parameter precision is ignored.

        EXAMPLES::

            sage: pari(1).incgamc(2)
            0.864664716763387
        """
        t0GEN(x)
        pari_catch_sig_on()
        return P.new_gen(incgamc(s.g, t0, pbw(precision)))

    def log(gen x, precision=0):
        r"""
        x.log(): natural logarithm of x.

        This function returns the principal branch of the natural logarithm
        of `x`, i.e., the branch such that
        `\Im(\log(x)) \in ]-\pi, \pi].` The result is
        complex (with imaginary part equal to `\pi`) if
        `x\in \RR` and `x<0`. In general, the algorithm uses
        the formula

        .. math::

                         \log(x) \simeq \frac{\pi}{2{\rm agm}(1,4/s)} - m\log(2),


        if `s=x 2^m` is large enough. (The result is exact to
        `B` bits provided that `s>2^{B/2}`.) At low
        accuracies, this function computes `\log` using the series
        expansion near `1`.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        Note that `p`-adic arguments can also be given as input,
        with the convention that `\log(p)=0`. Hence, in particular,
        `\exp(\log(x))/x` is not in general equal to `1`
        but instead to a `(p-1)`-st root of unity (or
        `\pm 1` if `p=2`) times a power of `p`.

        EXAMPLES::

            sage: pari(5).log()
            1.60943791243410
            sage: C.<i> = ComplexField()
            sage: pari(i).log()
            0.E-19 + 1.57079632679490*I
        """
        pari_catch_sig_on()
        return P.new_gen(glog(x.g, pbw(precision)))

    def lngamma(gen x, precision=0):
        r"""
        This method is deprecated, please use :meth:`.log_gamma` instead.

        See the :meth:`.log_gamma` method for documentation and examples.

        EXAMPLES::

            sage: pari(100).lngamma()
            doctest:...: DeprecationWarning: The method lngamma() is deprecated. Use log_gamma() instead.
            See http://trac.sagemath.org/6992 for details.
            359.134205369575
        """
        from sage.misc.superseded import deprecation
        deprecation(6992, "The method lngamma() is deprecated. Use log_gamma() instead.")
        return x.log_gamma(precision)

    def log_gamma(gen x, precision=0):
        r"""
        Logarithm of the gamma function of x.

        This function returns the principal branch of the logarithm of the
        gamma function of `x`. The function
        `\log(\Gamma(x))` is analytic on the complex plane with
        non-positive integers removed. This function can have much larger
        inputs than `\Gamma` itself.

        The `p`-adic analogue of this function is unfortunately not
        implemented.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(100).log_gamma()
            359.134205369575
        """
        pari_catch_sig_on()
        return P.new_gen(glngamma(x.g, pbw(precision)))

    def polylog(gen x, long m, flag=0, precision=0):
        """
        x.polylog(m,flag=0): m-th polylogarithm of x. flag is optional, and
        can be 0: default, 1: D_m -modified m-th polylog of x, 2:
        D_m-modified m-th polylog of x, 3: P_m-modified m-th polylog of
        x.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        TODO: Add more explanation, copied from the PARI manual.

        EXAMPLES::

            sage: pari(10).polylog(3)
            5.64181141475134 - 8.32820207698027*I
            sage: pari(10).polylog(3,0)
            5.64181141475134 - 8.32820207698027*I
            sage: pari(10).polylog(3,1)
            0.523778453502411
            sage: pari(10).polylog(3,2)
            -0.400459056163451
        """
        pari_catch_sig_on()
        return P.new_gen(polylog0(m, x.g, flag, pbw(precision)))

    def psi(gen x, precision=0):
        r"""
        x.psi(): psi-function at x.

        Return the `\psi`-function of `x`, i.e., the
        logarithmic derivative `\Gamma'(x)/\Gamma(x)`.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(1).psi()
            -0.577215664901533
        """
        pari_catch_sig_on()
        return P.new_gen(gpsi(x.g, pbw(precision)))

    def sin(gen x, precision=0):
        """
        x.sin(): The sine of x.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(1).sin()
            0.841470984807897
            sage: C.<i> = ComplexField()
            sage: pari(1+i).sin()
            1.29845758141598 + 0.634963914784736*I
        """
        pari_catch_sig_on()
        return P.new_gen(gsin(x.g, pbw(precision)))

    def sinh(gen x, precision=0):
        """
        The hyperbolic sine function.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(0).sinh()
            0.E-19
            sage: C.<i> = ComplexField()
            sage: pari(1+i).sinh()
            0.634963914784736 + 1.29845758141598*I
        """
        pari_catch_sig_on()
        return P.new_gen(gsh(x.g, pbw(precision)))

    def sqr(gen x):
        """
        x.sqr(): square of x. Faster than, and most of the time (but not
        always - see the examples) identical to x\*x.

        EXAMPLES::

            sage: pari(2).sqr()
            4

        For `2`-adic numbers, x.sqr() may not be identical to x\*x
        (squaring a `2`-adic number increases its precision)::

            sage: pari("1+O(2^5)").sqr()
            1 + O(2^6)
            sage: pari("1+O(2^5)")*pari("1+O(2^5)")
            1 + O(2^5)

        However::

            sage: x = pari("1+O(2^5)"); x*x
            1 + O(2^6)
        """
        pari_catch_sig_on()
        return P.new_gen(gsqr(x.g))


    def sqrt(gen x, precision=0):
        """
        x.sqrt(precision): The square root of x.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(2).sqrt()
            1.41421356237310
        """
        pari_catch_sig_on()
        return P.new_gen(gsqrt(x.g, pbw(precision)))

    def sqrtn(gen x, n, precision=0):
        r"""
        x.sqrtn(n): return the principal branch of the n-th root of x,
        i.e., the one such that
        `\arg(\sqrt(x)) \in ]-\pi/n, \pi/n]`. Also returns a second
        argument which is a suitable root of unity allowing one to recover
        all the other roots. If it was not possible to find such a number,
        then this second return value is 0. If the argument is present and
        no square root exists, return 0 instead of raising an error.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        .. note::

           intmods (modulo a prime) and `p`-adic numbers are
           allowed as arguments.

        INPUT:


        -  ``x`` - gen

        -  ``n`` - integer


        OUTPUT:


        -  ``gen`` - principal n-th root of x

        -  ``gen`` - root of unity z that gives the other
           roots


        EXAMPLES::

            sage: s, z = pari(2).sqrtn(5)
            sage: z
            0.309016994374947 + 0.951056516295154*I
            sage: s
            1.14869835499704
            sage: s^5
            2.00000000000000
            sage: z^5
            1.00000000000000 + 5.42101086 E-19*I        # 32-bit
            1.00000000000000 + 5.96311194867027 E-19*I  # 64-bit
            sage: (s*z)^5
            2.00000000000000 + 1.409462824 E-18*I       # 32-bit
            2.00000000000000 + 9.21571846612679 E-19*I  # 64-bit
        """
        # TODO: ???  lots of good examples in the PARI docs ???
        cdef GEN zetan
        t0GEN(n)
        pari_catch_sig_on()
        ans = P.new_gen_noclear(gsqrtn(x.g, t0, &zetan, pbw(precision)))
        return ans, P.new_gen(zetan)

    def tan(gen x, precision=0):
        """
        x.tan() - tangent of x

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(2).tan()
            -2.18503986326152
            sage: C.<i> = ComplexField()
            sage: pari(i).tan()
            0.E-19 + 0.761594155955765*I
        """
        pari_catch_sig_on()
        return P.new_gen(gtan(x.g, pbw(precision)))

    def tanh(gen x, precision=0):
        """
        x.tanh() - hyperbolic tangent of x

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(1).tanh()
            0.761594155955765
            sage: C.<i> = ComplexField()
            sage: z = pari(i); z
            0.E-19 + 1.00000000000000*I
            sage: result = z.tanh()
            sage: result.real() <= 1e-18
            True
            sage: result.imag()
            1.55740772465490
        """
        pari_catch_sig_on()
        return P.new_gen(gth(x.g, pbw(precision)))

    def teichmuller(gen x):
        r"""
        teichmuller(x): teichmuller character of p-adic number x.

        This is the unique `(p-1)`-st root of unity congruent to
        `x/p^{v_p(x)}` modulo `p`.

        EXAMPLES::

            sage: pari('2+O(7^5)').teichmuller()
            2 + 4*7 + 6*7^2 + 3*7^3 + O(7^5)
        """
        pari_catch_sig_on()
        return P.new_gen(teich(x.g))

    def theta(gen q, z, precision=0):
        """
        q.theta(z): Jacobi sine theta-function.

        If `q` or `z` is an exact argument, it is first
        converted to a real or complex number using the optional parameter
        precision (in bits). If the arguments are inexact (e.g. real), the
        smallest of their precisions is used in the computation, and the
        parameter precision is ignored.

        EXAMPLES::

            sage: pari(0.5).theta(2)
            1.63202590295260
        """
        t0GEN(z)
        pari_catch_sig_on()
        return P.new_gen(theta(q.g, t0, pbw(precision)))

    def thetanullk(gen q, long k, precision=0):
        """
        q.thetanullk(k): return the k-th derivative at z=0 of theta(q,z).

        If `q` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `q` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        EXAMPLES::

            sage: pari(0.5).thetanullk(1)
            0.548978532560341
        """
        pari_catch_sig_on()
        return P.new_gen(thetanullk(q.g, k, pbw(precision)))

    def weber(gen x, flag=0, precision=0):
        r"""
        x.weber(flag=0): One of Weber's f functions of x. flag is optional,
        and can be 0: default, function
        f(x)=exp(-i\*Pi/24)\*eta((x+1)/2)/eta(x) such that
        `j=(f^{24}-16)^3/f^{24}`, 1: function f1(x)=eta(x/2)/eta(x)
        such that `j=(f1^24+16)^3/f2^{24}`, 2: function
        f2(x)=sqrt(2)\*eta(2\*x)/eta(x) such that
        `j=(f2^{24}+16)^3/f2^{24}`.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        TODO: Add further explanation from PARI manual.

        EXAMPLES::

            sage: C.<i> = ComplexField()
            sage: pari(i).weber()
            1.18920711500272 + 0.E-19*I                 # 32-bit
            1.18920711500272 + 2.71050543121376 E-20*I  # 64-bit
            sage: pari(i).weber(1)
            1.09050773266526 + 0.E-19*I
            sage: pari(i).weber(2)
            1.09050773266526
        """
        pari_catch_sig_on()
        return P.new_gen(weber0(x.g, flag, pbw(precision)))

    def zeta(gen s, precision=0):
        """
        zeta(s): zeta function at s with s a complex or a p-adic number.

        If `s` is a complex number, this is the Riemann zeta
        function `\zeta(s)=\sum_{n\geq 1} n^{-s}`, computed either
        using the Euler-Maclaurin summation formula (if `s` is not
        an integer), or using Bernoulli numbers (if `s` is a
        negative integer or an even nonnegative integer), or using modular
        forms (if `s` is an odd nonnegative integer).

        If `s` is a `p`-adic number, this is the
        Kubota-Leopoldt zeta function, i.e. the unique continuous
        `p`-adic function on the `p`-adic integers that
        interpolates the values of `(1-p^{-k})\zeta(k)` at negative
        integers `k` such that `k\equiv 1\pmod{p-1}` if
        `p` is odd, and at odd `k` if `p=2`.

        If `x` is an exact argument, it is first converted to a
        real or complex number using the optional parameter precision (in
        bits). If `x` is inexact (e.g. real), its own precision is
        used in the computation, and the parameter precision is ignored.

        INPUT:


        -  ``s`` - gen (real, complex, or p-adic number)


        OUTPUT:


        -  ``gen`` - value of zeta at s.


        EXAMPLES::

            sage: pari(2).zeta()
            1.64493406684823
            sage: x = RR(pi)^2/6
            sage: pari(x)
            1.64493406684823
            sage: pari(3).zeta()
            1.20205690315959
            sage: pari('1+5*7+2*7^2+O(7^3)').zeta()
            4*7^-2 + 5*7^-1 + O(7^0)
        """
        pari_catch_sig_on()
        return P.new_gen(gzeta(s.g, pbw(precision)))

    ###########################################
    # 4: NUMBER THEORETICAL functions
    ###########################################

    def bezout(gen x, y):
        cdef gen u, v, g
        cdef GEN U, V, G
        t0GEN(y)
        pari_catch_sig_on()
        G = gbezout(x.g, t0, &U, &V)
        g = P.new_gen_noclear(G)
        u = P.new_gen_noclear(U)
        v = P.new_gen(V)
        return g, u, v

    def binomial(gen x, long k):
        """
        binomial(x, k): return the binomial coefficient "x choose k".

        INPUT:


        -  ``x`` - any PARI object (gen)

        -  ``k`` - integer


        EXAMPLES::

            sage: pari(6).binomial(2)
            15
            sage: pari('x+1').binomial(3)
            1/6*x^3 - 1/6*x
            sage: pari('2+x+O(x^2)').binomial(3)
            1/3*x + O(x^2)
        """
        pari_catch_sig_on()
        return P.new_gen(binomial(x.g, k))

    def contfrac(gen x, b=0, long lmax=0):
        """
        contfrac(x,b,lmax): continued fraction expansion of x (x rational,
        real or rational function). b and lmax are both optional, where b
        is the vector of numerators of the continued fraction, and lmax is
        a bound for the number of terms in the continued fraction
        expansion.
        """
        t0GEN(b)
        pari_catch_sig_on()
        return P.new_gen(contfrac0(x.g, t0, lmax))

    def contfracpnqn(gen x, b=0, long lmax=0):
        """
        contfracpnqn(x): [p_n,p_n-1; q_n,q_n-1] corresponding to the
        continued fraction x.
        """
        pari_catch_sig_on()
        return P.new_gen(pnqn(x.g))

    def ffgen(gen T, v=-1):
        r"""
        Return the generator `g=x \bmod T` of the finite field defined
        by the polynomial `T`.

        INPUT:

        - ``T`` -- a gen of type t_POL with coefficients of type t_INTMOD:
                   a polynomial over a prime finite field

        - ``v`` -- string: a variable name or -1 (optional)

        If `v` is a string, then `g` will be a polynomial in `v`, else the
        variable of the polynomial `T` is used.

        EXAMPLES::

            sage: x = GF(2)['x'].gen()
            sage: pari(x^2+x+2).ffgen()
            x
            sage: pari(x^2+x+1).ffgen('a')
            a
        """
        pari_catch_sig_on()
        return P.new_gen(ffgen(T.g, P.get_var(v)))

    def ffinit(gen p, long n, v=-1):
        r"""
        Return a monic irreducible polynomial `g` of degree `n` over the
        finite field of `p` elements.

        INPUT:

        - ``p`` -- a gen of type t_INT: a prime number

        - ``n`` -- integer: the degree of the polynomial

        - ``v`` -- string: a variable name or -1 (optional)

        If `v \geq 0', then `g` will be a polynomial in `v`, else the
        variable `x` is used.

        EXAMPLES::

            sage: pari(7).ffinit(11)
            Mod(1, 7)*x^11 + Mod(1, 7)*x^10 + Mod(4, 7)*x^9 + Mod(5, 7)*x^8 + Mod(1, 7)*x^7 + Mod(1, 7)*x^2 + Mod(1, 7)*x + Mod(6, 7)
            sage: pari(2003).ffinit(3)
            Mod(1, 2003)*x^3 + Mod(1, 2003)*x^2 + Mod(1993, 2003)*x + Mod(1995, 2003)
        """
        pari_catch_sig_on()
        return P.new_gen(ffinit(p.g, n, P.get_var(v)))

    def fibonacci(gen x):
        r"""
        Return the Fibonacci number of index x.

        EXAMPLES::

            sage: pari(18).fibonacci()
            2584
            sage: [pari(n).fibonacci() for n in range(10)]
            [0, 1, 1, 2, 3, 5, 8, 13, 21, 34]
        """
        pari_catch_sig_on()
        return P.new_gen(fibo(long(x)))


    def gcd(gen x, y, long flag=0):
        """
        gcd(x,y,flag=0): greatest common divisor of x and y. flag is
        optional, and can be 0: default, 1: use the modular gcd algorithm
        (x and y must be polynomials), 2 use the subresultant algorithm (x
        and y must be polynomials)
        """
        t0GEN(y)
        pari_catch_sig_on()
        return P.new_gen(ggcd0(x.g, t0))

    def issquare(gen x, find_root=False):
        """
        issquare(x,n): true(1) if x is a square, false(0) if not. If
        find_root is given, also returns the exact square root if it was
        computed.
        """
        cdef GEN G, t
        cdef gen g
        pari_catch_sig_on()
        if find_root:
            t = gissquareall(x.g, &G)
            v = bool(P.new_gen_noclear(t))
            if v:
                return v, P.new_gen(G)
            else:
                pari_catch_sig_off()
                return v, None
        else:
            return P.new_gen(gissquare(x.g))


    def issquarefree(gen self):
        """
        EXAMPLES::

            sage: pari(10).issquarefree()
            True
            sage: pari(20).issquarefree()
            False
        """
        pari_catch_sig_on()
        t = bool(issquarefree(self.g))
        pari_catch_sig_off()
        return t

    def lcm(gen x, y):
        """
        Return the least common multiple of x and y. EXAMPLES::

            sage: pari(10).lcm(15)
            30
        """
        t0GEN(y)
        pari_catch_sig_on()
        return P.new_gen(glcm(x.g, t0))

    def numdiv(gen n):
        """
        Return the number of divisors of the integer n.

        EXAMPLES::

            sage: pari(10).numdiv()
            4
        """
        pari_catch_sig_on()
        return P.new_gen(gnumbdiv(n.g))

    def phi(gen n):
        """
        Return the Euler phi function of n. EXAMPLES::

            sage: pari(10).phi()
            4
        """
        pari_catch_sig_on()
        return P.new_gen(geulerphi(n.g))

    def primepi(gen self):
        """
        Return the number of primes less than or equal to self.

        EXAMPLES::

            sage: pari(7).primepi()
            4
            sage: pari(100).primepi()
            25
            sage: pari(1000).primepi()
            168
            sage: pari(100000).primepi()
            9592
            sage: pari(0).primepi()
            0
            sage: pari(-15).primepi()
            0
            sage: pari(500509).primepi()
            41581
        """
        global num_primes
        pari_catch_sig_on()
        if self > num_primes:
            P.init_primes(self + 10)
        if signe(self.g) != 1:
            pari_catch_sig_off()
            return P.PARI_ZERO
        return P.new_gen(primepi(self.g))

    def sumdiv(gen n):
        """
        Return the sum of the divisors of `n`.

        EXAMPLES::

            sage: pari(10).sumdiv()
            18
        """
        pari_catch_sig_on()
        return P.new_gen(sumdiv(n.g))

    def sumdivk(gen n, long k):
        """
        Return the sum of the k-th powers of the divisors of n.

        EXAMPLES::

            sage: pari(10).sumdivk(2)
            130
        """
        pari_catch_sig_on()
        return P.new_gen(sumdivk(n.g, k))

    def xgcd(gen x, y):
        """
        Returns u,v,d such that d=gcd(x,y) and u\*x+v\*y=d.

        EXAMPLES::

            sage: pari(10).xgcd(15)
            (5, -1, 1)
        """
        return x.bezout(y)


    ##################################################
    # 5: Elliptic curve functions
    ##################################################

    def ellinit(self, int flag=0, precision=0):
        """
        Return the Pari elliptic curve object with Weierstrass coefficients
        given by self, a list with 5 elements.

        INPUT:


        -  ``self`` - a list of 5 coefficients

        -  ``flag (optional, default: 0)`` - if 0, ask for a
           Pari ell structure with 19 components; if 1, ask for a Pari sell
           structure with only the first 13 components

        -  ``precision (optional, default: 0)`` - the real
           precision to be used in the computation of the components of the
           Pari (s)ell structure; if 0, use the default 53 bits.

           .. note::

              the parameter precision in :meth:`.ellinit` controls not only
              the real precision of the resulting (s)ell structure,
              but also the precision of most subsequent computations
              with this elliptic curve.  You should therefore set it
              from the start to the value you require.


        OUTPUT:


        -  ``gen`` - either a Pari ell structure with 19
           components (if flag=0), or a Pari sell structure with 13 components
           (if flag=1)


        EXAMPLES: An elliptic curve with integer coefficients::

            sage: e = pari([0,1,0,1,0]).ellinit(); e
            [0, 1, 0, 1, 0, 4, 2, 0, -1, -32, 224, -48, 2048/3, [0.E-28, -0.500000000000000 - 0.866025403784439*I, -0.500000000000000 + 0.866025403784439*I]~, 3.37150070962519, -1.68575035481260 - 2.15651564749964*I, 1.37451455785745 - 1.084202173 E-19*I,      -0.687257278928726 + 0.984434956803824*I, 7.27069403586288] # 32-bit
            [0, 1, 0, 1, 0, 4, 2, 0, -1, -32, 224, -48, 2048/3, [0.E-38, -0.500000000000000 - 0.866025403784439*I, -0.500000000000000 + 0.866025403784439*I]~, 3.37150070962519, -1.68575035481260 - 2.15651564749964*I, 1.37451455785745 - 5.42101086242752 E-19*I, -0.687257278928726 + 0.984434956803824*I, 7.27069403586288] # 64-bit

        Its inexact components have the default precision of 53 bits::

            sage: RR(e[14])
            3.37150070962519

        We can compute this to higher precision::

            sage: R = RealField(150)
            sage: e = pari([0,1,0,1,0]).ellinit(precision=150)
            sage: R(e[14])
            3.3715007096251920857424073155981539790016018

        Using flag=1 returns a short elliptic curve Pari object::

            sage: pari([0,1,0,1,0]).ellinit(flag=1)
            [0, 1, 0, 1, 0, 4, 2, 0, -1, -32, 224, -48, 2048/3]

        The coefficients can be any ring elements that convert to Pari::

            sage: pari([0,1/2,0,-3/4,0]).ellinit(flag=1)
            [0, 1/2, 0, -3/4, 0, 2, -3/2, 0, -9/16, 40, -116, 117/4, 256000/117]
            sage: pari([0,0.5,0,-0.75,0]).ellinit(flag=1)
            [0, 0.500000000000000, 0, -0.750000000000000, 0, 2.00000000000000, -1.50000000000000, 0, -0.562500000000000, 40.0000000000000, -116.000000000000, 29.2500000000000, 2188.03418803419]
            sage: pari([0,I,0,1,0]).ellinit(flag=1)
            [0, I, 0, 1, 0, 4*I, 2, 0, -1, -64, 352*I, -80, 16384/5]
            sage: pari([0,x,0,2*x,1]).ellinit(flag=1)
            [0, x, 0, 2*x, 1, 4*x, 4*x, 4, -4*x^2 + 4*x, 16*x^2 - 96*x, -64*x^3 + 576*x^2 - 864, 64*x^4 - 576*x^3 + 576*x^2 - 432, (256*x^6 - 4608*x^5 + 27648*x^4 - 55296*x^3)/(4*x^4 - 36*x^3 + 36*x^2 - 27)]
        """
        pari_catch_sig_on()
        return P.new_gen(ellinit0(self.g, flag, pbw(precision)))

    def ellglobalred(self):
        """
        e.ellglobalred(): return information related to the global minimal
        model of the elliptic curve e.

        INPUT:


        -  ``e`` - elliptic curve (returned by ellinit)


        OUTPUT:


        -  ``gen`` - the (arithmetic) conductor of e

        -  ``gen`` - a vector giving the coordinate change over
           Q from e to its minimal integral model (see also ellminimalmodel)

        -  ``gen`` - the product of the local Tamagawa numbers
           of e


        EXAMPLES::

            sage: e = pari([0, 5, 2, -1, 1]).ellinit()
            sage: e.ellglobalred()
            [20144, [1, -2, 0, -1], 1]
            sage: e = pari(EllipticCurve('17a').a_invariants()).ellinit()
            sage: e.ellglobalred()
            [17, [1, 0, 0, 0], 4]
        """
        pari_catch_sig_on()
        return self.new_gen(ellglobalred(self.g))

    def elladd(self, z0, z1):
        """
        e.elladd(z0, z1): return the sum of the points z0 and z1 on this
        elliptic curve.

        INPUT:


        -  ``e`` - elliptic curve E

        -  ``z0`` - point on E

        -  ``z1`` - point on E


        OUTPUT: point on E

        EXAMPLES: First we create an elliptic curve::

            sage: e = pari([0, 1, 1, -2, 0]).ellinit()
            sage: str(e)[:65]   # first part of output
            '[0, 1, 1, -2, 0, 4, -4, 1, -3, 112, -856, 389, 1404928/389, [0.90'

        Next we add two points on the elliptic curve. Notice that the
        Python lists are automatically converted to PARI objects so you
        don't have to do that explicitly in your code.

        ::

            sage: e.elladd([1,0], [-1,1])
            [-3/4, -15/8]
        """
        t0GEN(z0); t1GEN(z1)
        pari_catch_sig_on()
        return self.new_gen(addell(self.g, t0, t1))

    def ellak(self, n):
        r"""
        e.ellak(n): Returns the coefficient `a_n` of the
        `L`-function of the elliptic curve e, i.e. the
        `n`-th Fourier coefficient of the weight 2 newform
        associated to e (according to Shimura-Taniyama).

            The curve `e` *must* be a medium or long vector of the type
            given by ellinit. For this function to work for every n and not
            just those prime to the conductor, e must be a minimal Weierstrass
            equation. If this is not the case, use the function ellminimalmodel
            first before using ellak (or you will get INCORRECT RESULTS!)


        INPUT:


        -  ``e`` - a PARI elliptic curve.

        -  ``n`` - integer.


        EXAMPLES::

            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.ellak(6)
            2
            sage: e.ellak(2005)
            2
            sage: e.ellak(-1)
            0
            sage: e.ellak(0)
            0
        """
        t0GEN(n)
        pari_catch_sig_on()
        return self.new_gen(akell(self.g, t0))


    def ellan(self, long n, python_ints=False):
        """
        Return the first `n` Fourier coefficients of the modular
        form attached to this elliptic curve. See ellak for more details.

        INPUT:


        -  ``n`` - a long integer

        -  ``python_ints`` - bool (default is False); if True,
           return a list of Python ints instead of a PARI gen wrapper.


        EXAMPLES::

            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.ellan(3)
            [1, -2, -1]
            sage: e.ellan(20)
            [1, -2, -1, 2, 1, 2, -2, 0, -2, -2, 1, -2, 4, 4, -1, -4, -2, 4, 0, 2]
            sage: e.ellan(-1)
            []
            sage: v = e.ellan(10, python_ints=True); v
            [1, -2, -1, 2, 1, 2, -2, 0, -2, -2]
            sage: type(v)
            <type 'list'>
            sage: type(v[0])
            <type 'int'>
        """
        pari_catch_sig_on()
        cdef GEN g
        if python_ints:
            g = anell(self.g, n)
            v = [gtolong(<GEN> g[i+1]) for i in range(glength(g))]
            (<PariInstance>pari).clear_stack()
            return v
        else:
            return self.new_gen(anell(self.g, n))

    def ellanalyticrank(self, long precision = 0):
        r"""
        Returns a 2-component vector with the order of vanishing at
        `s = 1` of the L-function of the elliptic curve and the value
        of the first non-zero derivative.

        EXAMPLE::

            sage: E = EllipticCurve('389a1')
            sage: pari(E).ellanalyticrank()
            [2, 1.51863300057685]
        """
        pari_catch_sig_on()
        return self.new_gen(ellanalyticrank(self.g, <GEN>0, pbw(precision)))

    def ellap(self, p):
        r"""
        e.ellap(p): Returns the prime-indexed coefficient `a_p` of the
        `L`-function of the elliptic curve `e`, i.e. the `p`-th Fourier
        coefficient of the newform attached to e.

        The computation uses the Shanks--Mestre method, or the SEA
        algorithm.

        .. WARNING::

            For this function to work for every n and not just those prime
            to the conductor, e must be a minimal Weierstrass equation.
            If this is not the case, use the function ellminimalmodel first
            before using ellap (or you will get INCORRECT RESULTS!)


        INPUT:


        -  ``e`` - a PARI elliptic curve.

        -  ``p`` - prime integer


        EXAMPLES::

            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.ellap(2)
            -2
            sage: e.ellap(2003)
            4
            sage: e.ellak(-1)
            0
        """
        t0GEN(p)
        pari_catch_sig_on()
        return self.new_gen(ellap(self.g, t0))


    def ellaplist(self, long n, python_ints=False):
        r"""
        e.ellaplist(n): Returns a PARI list of all the prime-indexed
        coefficients `a_p` (up to n) of the `L`-function
        of the elliptic curve `e`, i.e. the Fourier coefficients of
        the newform attached to `e`.

        INPUT:


        -  ``n`` - a long integer

        -  ``python_ints`` - bool (default is False); if True,
           return a list of Python ints instead of a PARI gen wrapper.


            The curve e must be a medium or long vector of the type given by
            ellinit. For this function to work for every n and not just those
            prime to the conductor, e must be a minimal Weierstrass equation.
            If this is not the case, use the function ellminimalmodel first
            before using ellaplist (or you will get INCORRECT RESULTS!)


        INPUT:


        -  ``e`` - a PARI elliptic curve.

        -  ``n`` - an integer


        EXAMPLES::

            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: v = e.ellaplist(10); v
            [-2, -1, 1, -2]
            sage: type(v)
            <type 'sage.libs.pari.gen.gen'>
            sage: v.type()
            't_VEC'
            sage: e.ellan(10)
            [1, -2, -1, 2, 1, 2, -2, 0, -2, -2]
            sage: v = e.ellaplist(10, python_ints=True); v
            [-2, -1, 1, -2]
            sage: type(v)
            <type 'list'>
            sage: type(v[0])
            <type 'int'>
        """
        # 1. make a table of primes up to n.
        pari_catch_sig_on()
        if n < 2:
            return self.new_gen(zerovec(0))
        cdef GEN g
        pari.init_primes(n+1)
        t0GEN(n)
        g = primes(gtolong(primepi(t0)))

        # 2. Replace each prime in the table by ellap of it.
        cdef long i

        if python_ints:
            v = [gtolong(ellap(self.g, <GEN> g[i+1])) \
                        for i in range(glength(g))]
            (<PariInstance>pari).clear_stack()
            return v
        else:
            for i from 0 <= i < glength(g):
                g[i+1] = <long> ellap(self.g, <GEN> g[i+1])
            return self.new_gen(g)


    def ellbil(self, z0, z1):
        """
        e.ellbil(z0, z1): return the value of the canonical bilinear form
        on z0 and z1.

        INPUT:


        -  ``e`` - elliptic curve (assumed integral given by a
           minimal model, as returned by ellminimalmodel)

        -  ``z0, z1`` - rational points on e


        EXAMPLES::

            sage: e = pari([0,1,1,-2,0]).ellinit().ellminimalmodel()[0]
            sage: e.ellbil([1, 0], [-1, 1])
            0.418188984498861
        """
        t0GEN(z0); t1GEN(z1)
        pari_catch_sig_on()
        # the prec argument has no effect
        return self.new_gen(bilhell(self.g, t0, t1, prec))

    def ellchangecurve(self, ch):
        """
        e.ellchangecurve(ch): return the new model (equation) for the
        elliptic curve e given by the change of coordinates ch.

        The change of coordinates is specified by a vector ch=[u,r,s,t]; if
        `x'` and `y'` are the new coordinates, then
        `x = u^2 x' + r` and `y = u^3 y' + su^2 x' + t`.

        INPUT:


        -  ``e`` - elliptic curve

        -  ``ch`` - change of coordinates vector with 4
           entries


        EXAMPLES::

            sage: e = pari([1,2,3,4,5]).ellinit()
            sage: e.ellglobalred()
            [10351, [1, -1, 0, -1], 1]
            sage: f = e.ellchangecurve([1,-1,0,-1])
            sage: f[:5]
            [1, -1, 0, 4, 3]
        """
        t0GEN(ch)
        pari_catch_sig_on()
        return self.new_gen(ellchangecurve(self.g, t0))

    def elleta(self):
        """
        e.elleta(): return the vector [eta1,eta2] of quasi-periods
        associated with the period lattice e.omega() of the elliptic curve
        e.

        EXAMPLES::

            sage: e = pari([0,0,0,-82,0]).ellinit()
            sage: e.elleta()
            [3.60546360143265, 3.60546360143265*I]
            sage: w1,w2 = e.omega()
            sage: eta1, eta2 = e.elleta()
            sage: w1*eta2-w2*eta1
            6.28318530717959*I
            sage: w1*eta2-w2*eta1 == pari(2*pi*I)
            True
        """
        pari_catch_sig_on()
        # the prec argument has no effect
        return self.new_gen(elleta(self.g, prec))

    def ellheight(self, a, flag=2, precision=0):
        """
        e.ellheight(a, flag=2): return the global Neron-Tate height of the
        point a on the elliptic curve e.

        INPUT:


        -  ``e`` - elliptic curve over `\QQ`,
           assumed to be in a standard minimal integral model (as given by
           ellminimalmodel)

        -  ``a`` - rational point on e

        -  ``flag (optional)`` - specifies which algorithm to
           be used for computing the archimedean local height:

           -  ``0`` - uses sigma- and theta-functions and a trick
               due to J. Silverman

           -  ``1`` - uses Tate's `4^n` algorithm

           -  ``2`` - uses Mestre's AGM algorithm (this is the
              default, being faster than the other two)

        -  ``precision (optional)`` - the precision of the
           result, in bits.


        Note that in order to achieve the desired precision, the
        elliptic curve must have been created using ellinit with the
        desired precision.

        EXAMPLES::

            sage: e = pari([0,1,1,-2,0]).ellinit().ellminimalmodel()[0]
            sage: e.ellheight([1,0])
            0.476711659343740
            sage: e.ellheight([1,0], flag=0)
            0.476711659343740
            sage: e.ellheight([1,0], flag=1)
            0.476711659343740
        """
        t0GEN(a)
        pari_catch_sig_on()
        return self.new_gen(ellheight0(self.g, t0, flag, pbw(precision)))

    def ellheightmatrix(self, x):
        """
        e.ellheightmatrix(x): return the height matrix for the vector x of
        points on the elliptic curve e.

        In other words, it returns the Gram matrix of x with respect to the
        height bilinear form on e (see ellbil).

        INPUT:


        -  ``e`` - elliptic curve over `\QQ`,
           assumed to be in a standard minimal integral model (as given by
           ellminimalmodel)

        -  ``x`` - vector of rational points on e


        EXAMPLES::

            sage: e = pari([0,1,1,-2,0]).ellinit().ellminimalmodel()[0]
            sage: e.ellheightmatrix([[1,0], [-1,1]])
            [0.476711659343740, 0.418188984498861; 0.418188984498861, 0.686667083305587]
        """
        t0GEN(x)
        pari_catch_sig_on()
        # the argument prec has no effect
        return self.new_gen(mathell(self.g, t0, prec))

    def ellisoncurve(self, x):
        """
        e.ellisoncurve(x): return True if the point x is on the elliptic
        curve e, False otherwise.

        If the point or the curve have inexact coefficients, an attempt is
        made to take this into account.

        EXAMPLES::

            sage: e = pari([0,1,1,-2,0]).ellinit()
            sage: e.ellisoncurve([1,0])
            True
            sage: e.ellisoncurve([1,1])
            False
            sage: e.ellisoncurve([1,0.00000000000000001])
            False
            sage: e.ellisoncurve([1,0.000000000000000001])
            True
            sage: e.ellisoncurve([0])
            True
        """
        t0GEN(x)
        pari_catch_sig_on()
        t = bool(oncurve(self.g, t0) == 1)
        pari_catch_sig_off()
        return t

    def elllocalred(self, p):
        r"""
        e.elllocalred(p): computes the data of local reduction at the prime
        p on the elliptic curve e

        For more details on local reduction and Kodaira types, see IV.8 and
        IV.9 in J. Silverman's book "Advanced topics in the arithmetic of
        elliptic curves".

        INPUT:


        -  ``e`` - elliptic curve with coefficients in `\ZZ`

        -  ``p`` - prime number


        OUTPUT:


        -  ``gen`` - the exponent of p in the arithmetic
           conductor of e

        -  ``gen`` - the Kodaira type of e at p, encoded as an
           integer:

        -  ``1`` - type `I_0`: good reduction,
           nonsingular curve of genus 1

        -  ``2`` - type `II`: rational curve with a
           cusp

        -  ``3`` - type `III`: two nonsingular rational
           curves intersecting tangentially at one point

        -  ``4`` - type `IV`: three nonsingular
           rational curves intersecting at one point

        -  ``5`` - type `I_1`: rational curve with a
           node

        -  ``6 or larger`` - think of it as `4+v`, then
           it is type `I_v`: `v` nonsingular rational curves
           arranged as a `v`-gon

        -  ``-1`` - type `I_0^*`: nonsingular rational
           curve of multiplicity two with four nonsingular rational curves of
           multiplicity one attached

        -  ``-2`` - type `II^*`: nine nonsingular
           rational curves in a special configuration

        -  ``-3`` - type `III^*`: eight nonsingular
           rational curves in a special configuration

        -  ``-4`` - type `IV^*`: seven nonsingular
           rational curves in a special configuration

        -  ``-5 or smaller`` - think of it as `-4-v`,
           then it is type `I_v^*`: chain of `v+1`
           nonsingular rational curves of multiplicity two, with two
           nonsingular rational curves of multiplicity one attached at either
           end

        -  ``gen`` - a vector with 4 components, giving the
           coordinate changes done during the local reduction; if the first
           component is 1, then the equation for e was already minimal at p

        -  ``gen`` - the local Tamagawa number `c_p`


        EXAMPLES:

        Type `I_0`::

            sage: e = pari([0,0,0,0,1]).ellinit()
            sage: e.elllocalred(7)
            [0, 1, [1, 0, 0, 0], 1]

        Type `II`::

            sage: e = pari(EllipticCurve('27a3').a_invariants()).ellinit()
            sage: e.elllocalred(3)
            [3, 2, [1, -1, 0, 1], 1]

        Type `III`::

            sage: e = pari(EllipticCurve('24a4').a_invariants()).ellinit()
            sage: e.elllocalred(2)
            [3, 3, [1, 1, 0, 1], 2]

        Type `IV`::

            sage: e = pari(EllipticCurve('20a2').a_invariants()).ellinit()
            sage: e.elllocalred(2)
            [2, 4, [1, 1, 0, 1], 3]

        Type `I_1`::

            sage: e = pari(EllipticCurve('11a2').a_invariants()).ellinit()
            sage: e.elllocalred(11)
            [1, 5, [1, 0, 0, 0], 1]

        Type `I_2`::

            sage: e = pari(EllipticCurve('14a4').a_invariants()).ellinit()
            sage: e.elllocalred(2)
            [1, 6, [1, 0, 0, 0], 2]

        Type `I_6`::

            sage: e = pari(EllipticCurve('14a1').a_invariants()).ellinit()
            sage: e.elllocalred(2)
            [1, 10, [1, 0, 0, 0], 2]

        Type `I_0^*`::

            sage: e = pari(EllipticCurve('32a3').a_invariants()).ellinit()
            sage: e.elllocalred(2)
            [5, -1, [1, 1, 1, 0], 1]

        Type `II^*`::

            sage: e = pari(EllipticCurve('24a5').a_invariants()).ellinit()
            sage: e.elllocalred(2)
            [3, -2, [1, 2, 1, 4], 1]

        Type `III^*`::

            sage: e = pari(EllipticCurve('24a2').a_invariants()).ellinit()
            sage: e.elllocalred(2)
            [3, -3, [1, 2, 1, 4], 2]

        Type `IV^*`::

            sage: e = pari(EllipticCurve('20a1').a_invariants()).ellinit()
            sage: e.elllocalred(2)
            [2, -4, [1, 0, 1, 2], 3]

        Type `I_1^*`::

            sage: e = pari(EllipticCurve('24a1').a_invariants()).ellinit()
            sage: e.elllocalred(2)
            [3, -5, [1, 0, 1, 2], 4]

        Type `I_6^*`::

            sage: e = pari(EllipticCurve('90c2').a_invariants()).ellinit()
            sage: e.elllocalred(3)
            [2, -10, [1, 96, 1, 316], 4]
        """
        t0GEN(p)
        pari_catch_sig_on()
        return self.new_gen(elllocalred(self.g, t0))

    def elllseries(self, s, A=1):
        """
        e.elllseries(s, A=1): return the value of the `L`-series of
        the elliptic curve e at the complex number s.

        This uses an `O(N^{1/2})` algorithm in the conductor N of
        e, so it is impractical for large conductors (say greater than
        `10^{12}`).

        INPUT:


        -  ``e`` - elliptic curve defined over `\QQ`

        -  ``s`` - complex number

        -  ``A (optional)`` - cutoff point for the integral,
           which must be chosen close to 1 for best speed.


        EXAMPLES::

            sage: e = pari([0,1,1,-2,0]).ellinit()
            sage: e.elllseries(2.1)
            0.402838047956645
            sage: e.elllseries(1)   # random, close to 0
            1.822829333527862 E-19
            sage: e.elllseries(-2)
            0

        The following example differs for the last digit on 32 vs. 64 bit
        systems

        ::

            sage: e.elllseries(2.1, A=1.1)
            0.402838047956645
        """
        t0GEN(s); t1GEN(A)
        pari_catch_sig_on()
        # the argument prec has no effect
        return self.new_gen(elllseries(self.g, t0, t1, prec))

    def ellminimalmodel(self):
        """
        ellminimalmodel(e): return the standard minimal integral model of
        the rational elliptic curve e and the corresponding change of
        variables. INPUT:


        -  ``e`` - gen (that defines an elliptic curve)


        OUTPUT:


        -  ``gen`` - minimal model

        -  ``gen`` - change of coordinates


        EXAMPLES::

            sage: e = pari([1,2,3,4,5]).ellinit()
            sage: F, ch = e.ellminimalmodel()
            sage: F[:5]
            [1, -1, 0, 4, 3]
            sage: ch
            [1, -1, 0, -1]
            sage: e.ellchangecurve(ch)[:5]
            [1, -1, 0, 4, 3]
        """
        cdef GEN x, y
        cdef gen model, change
        cdef pari_sp t
        pari_catch_sig_on()
        x = ellminimalmodel(self.g, &y)
        change = self.new_gen_noclear(y)
        model = self.new_gen(x)
        return model, change

    def ellorder(self, x):
        """
        e.ellorder(x): return the order of the point x on the elliptic
        curve e (return 0 if x is not a torsion point)

        INPUT:


        -  ``e`` - elliptic curve defined over `\QQ`

        -  ``x`` - point on e


        EXAMPLES::

            sage: e = pari(EllipticCurve('65a1').a_invariants()).ellinit()

        A point of order two::

            sage: e.ellorder([0,0])
            2

        And a point of infinite order::

            sage: e.ellorder([1,0])
            0
        """
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(orderell(self.g, t0))

    def ellordinate(self, x):
        """
        e.ellordinate(x): return the `y`-coordinates of the points
        on the elliptic curve e having x as `x`-coordinate.

        INPUT:


        -  ``e`` - elliptic curve

        -  ``x`` - x-coordinate (can be a complex or p-adic
           number, or a more complicated object like a power series)


        EXAMPLES::

            sage: e = pari([0,1,1,-2,0]).ellinit()
            sage: e.ellordinate(0)
            [0, -1]
            sage: C.<i> = ComplexField()
            sage: e.ellordinate(i)
            [0.582203589721741 - 1.38606082464177*I, -1.58220358972174 + 1.38606082464177*I]
            sage: e.ellordinate(1+3*5^1+O(5^3))
            [4*5 + 5^2 + O(5^3), 4 + 3*5^2 + O(5^3)]
            sage: e.ellordinate('z+2*z^2+O(z^4)')
            [-2*z - 7*z^2 - 23*z^3 + O(z^4), -1 + 2*z + 7*z^2 + 23*z^3 + O(z^4)]
        """
        t0GEN(x)
        pari_catch_sig_on()
        # the prec argument has no effect
        return self.new_gen(ellordinate(self.g, t0, prec))

    def ellpointtoz(self, P, long precision=0):
        """
        e.ellpointtoz(P): return the complex number (in the fundamental
        parallelogram) corresponding to the point P on the elliptic curve
        e, under the complex uniformization of e given by the Weierstrass
        p-function.

        The complex number z returned by this function lies in the
        parallelogram formed by the real and complex periods of e, as given
        by e.omega().

        EXAMPLES::

            sage: e = pari([0,0,0,1,0]).ellinit()
            sage: e.ellpointtoz([0,0])
            1.85407467730137

        The point at infinity is sent to the complex number 0::

            sage: e.ellpointtoz([0])
            0
        """
        t0GEN(P)
        pari_catch_sig_on()
        return self.new_gen(zell(self.g, t0, pbw(precision)))

    def ellpow(self, z, n):
        """
        e.ellpow(z, n): return `n` times the point `z` on the elliptic
        curve `e`.

        INPUT:


        -  ``e`` - elliptic curve

        -  ``z`` - point on `e`

        -  ``n`` - integer, or a complex quadratic integer of complex
           multiplication for `e`. Complex multiplication currently
           only works if `e` is defined over `Q`.


        EXAMPLES: We consider a curve with CM by `Z[i]`::

            sage: e = pari([0,0,0,3,0]).ellinit()
            sage: p = [1,2]  # Point of infinite order

        Multiplication by two::

            sage: e.ellpow([0,0], 2)
            [0]
            sage: e.ellpow(p, 2)
            [1/4, -7/8]

        Complex multiplication::

            sage: q = e.ellpow(p, 1+I); q
            [-2*I, 1 + I]
            sage: e.ellpow(q, 1-I)
            [1/4, -7/8]

        TESTS::

            sage: for D in [-7, -8, -11, -12, -16, -19, -27, -28]:  # long time (1s)
            ....:     hcpol = hilbert_class_polynomial(D)
            ....:     j = hcpol.roots(multiplicities=False)[0]
            ....:     t = (1728-j)/(27*j)
            ....:     E = EllipticCurve([4*t,16*t^2])
            ....:     P = E.point([0, 4*t])
            ....:     assert(E.j_invariant() == j)
            ....:     #
            ....:     # Compute some CM number and its minimal polynomial
            ....:     #
            ....:     cm = pari('cm = (3*quadgen(%s)+2)'%D)
            ....:     cm_minpoly = pari('minpoly(cm)')
            ....:     #
            ....:     # Evaluate cm_minpoly(cm)(P), which should be zero
            ....:     #
            ....:     e = pari(E)  # Convert E to PARI
            ....:     P2 = e.ellpow(P, cm_minpoly[2]*cm + cm_minpoly[1])
            ....:     P0 = e.elladd(e.ellpow(P, cm_minpoly[0]), e.ellpow(P2, cm))
            ....:     assert(P0 == E(0))
        """
        t0GEN(z); t1GEN(n)
        pari_catch_sig_on()
        return self.new_gen(powell(self.g, t0, t1))

    def ellrootno(self, p=1):
        """
        e.ellrootno(p): return the (local or global) root number of the
        `L`-series of the elliptic curve e

        If p is a prime number, the local root number at p is returned. If
        p is 1, the global root number is returned. Note that the global
        root number is the sign of the functional equation of the
        `L`-series, and therefore conjecturally equal to the parity
        of the rank of e.

        INPUT:


        -  ``e`` - elliptic curve over `\QQ`

        -  ``p (default = 1)`` - 1 or a prime number


        OUTPUT: 1 or -1

        EXAMPLES: Here is a curve of rank 3::

            sage: e = pari([0,0,0,-82,0]).ellinit()
            sage: e.ellrootno()
            -1
            sage: e.ellrootno(2)
            1
            sage: e.ellrootno(1009)
            1
        """
        cdef long rootno
        t0GEN(p)
        pari_catch_sig_on()
        rootno =  ellrootno(self.g, t0)
        pari_catch_sig_off()
        return rootno

    def ellsigma(self, z, flag=0):
        """
        e.ellsigma(z, flag=0): return the value at the complex point z of
        the Weierstrass `\sigma` function associated to the
        elliptic curve e.

        EXAMPLES::

            sage: e = pari([0,0,0,1,0]).ellinit()
            sage: C.<i> = ComplexField()
            sage: e.ellsigma(2+i)
            1.43490215804166 + 1.80307856719256*I
        """
        t0GEN(z)
        pari_catch_sig_on()
        # the prec argument has no effect
        return self.new_gen(ellsigma(self.g, t0, flag, prec))

    def ellsub(self, z0, z1):
        """
        e.ellsub(z0, z1): return z0-z1 on this elliptic curve.

        INPUT:


        -  ``e`` - elliptic curve E

        -  ``z0`` - point on E

        -  ``z1`` - point on E


        OUTPUT: point on E

        EXAMPLES::

            sage: e = pari([0, 1, 1, -2, 0]).ellinit()
            sage: e.ellsub([1,0], [-1,1])
            [0, 0]
        """
        t0GEN(z0); t1GEN(z1)
        pari_catch_sig_on()
        return self.new_gen(subell(self.g, t0, t1))

    def elltaniyama(self):
        pari_catch_sig_on()
        return self.new_gen(taniyama(self.g))

    def elltors(self, flag=0):
        """
        e.elltors(flag = 0): return information about the torsion subgroup
        of the elliptic curve e

        INPUT:


        -  ``e`` - elliptic curve over `\QQ`

        -  ``flag (optional)`` - specify which algorithm to
           use:

        -  ``0 (default)`` - use Doud's algorithm: bound
           torsion by computing the cardinality of e(GF(p)) for small primes
           of good reduction, then look for torsion points using Weierstrass
           parametrization and Mazur's classification

        -  ``1`` - use algorithm given by the Nagell-Lutz
           theorem (this is much slower)


        OUTPUT:


        -  ``gen`` - the order of the torsion subgroup, a.k.a.
           the number of points of finite order

        -  ``gen`` - vector giving the structure of the torsion
           subgroup as a product of cyclic groups, sorted in non-increasing
           order

        -  ``gen`` - vector giving points on e generating these
           cyclic groups


        EXAMPLES::

            sage: e = pari([1,0,1,-19,26]).ellinit()
            sage: e.elltors()
            [12, [6, 2], [[-2, 8], [3, -2]]]
        """
        pari_catch_sig_on()
        return self.new_gen(elltors0(self.g, flag))

    def ellzeta(self, z):
        """
        e.ellzeta(z): return the value at the complex point z of the
        Weierstrass `\zeta` function associated with the elliptic
        curve e.

        .. note::

           This function has infinitely many poles (one of which is at
           z=0); attempting to evaluate it too close to one of the
           poles will result in a PariError.

        INPUT:


        -  ``e`` - elliptic curve

        -  ``z`` - complex number


        EXAMPLES::

            sage: e = pari([0,0,0,1,0]).ellinit()
            sage: e.ellzeta(1)
            1.06479841295883 + 0.E-19*I                # 32-bit
            1.06479841295883 + 5.42101086242752 E-20*I # 64-bit
            sage: C.<i> = ComplexField()
            sage: e.ellzeta(i-1)
            -0.350122658523049 - 0.350122658523049*I
        """
        t0GEN(z)
        pari_catch_sig_on()
        # the prec argument has no effect
        return self.new_gen(ellzeta(self.g, t0, prec))

    def ellztopoint(self, z):
        """
        e.ellztopoint(z): return the point on the elliptic curve e
        corresponding to the complex number z, under the usual complex
        uniformization of e by the Weierstrass p-function.

        INPUT:


        -  ``e`` - elliptic curve

        -  ``z`` - complex number


        OUTPUT point on e

        EXAMPLES::

            sage: e = pari([0,0,0,1,0]).ellinit()
            sage: C.<i> = ComplexField()
            sage: e.ellztopoint(1+i)
            [0.E-19 - 1.02152286795670*I, -0.149072813701096 - 0.149072813701096*I] # 32-bit
            [...  - 1.02152286795670*I, -0.149072813701096 - 0.149072813701096*I] # 64-bit

        Complex numbers belonging to the period lattice of e are of course
        sent to the point at infinity on e::

            sage: e.ellztopoint(0)
            [0]
        """
        t0GEN(z)
        try:
            dprec = prec_words_to_dec(z.precision())
        except AttributeError:
            dprec = prec
        pari_catch_sig_on()
        # the prec argument has no effect
        return self.new_gen(pointell(self.g, t0, dprec))

    def omega(self):
        """
        e.omega(): return basis for the period lattice of the elliptic
        curve e.

        EXAMPLES::

            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.omega()
            [1.26920930427955, -0.634604652139777 - 1.45881661693850*I]
        """
        return self[14:16]

    def disc(self):
        """
        e.disc(): return the discriminant of the elliptic curve e.

        EXAMPLES::

            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.disc()
            -161051
            sage: _.factor()
            [-1, 1; 11, 5]
        """
        return self[11]

    def j(self):
        """
        e.j(): return the j-invariant of the elliptic curve e.

        EXAMPLES::

            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.j()
            -122023936/161051
            sage: _.factor()
            [-1, 1; 2, 12; 11, -5; 31, 3]
        """
        return self[12]

    def ellj(self):
        try:
            dprec = prec_words_to_dec(self.precision())
        except AttributeError:
            dprec = prec
        pari_catch_sig_on()
        return P.new_gen(jell(self.g, dprec))


    ###########################################
    # 6: Functions related to NUMBER FIELDS
    ###########################################
    def bnfcertify(self):
        r"""
        ``bnf`` being as output by ``bnfinit``, checks whether the result is
        correct, i.e. whether the calculation of the contents of ``self``
        are correct without assuming the Generalized Riemann Hypothesis.
        If it is correct, the answer is 1. If not, the program may output
        some error message or loop indefinitely.

        For more information about PARI and the Generalized Riemann
        Hypothesis, see [PariUsers], page 120.

        REFERENCES:

        .. [PariUsers] User's Guide to PARI/GP,
           http://pari.math.u-bordeaux.fr/pub/pari/manuals/2.5.1/users.pdf
        """
        cdef long n
        pari_catch_sig_on()
        n = bnfcertify(self.g)
        pari_catch_sig_off()
        return n

    def bnfinit(self, long flag=0, tech=None):
        if tech is None:
            pari_catch_sig_on()
            return P.new_gen(bnfinit0(self.g, flag, <GEN>0, prec))
        else:
            t0GEN(tech)
            pari_catch_sig_on()
            return P.new_gen(bnfinit0(self.g, flag, t0, prec))

    def bnfisintnorm(self, x):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(bnfisintnorm(self.g, t0))

    def bnfisnorm(self, x, long flag=0):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(bnfisnorm(self.g, t0, flag))

    def bnfisprincipal(self, x, long flag=1):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(bnfisprincipal0(self.g, t0, flag))

    def bnfnarrow(self):
        pari_catch_sig_on()
        return self.new_gen(buchnarrow(self.g))

    def bnfsunit(bnf, S, long precision=0):
        t0GEN(S)
        pari_catch_sig_on()
        return bnf.new_gen(bnfsunit(bnf.g, t0, pbw(precision)))

    def bnfunit(self):
        pari_catch_sig_on()
        return self.new_gen(bnf_get_fu(self.g))

    def bnfisunit(self, x):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(bnfisunit(self.g, t0))

    def dirzetak(self, n):
        t0GEN(n)
        pari_catch_sig_on()
        return self.new_gen(dirzetak(self.g, t0))

    def galoisapply(self, aut, x):
        t0GEN(aut)
        t1GEN(x)
        pari_catch_sig_on()
        return self.new_gen(galoisapply(self.g, t0, t1))

    def galoisinit(self, den=None):
        """
        galoisinit(K{,den}): calculate Galois group of number field K; see PARI manual
        for meaning of den
        """
        if den is None:
            pari_catch_sig_on()
            return self.new_gen(galoisinit(self.g, NULL))
        else:
            t0GEN(den)
            pari_catch_sig_on()
            return self.new_gen(galoisinit(self.g, t0))

    def galoispermtopol(self, perm):
        t0GEN(perm)
        pari_catch_sig_on()
        return self.new_gen(galoispermtopol(self.g, t0))

    def galoisfixedfield(self, perm, long flag=0, v=-1):
        t0GEN(perm);
        pari_catch_sig_on()
        return self.new_gen(galoisfixedfield(self.g, t0, flag, P.get_var(v)))

    def idealred(self, I, vdir=0):
        t0GEN(I); t1GEN(vdir)
        pari_catch_sig_on()
        return self.new_gen(idealred0(self.g, t0, t1 if vdir else NULL))

    def idealadd(self, x, y):
        t0GEN(x); t1GEN(y)
        pari_catch_sig_on()
        return self.new_gen(idealadd(self.g, t0, t1))

    def idealaddtoone(self, x, y):
        t0GEN(x); t1GEN(y)
        pari_catch_sig_on()
        return self.new_gen(idealaddtoone0(self.g, t0, t1))

    def idealappr(self, x, long flag=0):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(idealappr(self.g, t0))

    def idealcoprime(self, x, y):
        """
        Given two integral ideals x and y of a pari number field self,
        return an element a of the field (expressed in the integral
        basis of self) such that a*x is an integral ideal coprime to
        y.

        EXAMPLES::

            sage: F = NumberField(x^3-2, 'alpha')
            sage: nf = F._pari_()
            sage: x = pari('[1, -1, 2]~')
            sage: y = pari('[1, -1, 3]~')
            sage: nf.idealcoprime(x, y)
            [1, 0, 0]~

            sage: y = pari('[2, -2, 4]~')
            sage: nf.idealcoprime(x, y)
            [5/43, 9/43, -1/43]~
        """
        t0GEN(x); t1GEN(y)
        pari_catch_sig_on()
        return self.new_gen(idealcoprime(self.g, t0, t1))

    def idealdiv(self, x, y, long flag=0):
        t0GEN(x); t1GEN(y)
        pari_catch_sig_on()
        return self.new_gen(idealdiv0(self.g, t0, t1, flag))

    def idealfactor(self, x):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(idealfactor(self.g, t0))

    def idealhnf(self, a, b=None):
        t0GEN(a)
        if b is None:
            pari_catch_sig_on()
            return self.new_gen(idealhnf(self.g, t0))
        else:
            t1GEN(b)
            pari_catch_sig_on()
            return self.new_gen(idealhnf0(self.g, t0, t1))

    def idealintersection(self, x, y):
        t0GEN(x); t1GEN(y)
        pari_catch_sig_on()
        return self.new_gen(idealintersect(self.g, t0, t1))

    def ideallist(self, long bound, long flag = 4):
        """
        Vector of vectors `L` of all idealstar of all ideals of `norm <= bound`.

        The binary digits of flag mean:

         - 1: give generators;
         - 2: add units;
         - 4: (default) give only the ideals and not the bid.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^2 + 1)
            sage: L = K.pari_nf().ideallist(100)

        Now we have our list `L`. Entry `L[n-1]` contains all ideals of
        norm `n`::

            sage: L[0]   # One ideal of norm 1.
            [[1, 0; 0, 1]]
            sage: L[64]  # 4 ideals of norm 65.
            [[65, 8; 0, 1], [65, 47; 0, 1], [65, 18; 0, 1], [65, 57; 0, 1]]
        """
        pari_catch_sig_on()
        return self.new_gen(ideallist0(self.g, bound, flag))

    def ideallog(self, x, bid):
        """
        Return the discrete logarithm of the unit x in (ring of integers)/bid.

        INPUT:

        - ``self`` - a pari number field

        - ``bid``  - a big ideal structure (corresponding to an ideal I
          of self) output by idealstar

        - ``x``  - an element of self with valuation zero at all
          primes dividing I

        OUTPUT:

        - the discrete logarithm of x on the generators given in bid[2]

        EXAMPLE::

            sage: F = NumberField(x^3-2, 'alpha')
            sage: nf = F._pari_()
            sage: I = pari('[1, -1, 2]~')
            sage: bid = nf.idealstar(I)
            sage: x = pari('5')
            sage: nf.ideallog(x, bid)
            [25]~
        """
        t0GEN(x); t1GEN(bid)
        pari_catch_sig_on()
        return self.new_gen(ideallog(self.g, t0, t1))

    def idealmul(self, x, y, long flag=0):
        t0GEN(x); t1GEN(y)
        pari_catch_sig_on()
        if flag == 0:
            return self.new_gen(idealmul(self.g, t0, t1))
        else:
            return self.new_gen(idealmulred(self.g, t0, t1))

    def idealnorm(self, x):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(idealnorm(self.g, t0))

    def idealprimedec(nf, p):
        """
        Prime ideal decomposition of the prime number `p` in the number
        field `nf` as a vector of 5 component vectors `[p,a,e,f,b]`
        representing the prime ideals `p O_K + a O_K`, `e` ,`f` as usual,
        `a` as vector of components on the integral basis, `b` Lenstra's
        constant.

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: F = pari(K).idealprimedec(5); F
            [[5, [-2, 1]~, 1, 1, [2, 1]~], [5, [2, 1]~, 1, 1, [-2, 1]~]]
            sage: F[0].pr_get_p()
            5
        """
        t0GEN(p)
        pari_catch_sig_on()
        return nf.new_gen(idealprimedec(nf.g, t0))

    def idealstar(self, I, long flag=1):
        """
        Return the big ideal (bid) structure of modulus I.

        INPUT:

        - ``self`` - a pari number field

        - ``I`` -- an ideal of self, or a row vector whose first
          component is an ideal and whose second component
          is a row vector of r_1 0 or 1.

        - ``flag`` - determines the amount of computation and the shape
          of the output:

          - ``1`` (default): return a bid structure without
            generators

          - ``2``: return a bid structure with generators (slower)

          - ``0`` (deprecated): only outputs units of (ring of integers/I)
            as an abelian group, i.e as a 3-component
            vector [h,d,g]: h is the order, d is the vector
            of SNF cyclic components and g the corresponding
            generators. This flag is deprecated: it is in
            fact slightly faster to compute a true bid
            structure, which contains much more information.

        EXAMPLE::

            sage: F = NumberField(x^3-2, 'alpha')
            sage: nf = F._pari_()
            sage: I = pari('[1, -1, 2]~')
            sage: nf.idealstar(I)
            [[[43, 9, 5; 0, 1, 0; 0, 0, 1], [0]], [42, [42]], Mat([[43, [9, 1, 0]~, 1, 1, [-5, -9, 1]~], 1]), [[[[42], [[3, 0, 0]~], [[3, 0, 0]~], [Vecsmall([])], 1]], [[], [], []]], Mat(1)]
        """
        t0GEN(I)
        pari_catch_sig_on()
        return self.new_gen(idealstar0(self.g, t0, flag))

    def idealtwoelt(self, x, a=None):
        t0GEN(x)
        if a is None:
            pari_catch_sig_on()
            return self.new_gen(idealtwoelt0(self.g, t0, NULL))
        else:
            t1GEN(a)
            pari_catch_sig_on()
            return self.new_gen(idealtwoelt0(self.g, t0, t1))

    def idealval(self, x, p):
        cdef long v
        t0GEN(x); t1GEN(p)
        pari_catch_sig_on()
        v = idealval(self.g, t0, t1)
        pari_catch_sig_off()
        return v

    def elementval(self, x, p):
        cdef long v
        t0GEN(x); t1GEN(p)
        pari_catch_sig_on()
        v = nfval(self.g, t0, t1)
        pari_catch_sig_off()
        return v

    def modreverse(self):
        """
        modreverse(x): reverse polymod of the polymod x, if it exists.

        EXAMPLES:
        """
        pari_catch_sig_on()
        return self.new_gen(modreverse(self.g))

    def nfbasis(self, long flag=0, fa=0):
        """
        nfbasis(x, flag, fa): integral basis of the field QQ[a], where ``a`` is
        a root of the polynomial x.

        Binary digits of ``flag`` mean:

         - 1: assume that no square of a prime>primelimit divides the
              discriminant of ``x``.
         - 2: use round 2 algorithm instead of round 4.

        If present, ``fa`` provides the matrix of a partial factorization of
        the discriminant of ``x``, useful if one wants only an order maximal at
        certain primes only.

        EXAMPLES::

            sage: pari('x^3 - 17').nfbasis()
            [1, x, 1/3*x^2 - 1/3*x + 1/3]

        We test ``flag`` = 1, noting it gives a wrong result when the
        discriminant (-4 * `p`^2 * `q` in the example below) has a big square
        factor::

            sage: p = next_prime(10^10); q = next_prime(p)
            sage: x = polygen(QQ); f = x^2 + p^2*q
            sage: pari(f).nfbasis(1)   # Wrong result
            [1, x]
            sage: pari(f).nfbasis()    # Correct result
            [1, 1/10000000019*x]
            sage: pari(f).nfbasis(fa = "[2,2; %s,2]"%p)    # Correct result and faster
            [1, 1/10000000019*x]

        TESTS:

        ``flag`` = 2 should give the same result::

            sage: pari('x^3 - 17').nfbasis(flag = 2)
            [1, x, 1/3*x^2 - 1/3*x + 1/3]
        """
        global t0
        t0GEN(fa)
        if typ(t0) != t_MAT:
            t0 = <GEN>0
        pari_catch_sig_on()
        return self.new_gen(nfbasis0(self.g, flag, t0))

    def nfbasis_d(self, long flag=0, fa=0):
        """
        nfbasis_d(x): Return a basis of the number field defined over QQ
        by x and its discriminant.

        EXAMPLES::

            sage: F = NumberField(x^3-2,'alpha')
            sage: F._pari_()[0].nfbasis_d()
            ([1, y, y^2], -108)

        ::

            sage: G = NumberField(x^5-11,'beta')
            sage: G._pari_()[0].nfbasis_d()
            ([1, y, y^2, y^3, y^4], 45753125)

        ::

            sage: pari([-2,0,0,1]).Polrev().nfbasis_d()
            ([1, x, x^2], -108)
        """
        global t0
        cdef GEN disc
        t0GEN(fa)
        if typ(t0) != t_MAT:
            t0 = <GEN>0
        pari_catch_sig_on()
        B = self.new_gen_noclear(nfbasis(self.g, &disc, flag, t0))
        D = self.new_gen(disc);
        return B,D

    def nfbasistoalg(nf, x):
        r"""
        Transforms the column vector ``x`` on the integral basis into an
        algebraic number.

        INPUT:

         - ``nf`` -- a number field
         - ``x`` -- a column of rational numbers of length equal to the
           degree of ``nf`` or a single rational number

        OUTPUT:

         - A POLMOD representing the element of ``nf`` whose coordinates
           are ``x`` in the Z-basis of ``nf``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^3 - 17)
            sage: Kpari = K.pari_nf()
            sage: Kpari.getattr('zk')
            [1, 1/3*y^2 - 1/3*y + 1/3, y]
            sage: Kpari.nfbasistoalg(42)
            Mod(42, y^3 - 17)
            sage: Kpari.nfbasistoalg("[3/2, -5, 0]~")
            Mod(-5/3*y^2 + 5/3*y - 1/6, y^3 - 17)
            sage: Kpari.getattr('zk') * pari("[3/2, -5, 0]~")
            -5/3*y^2 + 5/3*y - 1/6
        """
        t0GEN(x)
        pari_catch_sig_on()
        return nf.new_gen(basistoalg(nf.g, t0))

    def nfbasistoalg_lift(nf, x):
        r"""
        Transforms the column vector ``x`` on the integral basis into a
        polynomial representing the algebraic number.

        INPUT:

         - ``nf`` -- a number field
         - ``x`` -- a column of rational numbers of length equal to the
           degree of ``nf`` or a single rational number

        OUTPUT:

         - ``nf.nfbasistoalg(x).lift()``

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^3 - 17)
            sage: Kpari = K.pari_nf()
            sage: Kpari.getattr('zk')
            [1, 1/3*y^2 - 1/3*y + 1/3, y]
            sage: Kpari.nfbasistoalg_lift(42)
            42
            sage: Kpari.nfbasistoalg_lift("[3/2, -5, 0]~")
            -5/3*y^2 + 5/3*y - 1/6
            sage: Kpari.getattr('zk') * pari("[3/2, -5, 0]~")
            -5/3*y^2 + 5/3*y - 1/6
        """
        t0GEN(x)
        pari_catch_sig_on()
        return nf.new_gen(gel(basistoalg(nf.g, t0), 2))

    def nfdisc(self, long flag=0, p=0):
        """
        nfdisc(x): Return the discriminant of the number field defined over
        QQ by x.

        EXAMPLES::

            sage: F = NumberField(x^3-2,'alpha')
            sage: F._pari_()[0].nfdisc()
            -108

        ::

            sage: G = NumberField(x^5-11,'beta')
            sage: G._pari_()[0].nfdisc()
            45753125

        ::

            sage: f = x^3-2
            sage: f._pari_()
            x^3 - 2
            sage: f._pari_().nfdisc()
            -108
        """
        cdef gen _p
        cdef GEN g
        if p != 0:
            _p = self.pari(p)
            g = _p.g
        else:
            g = <GEN>NULL
        pari_catch_sig_on()
        return self.new_gen(nfdisc0(self.g, flag, g))

    def nfeltdiveuc(self, x, y):
        """
        Given `x` and `y` in the number field ``self``, return `q` such
        that `x - q y` is "small".

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 5)
            sage: x = 10
            sage: y = a + 1
            sage: pari(k).nfeltdiveuc(pari(x), pari(y))
            [2, -2]~
        """
        t0GEN(x); t1GEN(y)
        pari_catch_sig_on()
        return self.new_gen(nfdiveuc(self.g, t0, t1))

    def nfeltreduce(self, x, I):
        """
        Given an ideal I in Hermite normal form and an element x of the pari
        number field self, finds an element r in self such that x-r belongs
        to the ideal and r is small.

        EXAMPLES::

            sage: k.<a> = NumberField(x^2 + 5)
            sage: I = k.ideal(a)
            sage: kp = pari(k)
            sage: kp.nfeltreduce(12, I.pari_hnf())
            [2, 0]~
            sage: 12 - k(kp.nfeltreduce(12, I.pari_hnf())) in I
            True
        """
        t0GEN(x); t1GEN(I)
        pari_catch_sig_on()
        return self.new_gen(nfreduce(self.g, t0, t1))

    def nffactor(self, x):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(nffactor(self.g, t0))

    def nfgenerator(self):
        f = self[0]
        x = f.variable()
        return x.Mod(f)

    def nfhilbert(self, a, b, p = None):
        """
        nfhilbert(nf,a,b,{p}): if p is omitted, global Hilbert symbol (a,b)
        in nf, that is 1 if X^2-aY^2-bZ^2 has a non-trivial solution (X,Y,Z)
        in nf, -1 otherwise. Otherwise compute the local symbol modulo the
        prime ideal p.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<t> = NumberField(x^3 - x + 1)
            sage: pari(K).nfhilbert(t, t+2)  # not tested, known bug #11868
            -1
            sage: pari(K).nfhilbert(pari(t), pari(t+2))
            -1
            sage: P = K.ideal(t^2 + t - 2)   # Prime ideal above 5
            sage: pari(K).nfhilbert(pari(t), pari(t+2), P.pari_prime())
            -1
            sage: P = K.ideal(t^2 + 3*t - 1) # Prime ideal above 23, ramified
            sage: pari(K).nfhilbert(pari(t), pari(t+2), P.pari_prime())
            1
        """
        cdef long r
        t0GEN(a)
        t1GEN(b)
        if p:
            t2GEN(p)
            pari_catch_sig_on()
            r = nfhilbert0(self.g, t0, t1, t2)
        else:
            pari_catch_sig_on()
            r = nfhilbert(self.g, t0, t1)
        P.clear_stack()
        return r


    def nfhnf(self,x):
        """
        nfhnf(nf,x) : given a pseudo-matrix (A, I) or an integral pseudo-matrix (A,I,J), finds a
        pseudo-basis in Hermite normal form of the module it generates.

        A pseudo-matrix is a 2-component row vector (A, I) where A is a relative m x n matrix and
        I an ideal list of length n. An integral pseudo-matrix is a 3-component row vector (A, I, J).

        .. NOTE::

            The definition of a pseudo-basis ([Cohen]_):
            Let M be a finitely generated, torsion-free R-module, and set V = KM.  If `\mathfrak{a}_i` are
            fractional ideals of R and `w_i` are elements of V, we say that
            `(w_i, \mathfrak{a}_k)_{1 \leq i \leq k}`
            is a pseudo-basis of M if
            `M = \mathfrak{a}_1 w_1 \oplus \cdots \oplus \mathfrak{a}_k w_k.`

        REFERENCES:

        .. [Cohen] Cohen, "Advanced Topics in Computational Number Theory"

        EXAMPLES::

            sage: F.<a> = NumberField(x^2-x-1)
            sage: Fp = pari(F)
            sage: A = matrix(F,[[1,2,a,3],[3,0,a+2,0],[0,0,a,2],[3+a,a,0,1]])
            sage: I = [F.ideal(-2*a+1),F.ideal(7), F.ideal(3),F.ideal(1)]
            sage: Fp.nfhnf([pari(A),[pari(P) for P in I]])
            [[1, [-969/5, -1/15]~, [15, -2]~, [-1938, -3]~; 0, 1, 0, 0; 0, 0, 1, 0;
            0, 0, 0, 1], [[3997, 1911; 0, 7], [15, 6; 0, 3], [1, 0; 0, 1], [1, 0; 0,
            1]]]
            sage: K.<b> = NumberField(x^3-2)
            sage: Kp = pari(K)
            sage: A = matrix(K,[[1,0,0,5*b],[1,2*b^2,b,57],[0,2,1,b^2-3],[2,0,0,b]])
            sage: I = [K.ideal(2),K.ideal(3+b^2),K.ideal(1),K.ideal(1)]
            sage: Kp.nfhnf([pari(A),[pari(P) for P in I]])
            [[1, -225, 72, -31; 0, 1, [0, -1, 0]~, [0, 0, -1/2]~; 0, 0, 1, [0, 0,
            -1/2]~; 0, 0, 0, 1], [[1116, 756, 612; 0, 18, 0; 0, 0, 18], [2, 0, 0; 0,
            2, 0; 0, 0, 2], [1, 0, 0; 0, 1, 0; 0, 0, 1], [2, 0, 0; 0, 1, 0; 0, 0,
            1]]]

        An example where the ring of integers of the number field is not a PID::

            sage: K.<b> = NumberField(x^2+5)
            sage: Kp = pari(K)
            sage: A = matrix(K,[[1,0,0,5*b],[1,2*b^2,b,57],[0,2,1,b^2-3],[2,0,0,b]])
            sage: I = [K.ideal(2),K.ideal(3+b^2),K.ideal(1),K.ideal(1)]
            sage: Kp.nfhnf([pari(A),[pari(P) for P in I]])
            [[1, [15, 6]~, [0, -54]~, [113, 72]~; 0, 1, [-4, -1]~, [0, -1]~; 0, 0,
            1, 0; 0, 0, 0, 1], [[360, 180; 0, 180], [6, 4; 0, 2], [1, 0; 0, 1], [1,
            0; 0, 1]]]
            sage: A = matrix(K,[[1,0,0,5*b],[1,2*b,b,57],[0,2,1,b-3],[2,0,b,b]])
            sage: I = [K.ideal(2).factor()[0][0],K.ideal(3+b),K.ideal(1),K.ideal(1)]
            sage: Kp.nfhnf([pari(A),[pari(P) for P in I]])
            [[1, [7605, 4]~, [5610, 5]~, [7913, -6]~; 0, 1, 0, -1; 0, 0, 1, 0; 0, 0,
            0, 1], [[19320, 13720; 0, 56], [2, 1; 0, 1], [1, 0; 0, 1], [1, 0; 0,
            1]]]

        AUTHORS:

        - Aly Deines (2012-09-19)
        """
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(nfhnf(self.g,t0))




    def nfinit(self, long flag=0, long precision=0):
        """
        nfinit(pol, {flag=0}): ``pol`` being a nonconstant irreducible
        polynomial, gives a vector containing all the data necessary for PARI
        to compute in this number field.

        ``flag`` is optional and can be set to:
         - 0: default
         - 1: do not compute different
         - 2: first use polred to find a simpler polynomial
         - 3: outputs a two-element vector [nf,Mod(a,P)], where nf is as in 2
              and Mod(a,P) is a polmod equal to Mod(x,pol) and P=nf.pol

        EXAMPLES::

            sage: pari('x^3 - 17').nfinit()
            [x^3 - 17, [1, 1], -867, 3, [[1, 1.68006..., 2.57128...; 1, -0.340034... + 2.65083...*I, -1.28564... - 2.22679...*I], [1, 1.68006..., 2.57128...; 1, 2.31080..., -3.51243...; 1, -2.99087..., 0.941154...], [1, 2, 3; 1, 2, -4; 1, -3, 1], [3, 1, 0; 1, -11, 17; 0, 17, 0], [51, 0, 16; 0, 17, 3; 0, 0, 1], [17, 0, -1; 0, 0, 3; -1, 3, 2], [51, [-17, 6, -1; 0, -18, 3; 1, 0, -16]]], [2.57128..., -1.28564... - 2.22679...*I], [1, 1/3*x^2 - 1/3*x + 1/3, x], [1, 0, -1; 0, 0, 3; 0, 1, 1], [1, 0, 0, 0, -4, 6, 0, 6, -1; 0, 1, 0, 1, 1, -1, 0, -1, 3; 0, 0, 1, 0, 2, 0, 1, 0, 1]]

        TESTS:

        This example only works after increasing precision::

            sage: pari('x^2 + 10^100 + 1').nfinit(precision=64)
            Traceback (most recent call last):
            ...
            PariError: precision too low (10)
            sage: pari('x^2 + 10^100 + 1').nfinit()
            [...]

        Throw a PARI error which is not precer::

            sage: pari('1.0').nfinit()
            Traceback (most recent call last):
            ...
            PariError: incorrect type (11)
        """

        # If explicit precision is given, use only that
        if precision:
            return self._nfinit_with_prec(flag, precision)

        # Otherwise, start with 64 bits of precision and increase as needed:
        precision = 64
        while True:
            try:
                return self._nfinit_with_prec(flag, precision)
            except PariError, err:
                if err.errnum() == precer:
                    precision *= 2
                else:
                    raise err

    # NOTE: because of the way pari_catch_sig_on() and Cython exceptions work, this
    # function MUST NOT be folded into nfinit() above. It has to be a
    # seperate function.
    def _nfinit_with_prec(self, long flag, long precision):
        """
        See ``self.nfinit()``.
        """
        pari_catch_sig_on()
        return P.new_gen(nfinit0(self.g, flag, pbw(precision)))

    def nfisisom(self, gen other):
        """
        nfisisom(x, y): Determine if the number fields defined by x and y
        are isomorphic. According to the PARI documentation, this is much
        faster if at least one of x or y is a number field. If they are
        isomorphic, it returns an embedding for the generators. If not,
        returns 0.

        EXAMPLES::

            sage: F = NumberField(x^3-2,'alpha')
            sage: G = NumberField(x^3-2,'beta')
            sage: F._pari_().nfisisom(G._pari_())
            [y]

        ::

            sage: GG = NumberField(x^3-4,'gamma')
            sage: F._pari_().nfisisom(GG._pari_())
            [1/2*y^2]

        ::

            sage: F._pari_().nfisisom(GG.pari_nf())
            [1/2*y^2]

        ::

            sage: F.pari_nf().nfisisom(GG._pari_()[0])
            [y^2]

        ::

            sage: H = NumberField(x^2-2,'alpha')
            sage: F._pari_().nfisisom(H._pari_())
            0
        """
        pari_catch_sig_on()
        return P.new_gen(nfisisom(self.g, other.g))

    def nfrootsof1(self):
        """
        nf.nfrootsof1()

        number of roots of unity and primitive root of unity in the number
        field nf.

        EXAMPLES::

            sage: nf = pari('x^2 + 1').nfinit()
            sage: nf.nfrootsof1()
            [4, -x]
        """
        pari_catch_sig_on()
        return P.new_gen(rootsof1(self.g))

    def nfsubfields(self, long d=0):
        """
        Find all subfields of degree d of number field nf (all subfields if
        d is null or omitted). Result is a vector of subfields, each being
        given by [g,h], where g is an absolute equation and h expresses one
        of the roots of g in terms of the root x of the polynomial defining
        nf.

        INPUT:


        -  ``self`` - nf number field

        -  ``d`` - C long integer
        """
        pari_catch_sig_on()
        return self.new_gen(nfsubfields(self.g, d))

    def rnfcharpoly(self, T, a, v='x'):
        t0GEN(T); t1GEN(a); t2GEN(v)
        pari_catch_sig_on()
        return self.new_gen(rnfcharpoly(self.g, t0, t1, gvar(t2)))

    def rnfdisc(self, x):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(rnfdiscf(self.g, t0))

    def rnfeltabstorel(self, x):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(rnfelementabstorel(self.g, t0))

    def rnfeltreltoabs(self, x):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(rnfelementreltoabs(self.g, t0))

    def rnfequation(self, poly, long flag=0):
        t0GEN(poly)
        pari_catch_sig_on()
        return self.new_gen(rnfequation0(self.g, t0, flag))

    def rnfidealabstorel(self, x):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(rnfidealabstorel(self.g, t0))

    def rnfidealdown(self, x):
        r"""
        rnfidealdown(rnf,x): finds the intersection of the ideal x with the base field.

        EXAMPLES:
            sage: x = ZZ['xx1'].0; pari(x)
            xx1
            sage: y = ZZ['yy1'].0; pari(y)
            yy1
            sage: nf = pari(y^2 - 6*y + 24).nfinit()
            sage: rnf = nf.rnfinit(x^2 - pari(y))

        This is the relative HNF of the inert ideal (2) in rnf:

            sage: P = pari('[[[1, 0]~, [0, 0]~; [0, 0]~, [1, 0]~], [[2, 0; 0, 2], [2, 0; 0, 1/2]]]')

        And this is the HNF of the inert ideal (2) in nf:

            sage: rnf.rnfidealdown(P)
            [2, 0; 0, 2]
        """
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(rnfidealdown(self.g, t0))

    def rnfidealhnf(self, x):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(rnfidealhermite(self.g, t0))

    def rnfidealnormrel(self, x):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(rnfidealnormrel(self.g, t0))

    def rnfidealreltoabs(self, x):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(rnfidealreltoabs(self.g, t0))

    def rnfidealtwoelt(self, x):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(rnfidealtwoelement(self.g, t0))

    def rnfinit(self, poly):
        """
        EXAMPLES: We construct a relative number field.

        ::

            sage: f = pari('y^3+y+1')
            sage: K = f.nfinit()
            sage: x = pari('x'); y = pari('y')
            sage: g = x^5 - x^2 + y
            sage: L = K.rnfinit(g)
        """
        t0GEN(poly)
        pari_catch_sig_on()
        return P.new_gen(rnfinit(self.g, t0))

    def rnfisfree(self, poly):
        t0GEN(poly)
        pari_catch_sig_on()
        r = rnfisfree(self.g, t0)
        pari_catch_sig_off()
        return r

    def quadhilbert(self):
        r"""
        Returns a polynomial over `\QQ` whose roots generate the
        Hilbert class field of the quadratic field of discriminant
        ``self`` (which must be fundamental).

        EXAMPLES::

            sage: pari(-23).quadhilbert()
            x^3 - x^2 + 1
            sage: pari(145).quadhilbert()
            x^4 - x^3 - 3*x^2 + x + 1
            sage: pari(-12).quadhilbert()   # Not fundamental
            Traceback (most recent call last):
            ...
            PariError:  (5)
        """
        pari_catch_sig_on()
        # Precision argument is only used for real quadratic extensions
        # and will be automatically increased by PARI if needed.
        return P.new_gen(quadhilbert(self.g, DEFAULTPREC))


    ##################################################
    # 7: POLYNOMIALS and power series
    ##################################################
    def reverse(self):
        """
        Return the polynomial obtained by reversing the coefficients of
        this polynomial.
        """
        return self.Vec().Polrev()

    def content(self):
        """
        Greatest common divisor of all the components of ``self``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: pari(2*x^2 + 2).content()
            2
            sage: pari("4*x^3 - 2*x/3 + 2/5").content()
            2/15
        """
        pari_catch_sig_on()
        return self.new_gen(content(self.g))

    def deriv(self, v=-1):
        pari_catch_sig_on()
        return self.new_gen(deriv(self.g, self.get_var(v)))

    def eval(self, x):
        t0GEN(x)
        pari_catch_sig_on()
        return self.new_gen(poleval(self.g, t0))

    def __call__(self, x):
        return self.eval(x)

    def factornf(self, t):
        """
        Factorization of the polynomial ``self`` over the number field
        defined by the polynomial ``t``.  This does not require that `t`
        is integral, nor that the discriminant of the number field can be
        factored.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2 - 1/8)
            sage: pari(x^2 - 2).factornf(K.pari_polynomial("a"))
            [x + Mod(-4*a, 8*a^2 - 1), 1; x + Mod(4*a, 8*a^2 - 1), 1]
        """
        t0GEN(t)
        pari_catch_sig_on()
        return self.new_gen(polfnf(self.g, t0))

    def factorpadic(self, p, long r=20, long flag=0):
        """
        self.factorpadic(p,r=20,flag=0): p-adic factorization of the
        polynomial x to precision r. flag is optional and may be set to 0
        (use round 4) or 1 (use Buchmann-Lenstra)
        """
        t0GEN(p)
        pari_catch_sig_on()
        return self.new_gen(factorpadic0(self.g, t0, r, flag))

    def factormod(self, p, long flag=0):
        """
        x.factormod(p,flag=0): factorization mod p of the polynomial x
        using Berlekamp. flag is optional, and can be 0: default or 1:
        simple factormod, same except that only the degrees of the
        irreducible factors are given.
        """
        t0GEN(p)
        pari_catch_sig_on()
        return self.new_gen(factormod0(self.g, t0, flag))

    def intformal(self, y=-1):
        """
        x.intformal(y): formal integration of x with respect to the main
        variable of y, or to the main variable of x if y is omitted
        """
        pari_catch_sig_on()
        return self.new_gen(integ(self.g, self.get_var(y)))

    def padicappr(self, a):
        """
        x.padicappr(a): p-adic roots of the polynomial x congruent to a mod
        p
        """
        t0GEN(a)
        pari_catch_sig_on()
        return self.new_gen(padicappr(self.g, t0))

    def newtonpoly(self, p):
        """
        x.newtonpoly(p): Newton polygon of polynomial x with respect to the
        prime p.

        EXAMPLES::

            sage: x = pari('y^8+6*y^6-27*y^5+1/9*y^2-y+1')
            sage: x.newtonpoly(3)
            [1, 1, -1/3, -1/3, -1/3, -1/3, -1/3, -1/3]
        """
        t0GEN(p)
        pari_catch_sig_on()
        return self.new_gen(newtonpoly(self.g, t0))

    def polcoeff(self, long n, var=-1):
        """
        EXAMPLES::

            sage: f = pari("x^2 + y^3 + x*y")
            sage: f
            x^2 + y*x + y^3
            sage: f.polcoeff(1)
            y
            sage: f.polcoeff(3)
            0
            sage: f.polcoeff(3, "y")
            1
            sage: f.polcoeff(1, "y")
            x
        """
        pari_catch_sig_on()
        return self.new_gen(polcoeff0(self.g, n, self.get_var(var)))

    def polcompositum(self, pol2, long flag=0):
        t0GEN(pol2)
        pari_catch_sig_on()
        return self.new_gen(polcompositum0(self.g, t0, flag))

    def poldegree(self, var=-1):
        """
        f.poldegree(var=x): Return the degree of this polynomial.
        """
        pari_catch_sig_on()
        n = poldegree(self.g, self.get_var(var))
        pari_catch_sig_off()
        return n

    def poldisc(self, var=-1):
        """
        f.poldist(var=x): Return the discriminant of this polynomial.
        """
        pari_catch_sig_on()
        return self.new_gen(poldisc0(self.g, self.get_var(var)))

    def poldiscreduced(self):
        pari_catch_sig_on()
        return self.new_gen(reduceddiscsmith(self.g))

    def polgalois(self):
        """
        f.polgalois(): Galois group of the polynomial f
        """
        pari_catch_sig_on()
        return self.new_gen(polgalois(self.g, prec))

    def nfgaloisconj(self, long flag=0, denom=None, long precision=0):
        r"""
        Edited from the pari documentation:

        nfgaloisconj(nf): list of conjugates of a root of the
        polynomial x=nf.pol in the same number field.

        Uses a combination of Allombert's algorithm and nfroots.

        EXAMPLES::

            sage: x = QQ['x'].0; nf = pari(x^2 + 2).nfinit()
            sage: nf.nfgaloisconj()
            [-x, x]~
            sage: nf = pari(x^3 + 2).nfinit()
            sage: nf.nfgaloisconj()
            [x]~
            sage: nf = pari(x^4 + 2).nfinit()
            sage: nf.nfgaloisconj()
            [-x, x]~
        """
        global t0
        if denom is not None:
            t0GEN(denom)
        else:
            t0 = NULL
        pari_catch_sig_on()
        return self.new_gen(galoisconj0(self.g, flag, t0, pbw(precision)))

    def nfroots(self, poly):
        r"""
        Return the roots of `poly` in the number field self without
        multiplicity.

        EXAMPLES::

            sage: y = QQ['yy'].0; _ = pari(y) # pari has variable ordering rules
            sage: x = QQ['zz'].0; nf = pari(x^2 + 2).nfinit()
            sage: nf.nfroots(y^2 + 2)
            [Mod(-zz, zz^2 + 2), Mod(zz, zz^2 + 2)]
            sage: nf = pari(x^3 + 2).nfinit()
            sage: nf.nfroots(y^3 + 2)
            [Mod(zz, zz^3 + 2)]
            sage: nf = pari(x^4 + 2).nfinit()
            sage: nf.nfroots(y^4 + 2)
            [Mod(-zz, zz^4 + 2), Mod(zz, zz^4 + 2)]
        """
        t0GEN(poly)
        pari_catch_sig_on()
        return self.new_gen(nfroots(self.g, t0))

    def polhensellift(self, y, p, long e):
        """
        self.polhensellift(y, p, e): lift the factorization y of self
        modulo p to a factorization modulo `p^e` using Hensel lift.
        The factors in y must be pairwise relatively prime modulo p.
        """
        t0GEN(y)
        t1GEN(p)
        pari_catch_sig_on()
        return self.new_gen(polhensellift(self.g, t0, t1, e))

    def polisirreducible(self):
        """
        f.polisirreducible(): Returns True if f is an irreducible
        non-constant polynomial, or False if f is reducible or constant.
        """
        pari_catch_sig_on()
        return bool(self.new_gen(gisirreducible(self.g)))


    def pollead(self, v=-1):
        """
        self.pollead(v): leading coefficient of polynomial or series self,
        or self itself if self is a scalar. Error otherwise. With respect
        to the main variable of self if v is omitted, with respect to the
        variable v otherwise
        """
        pari_catch_sig_on()
        return self.new_gen(pollead(self.g, self.get_var(v)))

    def polrecip(self):
        pari_catch_sig_on()
        return self.new_gen(polrecip(self.g))

    def polred(self, flag=0, fa=None):
        if fa is None:
            pari_catch_sig_on()
            return self.new_gen(polred0(self.g, flag, NULL))
        else:
            t0GEN(fa)
            pari_catch_sig_on()
            return self.new_gen(polred0(self.g, flag, t0))

    def polredabs(self, flag=0):
        pari_catch_sig_on()
        return self.new_gen(polredabs0(self.g, flag))

    def polredbest(self, flag=0):
        pari_catch_sig_on()
        return self.new_gen(polredbest(self.g, flag))

    def polresultant(self, y, var=-1, flag=0):
        t0GEN(y)
        pari_catch_sig_on()
        return self.new_gen(polresultant0(self.g, t0, self.get_var(var), flag))

    def polroots(self, flag=0, precision=0):
        """
        polroots(x,flag=0): complex roots of the polynomial x. flag is
        optional, and can be 0: default, uses Schonhage's method modified
        by Gourdon, or 1: uses a modified Newton method.
        """
        pari_catch_sig_on()
        return self.new_gen(roots0(self.g, flag, pbw(precision)))

    def polrootsmod(self, p, flag=0):
        t0GEN(p)
        pari_catch_sig_on()
        return self.new_gen(rootmod0(self.g, t0, flag))

    def polrootspadic(self, p, r=20):
        t0GEN(p)
        pari_catch_sig_on()
        return self.new_gen(rootpadic(self.g, t0, r))

    def polrootspadicfast(self, p, r=20):
        t0GEN(p)
        pari_catch_sig_on()
        return self.new_gen(rootpadicfast(self.g, t0, r))

    def polsturm(self, a, b):
        t0GEN(a)
        t1GEN(b)
        pari_catch_sig_on()
        n = sturmpart(self.g, t0, t1)
        pari_catch_sig_off()
        return n

    def polsturm_full(self):
        pari_catch_sig_on()
        n = sturmpart(self.g, NULL, NULL)
        pari_catch_sig_off()
        return n

    def polsylvestermatrix(self, g):
        t0GEN(g)
        pari_catch_sig_on()
        return self.new_gen(sylvestermatrix(self.g, t0))

    def polsym(self, long n):
        pari_catch_sig_on()
        return self.new_gen(polsym(self.g, n))

    def serconvol(self, g):
        t0GEN(g)
        pari_catch_sig_on()
        return self.new_gen(convol(self.g, t0))

    def serlaplace(self):
        pari_catch_sig_on()
        return self.new_gen(laplace(self.g))

    def serreverse(self):
        """
        serreverse(f): reversion of the power series f.

        If f(t) is a series in t with valuation 1, find the series g(t)
        such that g(f(t)) = t.

        EXAMPLES::

            sage: f = pari('x+x^2+x^3+O(x^4)'); f
            x + x^2 + x^3 + O(x^4)
            sage: g = f.serreverse(); g
            x - x^2 + x^3 + O(x^4)
            sage: f.subst('x',g)
            x + O(x^4)
            sage: g.subst('x',f)
            x + O(x^4)
        """
        pari_catch_sig_on()
        return self.new_gen(recip(self.g))

    def thueinit(self, flag=0):
        pari_catch_sig_on()
        return self.new_gen(thueinit(self.g, flag, prec))


    def rnfisnorminit(self, polrel, long flag=2):
        t0GEN(polrel)
        pari_catch_sig_on()
        return self.new_gen(rnfisnorminit(self.g, t0, flag))

    def rnfisnorm(self, T, long flag=0):
        t0GEN(T)
        pari_catch_sig_on()
        return self.new_gen(rnfisnorm(t0, self.g, flag))

    ###########################################
    # 8: Vectors, matrices, LINEAR ALGEBRA and sets
    ###########################################

    def vecextract(self, y, z=None):
        r"""
        self.vecextract(y,z): extraction of the components of the matrix or
        vector x according to y and z. If z is omitted, y designates
        columns, otherwise y corresponds to rows and z to columns. y and z
        can be vectors (of indices), strings (indicating ranges as
        in"1..10") or masks (integers whose binary representation indicates
        the indices to extract, from left to right 1, 2, 4, 8, etc.)

        .. note::

           This function uses the PARI row and column indexing, so the
           first row or column is indexed by 1 instead of 0.
        """
        t0GEN(y)
        if z is None:
            pari_catch_sig_on()
            return P.new_gen(shallowextract(self.g, t0))
        else:
            t1GEN(z)
            pari_catch_sig_on()
            return P.new_gen(extract0(self.g, t0, t1))

    def ncols(self):
        """
        Return the number of columns of self.

        EXAMPLES::

            sage: pari('matrix(19,8)').ncols()
            8
        """
        cdef long n
        pari_catch_sig_on()
        n = glength(self.g)
        pari_catch_sig_off()
        return n

    def nrows(self):
        """
        Return the number of rows of self.

        EXAMPLES::

            sage: pari('matrix(19,8)').nrows()
            19
        """
        cdef long n
        pari_catch_sig_on()
        # if this matrix has no columns
        # then it has no rows.
        if self.ncols() == 0:
            pari_catch_sig_off()
            return 0
        n = glength(<GEN>(self.g[1]))
        pari_catch_sig_off()
        return n

    def mattranspose(self):
        """
        Transpose of the matrix self.

        EXAMPLES::

            sage: pari('[1,2,3; 4,5,6;  7,8,9]').mattranspose()
            [1, 4, 7; 2, 5, 8; 3, 6, 9]
        """
        pari_catch_sig_on()
        return self.new_gen(gtrans(self.g)).Mat()

    def matadjoint(self):
        """
        matadjoint(x): adjoint matrix of x.

        EXAMPLES::

            sage: pari('[1,2,3; 4,5,6;  7,8,9]').matadjoint()
            [-3, 6, -3; 6, -12, 6; -3, 6, -3]
            sage: pari('[a,b,c; d,e,f; g,h,i]').matadjoint()
            [(i*e - h*f), (-i*b + h*c), (f*b - e*c); (-i*d + g*f), i*a - g*c, -f*a + d*c; (h*d - g*e), -h*a + g*b, e*a - d*b]
        """
        pari_catch_sig_on()
        return self.new_gen(adj(self.g)).Mat()

    def qflll(self, long flag=0):
        """
        qflll(x,flag=0): LLL reduction of the vectors forming the matrix x
        (gives the unimodular transformation matrix). The columns of x must
        be linearly independent, unless specified otherwise below. flag is
        optional, and can be 0: default, 1: assumes x is integral, columns
        may be dependent, 2: assumes x is integral, returns a partially
        reduced basis, 4: assumes x is integral, returns [K,I] where K is
        the integer kernel of x and I the LLL reduced image, 5: same as 4
        but x may have polynomial coefficients, 8: same as 0 but x may have
        polynomial coefficients.
        """
        pari_catch_sig_on()
        return self.new_gen(qflll0(self.g,flag)).Mat()

    def qflllgram(self, long flag=0):
        """
        qflllgram(x,flag=0): LLL reduction of the lattice whose gram matrix
        is x (gives the unimodular transformation matrix). flag is optional
        and can be 0: default,1: lllgramint algorithm for integer matrices,
        4: lllgramkerim giving the kernel and the LLL reduced image, 5:
        lllgramkerimgen same when the matrix has polynomial coefficients,
        8: lllgramgen, same as qflllgram when the coefficients are
        polynomials.
        """
        pari_catch_sig_on()
        return self.new_gen(qflllgram0(self.g,flag)).Mat()

    def lllgram(self):
        return self.qflllgram(0)

    def lllgramint(self):
        return self.qflllgram(1)

    def qfminim(self, B, max, long flag=0):
        """
        qfminim(x,bound,maxnum,flag=0): number of vectors of square norm =
        bound, maximum norm and list of vectors for the integral and
        definite quadratic form x; minimal non-zero vectors if bound=0.
        flag is optional, and can be 0: default; 1: returns the first
        minimal vector found (ignore maxnum); 2: as 0 but uses a more
        robust, slower implementation, valid for non integral quadratic
        forms.
        """
        t0GEN(B)
        t1GEN(max)
        pari_catch_sig_on()
        return self.new_gen(qfminim0(self.g,t0,t1,flag,precdl))

    def qfrep(self, B, long flag=0):
        """
        qfrep(x,B,flag=0): vector of (half) the number of vectors of norms
        from 1 to B for the integral and definite quadratic form x. Binary
        digits of flag mean 1: count vectors of even norm from 1 to 2B, 2:
        return a t_VECSMALL instead of a t_VEC.
        """
        t0GEN(B)
        pari_catch_sig_on()
        return self.new_gen(qfrep0(self.g,t0,flag))

    def matsolve(self, B):
        """
        matsolve(B): Solve the linear system Mx=B for an invertible matrix
        M

        matsolve(B) uses Gaussian elimination to solve Mx=B, where M is
        invertible and B is a column vector.

        The corresponding pari library routine is gauss. The gp-interface
        name matsolve has been given preference here.

        INPUT:


        -  ``B`` - a column vector of the same dimension as the
           square matrix self


        EXAMPLES::

            sage: pari('[1,1;1,-1]').matsolve(pari('[1;0]'))
            [1/2; 1/2]
        """
        t0GEN(B)
        pari_catch_sig_on()
        return self.new_gen(gauss(self.g,t0))

    def matsolvemod(self, D, B, long flag = 0):
        r"""
        For column vectors `D=(d_i)` and `B=(b_i)`, find a small integer
        solution to the system of linear congruences

        .. math::

            R_ix=b_i\text{ (mod }d_i),

        where `R_i` is the ith row of ``self``. If `d_i=0`, the equation is
        considered over the integers. The entries of ``self``, ``D``, and
        ``B`` should all be integers (those of ``D`` should also be
        non-negative).

        If ``flag`` is 1, the output is a two-component row vector whose first
        component is a solution and whose second component is a matrix whose
        columns form a basis of the solution set of the homogeneous system.

        For either value of ``flag``, the output is 0 if there is no solution.

        Note that if ``D`` or ``B`` is an integer, then it will be considered
        as a vector all of whose entries are that integer.

        EXAMPLES::

            sage: D = pari('[3,4]~')
            sage: B = pari('[1,2]~')
            sage: M = pari('[1,2;3,4]')
            sage: M.matsolvemod(D, B)
            [-2, 0]~
            sage: M.matsolvemod(3, 1)
            [-1, 1]~
            sage: M.matsolvemod(pari('[3,0]~'), pari('[1,2]~'))
            [6, -4]~
            sage: M2 = pari('[1,10;9,18]')
            sage: M2.matsolvemod(3, pari('[2,3]~'), 1)
            [[0, -1]~, [-1, -2; 1, -1]]
            sage: M2.matsolvemod(9, pari('[2,3]~'))
            0
            sage: M2.matsolvemod(9, pari('[2,45]~'), 1)
            [[1, 1]~, [-1, -4; 1, -5]]
        """
        t0GEN(D)
        t1GEN(B)
        pari_catch_sig_on()
        return self.new_gen(matsolvemod0(self.g, t0, t1, flag))

    def matker(self, long flag=0):
        """
        Return a basis of the kernel of this matrix.

        INPUT:


        -  ``flag`` - optional; may be set to 0: default;
           non-zero: x is known to have integral entries.


        EXAMPLES::

            sage: pari('[1,2,3;4,5,6;7,8,9]').matker()
            [1; -2; 1]

        With algorithm 1, even if the matrix has integer entries the kernel
        need not be saturated (which is weird)::

            sage: pari('[1,2,3;4,5,6;7,8,9]').matker(1)
            [3; -6; 3]
            sage: pari('matrix(3,3,i,j,i)').matker()
            [-1, -1; 1, 0; 0, 1]
            sage: pari('[1,2,3;4,5,6;7,8,9]*Mod(1,2)').matker()
            [Mod(1, 2); Mod(0, 2); Mod(1, 2)]
        """
        pari_catch_sig_on()
        return self.new_gen(matker0(self.g, flag))

    def matkerint(self, long flag=0):
        """
        Return the integer kernel of a matrix.

        This is the LLL-reduced Z-basis of the kernel of the matrix x with
        integral entries.

        INPUT:


        -  ``flag`` - optional, and may be set to 0: default,
           uses a modified LLL, 1: uses matrixqz.


        EXAMPLES::

            sage: pari('[2,1;2,1]').matker()
            [-1/2; 1]
            sage: pari('[2,1;2,1]').matkerint()
            [1; -2]
            sage: pari('[2,1;2,1]').matkerint(1)
            [1; -2]
        """
        pari_catch_sig_on()
        return self.new_gen(matkerint0(self.g, flag))

    def matdet(self, long flag=0):
        """
        Return the determinant of this matrix.

        INPUT:


        -  ``flag`` - (optional) flag 0: using Gauss-Bareiss.
           1: use classical Gaussian elimination (slightly better for integer
           entries)


        EXAMPLES::

            sage: pari('[1,2; 3,4]').matdet(0)
            -2
            sage: pari('[1,2; 3,4]').matdet(1)
            -2
        """
        pari_catch_sig_on()
        return self.new_gen(det0(self.g, flag))

    def trace(self):
        """
        Return the trace of this PARI object.

        EXAMPLES::

            sage: pari('[1,2; 3,4]').trace()
            5
        """
        pari_catch_sig_on()
        return self.new_gen(gtrace(self.g))

    def mathnf(self, flag=0):
        """
        A.mathnf(flag=0): (upper triangular) Hermite normal form of A,
        basis for the lattice formed by the columns of A.

        INPUT:


        -  ``flag`` - optional, value range from 0 to 4 (0 if
           omitted), meaning : 0: naive algorithm

        -  ``1: Use Batut's algorithm`` - output 2-component
           vector [H,U] such that H is the HNF of A, and U is a unimodular
           matrix such that xU=H. 3: Use Batut's algorithm. Output [H,U,P]
           where P is a permutation matrix such that P A U = H. 4: As 1, using
           a heuristic variant of LLL reduction along the way.


        EXAMPLES::

            sage: pari('[1,2,3; 4,5,6;  7,8,9]').mathnf()
            [6, 1; 3, 1; 0, 1]
        """
        pari_catch_sig_on()
        return self.new_gen(mathnf0(self.g, flag))

    def mathnfmod(self, d):
        """
        Returns the Hermite normal form if d is a multiple of the
        determinant

        Beware that PARI's concept of a Hermite normal form is an upper
        triangular matrix with the same column space as the input matrix.

        INPUT:


        -  ``d`` - multiple of the determinant of self


        EXAMPLES::

                   sage: M=matrix([[1,2,3],[4,5,6],[7,8,11]])
            sage: d=M.det()
            sage: pari(M).mathnfmod(d)
                   [6, 4, 3; 0, 1, 0; 0, 0, 1]

        Note that d really needs to be a multiple of the discriminant, not
        just of the exponent of the cokernel::

                   sage: M=matrix([[1,0,0],[0,2,0],[0,0,6]])
            sage: pari(M).mathnfmod(6)
            [1, 0, 0; 0, 1, 0; 0, 0, 6]
            sage: pari(M).mathnfmod(12)
            [1, 0, 0; 0, 2, 0; 0, 0, 6]
        """
        t0GEN(d)
        pari_catch_sig_on()
        return self.new_gen(hnfmod(self.g, t0))

    def mathnfmodid(self, d):
        """
        Returns the Hermite Normal Form of M concatenated with d\*Identity

        Beware that PARI's concept of a Hermite normal form is a maximal
        rank upper triangular matrix with the same column space as the
        input matrix.

        INPUT:


        -  ``d`` - Determines


        EXAMPLES::

                   sage: M=matrix([[1,0,0],[0,2,0],[0,0,6]])
            sage: pari(M).mathnfmodid(6)
                   [1, 0, 0; 0, 2, 0; 0, 0, 6]

        This routine is not completely equivalent to mathnfmod::

            sage: pari(M).mathnfmod(6)
            [1, 0, 0; 0, 1, 0; 0, 0, 6]
        """
        t0GEN(d)
        pari_catch_sig_on()
        return self.new_gen(hnfmodid(self.g, t0))

    def matsnf(self, flag=0):
        """
        x.matsnf(flag=0): Smith normal form (i.e. elementary divisors) of
        the matrix x, expressed as a vector d. Binary digits of flag mean
        1: returns [u,v,d] where d=u\*x\*v, otherwise only the diagonal d
        is returned, 2: allow polynomial entries, otherwise assume x is
        integral, 4: removes all information corresponding to entries equal
        to 1 in d.

        EXAMPLES::

            sage: pari('[1,2,3; 4,5,6;  7,8,9]').matsnf()
            [0, 3, 1]
        """
        pari_catch_sig_on()
        return self.new_gen(matsnf0(self.g, flag))

    def matfrobenius(self, flag=0):
        r"""
        M.matfrobenius(flag=0): Return the Frobenius form of the square
        matrix M. If flag is 1, return only the elementary divisors (a list
        of polynomials). If flag is 2, return a two-components vector [F,B]
        where F is the Frobenius form and B is the basis change so that
        `M=B^{-1} F B`.

        EXAMPLES::

            sage: a = pari('[1,2;3,4]')
            sage: a.matfrobenius()
            [0, 2; 1, 5]
            sage: a.matfrobenius(flag=1)
            [x^2 - 5*x - 2]
            sage: a.matfrobenius(2)
            [[0, 2; 1, 5], [1, -1/3; 0, 1/3]]
            sage: v = a.matfrobenius(2)
            sage: v[0]
            [0, 2; 1, 5]
            sage: v[1]^(-1)*v[0]*v[1]
            [1, 2; 3, 4]

        We let t be the matrix of `T_2` acting on modular symbols
        of level 43, which was computed using
        ``ModularSymbols(43,sign=1).T(2).matrix()``::

            sage: t = pari('[3, -2, 0, 0; 0, -2, 0, 1; 0, -1, -2, 2; 0, -2, 0, 2]')
            sage: t.matfrobenius()
            [0, 0, 0, -12; 1, 0, 0, -2; 0, 1, 0, 8; 0, 0, 1, 1]
            sage: t.charpoly('x')
            x^4 - x^3 - 8*x^2 + 2*x + 12
            sage: t.matfrobenius(1)
            [x^4 - x^3 - 8*x^2 + 2*x + 12]

        AUTHORS:

        - Martin Albrect (2006-04-02)
        """
        pari_catch_sig_on()
        return self.new_gen(matfrobenius(self.g, flag, 0))


    ###########################################
    # polarit2.c
    ###########################################
    def factor(gen self, limit=-1, bint proof=1):
        """
        Return the factorization of x.

        INPUT:


        -  ``limit`` - (default: -1) is optional and can be set
           whenever x is of (possibly recursive) rational type. If limit is
           set return partial factorization, using primes up to limit (up to
           primelimit if limit=0).


        proof - (default: True) optional. If False (not the default),
        returned factors `<10^{15}` may only be pseudoprimes.

        .. note::

           In the standard PARI/GP interpreter and C-library the
           factor command *always* has proof=False, so beware!

        EXAMPLES::

            sage: pari('x^10-1').factor()
            [x - 1, 1; x + 1, 1; x^4 - x^3 + x^2 - x + 1, 1; x^4 + x^3 + x^2 + x + 1, 1]
            sage: pari(2^100-1).factor()
            [3, 1; 5, 3; 11, 1; 31, 1; 41, 1; 101, 1; 251, 1; 601, 1; 1801, 1; 4051, 1; 8101, 1; 268501, 1]
            sage: pari(2^100-1).factor(proof=False)
            [3, 1; 5, 3; 11, 1; 31, 1; 41, 1; 101, 1; 251, 1; 601, 1; 1801, 1; 4051, 1; 8101, 1; 268501, 1]

        We illustrate setting a limit::

            sage: pari(next_prime(10^50)*next_prime(10^60)*next_prime(10^4)).factor(10^5)
            [10007, 1; 100000000000000000000000000000000000000000000000151000000000700000000000000000000000000000000000000000000001057, 1]

        PARI doesn't have an algorithm for factoring multivariate
        polynomials::

            sage: pari('x^3 - y^3').factor()
            Traceback (most recent call last):
            ...
            PariError:  (7)
        """
        cdef int r
        if limit == -1 and typ(self.g) == t_INT and proof:
            pari_catch_sig_on()
            r = factorint_withproof_sage(&t0, self.g, ten_to_15)
            z = P.new_gen(t0)
            if not r:
                return z
            else:
                return _factor_int_when_pari_factor_failed(self, z)
        pari_catch_sig_on()
        return P.new_gen(factor0(self.g, limit))


    ###########################################
    # misc (classify when I know where they go)
    ###########################################

    def hilbert(x, y, p):
        cdef long ret
        t0GEN(y)
        t1GEN(p)
        pari_catch_sig_on()
        ret = hilbert0(x.g, t0, t1)
        pari_catch_sig_off()
        return ret

    def chinese(self, y):
        t0GEN(y)
        pari_catch_sig_on()
        return P.new_gen(chinese(self.g, t0))

    def order(self):
        pari_catch_sig_on()
        return P.new_gen(order(self.g))

    def znprimroot(self):
        """
        Return a primitive root modulo self, whenever it exists.

        This is a generator of the group `(\ZZ/n\ZZ)^*`, whenever
        this group is cyclic, i.e. if `n=4` or `n=p^k` or
        `n=2p^k`, where `p` is an odd prime and `k`
        is a natural number.

        INPUT:


        -  ``self`` - positive integer equal to 4, or a power
           of an odd prime, or twice a power of an odd prime


        OUTPUT: gen

        EXAMPLES::

            sage: pari(4).znprimroot()
            Mod(3, 4)
            sage: pari(10007^3).znprimroot()
            Mod(5, 1002101470343)
            sage: pari(2*109^10).znprimroot()
            Mod(236736367459211723407, 473472734918423446802)
        """
        pari_catch_sig_on()
        return P.new_gen(znprimroot0(self.g))

    def __abs__(self):
        return self.abs()

    def norm(gen self):
        pari_catch_sig_on()
        return P.new_gen(gnorm(self.g))

    def nextprime(gen self, bint add_one=0):
        """
        nextprime(x): smallest pseudoprime greater than or equal to `x`.
        If ``add_one`` is non-zero, return the smallest pseudoprime
        strictly greater than `x`.

        EXAMPLES::

            sage: pari(1).nextprime()
            2
            sage: pari(2).nextprime()
            2
            sage: pari(2).nextprime(add_one = 1)
            3
            sage: pari(2^100).nextprime()
            1267650600228229401496703205653
        """
        pari_catch_sig_on()
        if add_one:
            return P.new_gen(gnextprime(gaddsg(1,self.g)))
        return P.new_gen(gnextprime(self.g))

    def change_variable_name(self, var):
        """
        In ``self``, which must be a ``t_POL`` or ``t_SER``, set the
        variable to ``var``.  If the variable of ``self`` is already
        ``var``, then return ``self``.

        .. WARNING::

            You should be careful with variable priorities when
            applying this on a polynomial or series of which the
            coefficients have polynomial components.  To be safe, only
            use this function on polynomials with integer or rational
            coefficients.  For a safer alternative, use :meth:`subst`.

        EXAMPLES::

            sage: f = pari('x^3 + 17*x + 3')
            sage: f.change_variable_name("y")
            y^3 + 17*y + 3
            sage: f = pari('1 + 2*y + O(y^10)')
            sage: f.change_variable_name("q")
            1 + 2*q + O(q^10)
            sage: f.change_variable_name("y") is f
            True

        In PARI, ``I`` refers to the square root of -1, so it cannot be
        used as variable name.  Note the difference with :meth:`subst`::

            sage: f = pari('x^2 + 1')
            sage: f.change_variable_name("I")
            Traceback (most recent call last):
            ...
            PariError:  (5)
            sage: f.subst("x", "I")
            0
        """
        pari_catch_sig_on()
        cdef long n = P.get_var(var)
        pari_catch_sig_off()
        if varn(self.g) == n:
            return self
        if typ(self.g) != t_POL and typ(self.g) != t_SER:
            raise TypeError, "set_variable() only works for polynomials or power series"
        # Copy self and then change the variable in place
        cdef gen newg = P.new_gen_noclear(self.g)
        setvarn(newg.g, n)
        return newg

    def subst(self, var, z):
        """
        In ``self``, replace the variable ``var`` by the expression `z`.

        EXAMPLES::

            sage: x = pari("x"); y = pari("y")
            sage: f = pari('x^3 + 17*x + 3')
            sage: f.subst(x, y)
            y^3 + 17*y + 3
            sage: f.subst(x, "z")
            z^3 + 17*z + 3
            sage: f.subst(x, "z")^2
            z^6 + 34*z^4 + 6*z^3 + 289*z^2 + 102*z + 9
            sage: f.subst(x, "x+1")
            x^3 + 3*x^2 + 20*x + 21
            sage: f.subst(x, "xyz")
            xyz^3 + 17*xyz + 3
            sage: f.subst(x, "xyz")^2
            xyz^6 + 34*xyz^4 + 6*xyz^3 + 289*xyz^2 + 102*xyz + 9
        """
        cdef long n
        n = P.get_var(var)
        t0GEN(z)
        pari_catch_sig_on()
        return P.new_gen(gsubst(self.g, n, t0))

    def substpol(self, y, z):
        t0GEN(y)
        t1GEN(z)
        pari_catch_sig_on()
        return self.new_gen(gsubstpol(self.g, t0, t1))

    def nf_subst(self, z):
        """
        Given a PARI number field ``self``, return the same PARI
        number field but in the variable ``z``.

        INPUT:

        - ``self`` -- A PARI number field being the output of ``nfinit()``,
                      ``bnfinit()`` or ``bnrinit()``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K = NumberField(x^2 + 5, 'a')

        We can substitute in a PARI ``nf`` structure::

            sage: Kpari = K.pari_nf()
            sage: Kpari.nf_get_pol()
            y^2 + 5
            sage: Lpari = Kpari.nf_subst('a')
            sage: Lpari.nf_get_pol()
            a^2 + 5

        We can also substitute in a PARI ``bnf`` structure::

            sage: Kpari = K.pari_bnf()
            sage: Kpari.nf_get_pol()
            y^2 + 5
            sage: Kpari.bnf_get_cyc()  # Structure of class group
            [2]
            sage: Lpari = Kpari.nf_subst('a')
            sage: Lpari.nf_get_pol()
            a^2 + 5
            sage: Lpari.bnf_get_cyc()  # We still have a bnf after substituting
            [2]
        """
        cdef GEN nf = self.get_nf()
        t0GEN(z)
        pari_catch_sig_on()
        return P.new_gen(gsubst(self.g, nf_get_varn(nf), t0))

    def taylor(self, v=-1):
        pari_catch_sig_on()
        return self.new_gen(tayl(self.g, self.get_var(v), precdl))

    def thue(self, rhs, ne):
        t0GEN(rhs)
        t1GEN(ne)
        pari_catch_sig_on()
        return self.new_gen(thue(self.g, t0, t1))

    def charpoly(self, var=-1, flag=0):
        """
        charpoly(A,v=x,flag=0): det(v\*Id-A) = characteristic polynomial of
        A using the comatrix. flag is optional and may be set to 1 (use
        Lagrange interpolation) or 2 (use Hessenberg form), 0 being the
        default.
        """
        pari_catch_sig_on()
        return P.new_gen(charpoly0(self.g, P.get_var(var), flag))


    def kronecker(gen self, y):
        t0GEN(y)
        pari_catch_sig_on()
        return P.new_gen(gkronecker(self.g, t0))


    def type(gen self):
        """
        Return the Pari type of self as a string.

        .. note::

           In Cython, it is much faster to simply use typ(self.g) for
           checking Pari types.

        EXAMPLES::

            sage: pari(7).type()
            't_INT'
            sage: pari('x').type()
            't_POL'
        """
        # The following original code leaks memory:
        #        return str(type_name(typ(self.g)))
        #
        # This code is the usual workaround:
        #        cdef char* s= <char*>type_name(typ(self.g))
        #        t=str(s)
        #        free(s)
        #        return(t)
        # However, it causes segfaults with t_INTs on some
        # machines, and errors about freeing non-aligned
        # pointers on others. So we settle for the following
        # fast but ugly code. Note that should the list of
        # valid Pari types ever be updated, this code would
        # need to be updated accordingly.
        #
        cdef long t = typ(self.g)

        if   t == t_INT:      return 't_INT'
        elif t == t_REAL:     return 't_REAL'
        elif t == t_INTMOD:   return 't_INTMOD'
        elif t == t_FRAC:     return 't_FRAC'
        elif t == t_FFELT:    return 't_FFELT'
        elif t == t_COMPLEX:  return 't_COMPLEX'
        elif t == t_PADIC:    return 't_PADIC'
        elif t == t_QUAD:     return 't_QUAD'
        elif t == t_POLMOD:   return 't_POLMOD'
        elif t == t_POL:      return 't_POL'
        elif t == t_SER:      return 't_SER'
        elif t == t_RFRAC:    return 't_RFRAC'
        elif t == t_QFR:      return 't_QFR'
        elif t == t_QFI:      return 't_QFI'
        elif t == t_VEC:      return 't_VEC'
        elif t == t_COL:      return 't_COL'
        elif t == t_MAT:      return 't_MAT'
        elif t == t_LIST:     return 't_LIST'
        elif t == t_STR:      return 't_STR'
        elif t == t_VECSMALL: return 't_VECSMALL'
        elif t == t_CLOSURE:  return 't_CLOSURE'
        else:
            raise TypeError, "Unknown Pari type: %s"%t


    def polinterpolate(self, ya, x):
        """
        self.polinterpolate(ya,x,e): polynomial interpolation at x
        according to data vectors self, ya (i.e. return P such that
        P(self[i]) = ya[i] for all i). Also return an error estimate on the
        returned value.
        """
        t0GEN(ya)
        t1GEN(x)
        cdef GEN dy, g
        pari_catch_sig_on()
        g = polint(self.g, t0, t1, &dy)
        dif = self.new_gen_noclear(dy)
        return self.new_gen(g), dif

    def algdep(self, long n):
        """
        EXAMPLES::

            sage: n = pari.set_real_precision(210)
            sage: w1 = pari('z1=2-sqrt(26); (z1+I)/(z1-I)')
            sage: f = w1.algdep(12); f
            545*x^11 - 297*x^10 - 281*x^9 + 48*x^8 - 168*x^7 + 690*x^6 - 168*x^5 + 48*x^4 - 281*x^3 - 297*x^2 + 545*x
            sage: f(w1).abs() < 1.0e-200
            True
            sage: f.factor()
            [x, 1; x + 1, 2; x^2 + 1, 1; x^2 + x + 1, 1; 545*x^4 - 1932*x^3 + 2790*x^2 - 1932*x + 545, 1]
            sage: pari.set_real_precision(n)
            210
        """
        pari_catch_sig_on()
        return self.new_gen(algdep(self.g, n))

    def concat(self, y):
        t0GEN(y)
        pari_catch_sig_on()
        return self.new_gen(concat(self.g, t0))

    def lindep(self, flag=0):
        pari_catch_sig_on()
        return self.new_gen(lindep0(self.g, flag))

    def listinsert(self, obj, long n):
        t0GEN(obj)
        pari_catch_sig_on()
        return self.new_gen(listinsert(self.g, t0, n))

    def listput(self, obj, long n):
        t0GEN(obj)
        pari_catch_sig_on()
        return self.new_gen(listput(self.g, t0, n))



    def elleisnum(self, long k, int flag=0):
        """
        om.elleisnum(k, flag=0): om=[om1,om2] being a 2-component vector
        giving a basis of a lattice L and k an even positive integer,
        computes the numerical value of the Eisenstein series of weight k.
        When flag is non-zero and k=4 or 6, this gives g2 or g3 with the
        correct normalization.

        INPUT:


        -  ``om`` - gen, 2-component vector giving a basis of a
           lattice L

        -  ``k`` - int (even positive)

        -  ``flag`` - int (default 0)


        OUTPUT:


        -  ``gen`` - numerical value of E_k


        EXAMPLES::

            sage: e = pari([0,1,1,-2,0]).ellinit()
            sage: om = e.omega()
            sage: om
            [2.49021256085506, -1.97173770155165*I]
            sage: om.elleisnum(2) # was:  -5.28864933965426
            10.0672605281120
            sage: om.elleisnum(4)
            112.000000000000
            sage: om.elleisnum(100)
            2.15314248576078 E50
        """
        pari_catch_sig_on()
        # the argument prec has no effect
        return self.new_gen(elleisnum(self.g, k, flag, prec))

    def ellwp(gen self, z='z', long n=20, long flag=0):
        """
        Return the value or the series expansion of the Weierstrass
        `P`-function at `z` on the lattice `self` (or the lattice
        defined by the elliptic curve `self`).

        INPUT:

        -  ``self`` -- an elliptic curve created using ``ellinit`` or a
           list ``[om1, om2]`` representing generators for a lattice.

        -  ``z`` -- (default: 'z') a complex number or a variable name
           (as string or PARI variable).

        -  ``n`` -- (default: 20) if 'z' is a variable, compute the
           series expansion up to at least `O(z^n)`.

        -  ``flag`` -- (default = 0): If ``flag`` is 0, compute only
           `P(z)`.  If ``flag`` is 1, compute `[P(z), P'(z)]`.

        OUTPUT:

        - `P(z)` (if ``flag`` is 0) or `[P(z), P'(z)]` (if ``flag`` is 1).
           numbers

        EXAMPLES:

        We first define the elliptic curve X_0(11)::

            sage: E = pari([0,-1,1,-10,-20]).ellinit()

        Compute P(1)::

            sage: E.ellwp(1)
            13.9658695257485 + 0.E-18*I

        Compute P(1+i), where i = sqrt(-1)::

            sage: C.<i> = ComplexField()
            sage: E.ellwp(pari(1+i))
            -1.11510682565555 + 2.33419052307470*I
            sage: E.ellwp(1+i)
            -1.11510682565555 + 2.33419052307470*I

        The series expansion, to the default `O(z^20)` precision::

            sage: E.ellwp()
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8 + 1202285717/928746000*z^10 + 2403461/2806650*z^12 + 30211462703/43418875500*z^14 + 3539374016033/7723451736000*z^16 + 413306031683977/1289540602350000*z^18 + O(z^20)

        Compute the series for wp to lower precision::

            sage: E.ellwp(n=4)
            z^-2 + 31/15*z^2 + O(z^4)

        Next we use the version where the input is generators for a
        lattice::

            sage: pari([1.2692, 0.63 + 1.45*i]).ellwp(1)
            13.9656146936689 + 0.000644829272810...*I

        With flag=1, compute the pair P(z) and P'(z)::

            sage: E.ellwp(1, flag=1)
            [13.9658695257485 + 0.E-18*I, 50.5619300880073 ... E-18*I]
        """
        t0GEN(z)
        pari_catch_sig_on()
        cdef long dprec
        dprec = gprecision(t0)
        if dprec:
            dprec = prec_words_to_dec(dprec)
        else:
            dprec = prec
        return self.new_gen(ellwp0(self.g, t0, flag, n+2, dprec))

    def ellchangepoint(self, y):
        """
        self.ellchangepoint(y): change data on point or vector of points
        self on an elliptic curve according to y=[u,r,s,t]

        EXAMPLES::

            sage: e = pari([0,1,1,-2,0]).ellinit()
            sage: x = pari([1,0])
            sage: e.ellisoncurve([1,4])
            False
            sage: e.ellisoncurve(x)
            True
            sage: f = e.ellchangecurve([1,2,3,-1])
            sage: f[:5]   # show only first five entries
            [6, -2, -1, 17, 8]
            sage: x.ellchangepoint([1,2,3,-1])
            [-1, 4]
            sage: f.ellisoncurve([-1,4])
            True
        """
        t0GEN(y)
        pari_catch_sig_on()
        return self.new_gen(ellchangepoint(self.g, t0))

    def debug(gen self, long depth = -1):
        r"""
        Show the internal structure of self (like the ``\x`` command in gp).

        EXAMPLE::

            sage: pari('[1/2, 1.0*I]').debug()  # random addresses
            [&=0000000004c5f010] VEC(lg=3):2200000000000003 0000000004c5eff8 0000000004c5efb0
              1st component = [&=0000000004c5eff8] FRAC(lg=3):0800000000000003 0000000004c5efe0 0000000004c5efc8
                num = [&=0000000004c5efe0] INT(lg=3):0200000000000003 (+,lgefint=3):4000000000000003 0000000000000001
                den = [&=0000000004c5efc8] INT(lg=3):0200000000000003 (+,lgefint=3):4000000000000003 0000000000000002
              2nd component = [&=0000000004c5efb0] COMPLEX(lg=3):0c00000000000003 00007fae8a2eb840 0000000004c5ef90
                real = gen_0
                imag = [&=0000000004c5ef90] REAL(lg=4):0400000000000004 (+,expo=0):6000000000000000 8000000000000000 0000000000000000
        """
        pari_catch_sig_on()
        dbgGEN(self.g, depth)
        pari_catch_sig_off()
        return


    ##################################################
    # Technical functions that can be used by other
    # classes that derive from gen.
    ##################################################
    cdef gen pari(self, object x):
        return pari(x)

    cdef gen new_gen(self, GEN x):
        return P.new_gen(x)

    cdef gen new_gen_noclear(self, GEN x):
        return P.new_gen_noclear(x)

    cdef GEN _deepcopy_to_python_heap(self, GEN x, pari_sp* address):
        return P.deepcopy_to_python_heap(x, address)

    cdef long get_var(self, v):
        return P.get_var(v)



cdef unsigned long num_primes

# Callbacks from PARI to print stuff using sys.stdout.write() instead
# of C library functions like puts().
cdef PariOUT sage_pariOut

cdef void sage_putchar(char c):
    cdef char str[2]
    str[0] = c
    str[1] = 0
    sys.stdout.write(str)
    return

cdef void sage_puts(char* s):
    sys.stdout.write(s)
    return

cdef void sage_flush():
    sys.stdout.flush()
    return


cdef class PariInstance(sage.structure.parent_base.ParentWithBase):
    def __init__(self, long size=16000000, unsigned long maxprime=500000):
        """
        Initialize the PARI system.

        INPUT:


        -  ``size`` - long, the number of bytes for the initial
           PARI stack (see note below)

        -  ``maxprime`` - unsigned long, upper limit on a
           precomputed prime number table (default: 500000)


        .. note::

           In py_pari, the PARI stack is different than in gp or the
           PARI C library. In Python, instead of the PARI stack
           holding the results of all computations, it *only* holds
           the results of an individual computation. Each time a new
           Python/PARI object is computed, it it copied to its own
           space in the Python heap, and the memory it occupied on the
           PARI stack is freed. Thus it is not necessary to make the
           stack very large. Also, unlike in PARI, if the stack does
           overflow, in most cases the PARI stack is automatically
           increased and the relevant step of the computation rerun.

           This design obviously involves some performance penalties
           over the way PARI works, but it scales much better and is
           far more robust for large projects.

        .. note::

           If you do not want prime numbers, put ``maxprime=2``, but be
           careful because many PARI functions require this table. If
           you get the error message "not enough precomputed primes",
           increase this parameter.
        """
        if bot:
            return  # pari already initialized.

        global num_primes, avma, top, bot, prec

        # The size here doesn't really matter, because we will allocate
        # our own stack anyway. We ask PARI not to set up signal handlers.
        pari_init_opts(10000, maxprime, INIT_JMPm | INIT_DFTm)
        num_primes = maxprime

        # NOTE: pari_catch_sig_on() can only come AFTER pari_init_opts()!
        pari_catch_sig_on()

        # Free the PARI stack and allocate our own (using Cython)
        pari_free(<void*>bot); bot = 0
        init_stack(size)

        GP_DATA.fmt.prettyp = 0

        # how do I get the following to work? seems to be a circular import
        #from sage.rings.real_mpfr import RealField
        #prec_bits = RealField().prec()
        prec = prec_bits_to_words(53)
        GP_DATA.fmt.sigd = prec_bits_to_dec(53)

        # Set printing functions
        global pariOut
        pariOut = &sage_pariOut
        pariOut.putch = sage_putchar
        pariOut.puts = sage_puts
        pariOut.flush = sage_flush
        pari_catch_sig_off()

    def __dealloc__(self):
        """
        Deallocation of the Pari instance.

        NOTE:

        Usually this deallocation happens only when Sage quits.
        We do not provide a direct test, since usually there
        is only one Pari instance, and when artificially creating
        another instance, C-data are shared.

        The fact that Sage does not crash when quitting is an
        indirect doctest. See the discussion at :trac:`13741`.

        """
        if bot:
            sage_free(<void*>bot)
        global top, bot
        top = 0
        bot = 0
        pari_close()

    def __repr__(self):
        return "Interface to the PARI C library"

    def __hash__(self):
        return 907629390   # hash('pari')

    cdef has_coerce_map_from_c_impl(self, x):
        return True

    def __richcmp__(left, right, int op):
        """
        EXAMPLES::

            sage: pari == pari
            True
            sage: pari == gp
            False
            sage: pari == 5
            False
        """
        return (<Parent>left)._richcmp(right, op)

    def default(self, variable, value=None):
        if not value is None:
            return self('default(%s, %s)'%(variable, value))
        return self('default(%s)'%variable)

    def set_debug_level(self, level):
        """
        Set the debug PARI C library variable.
        """
        self.default('debug', int(level))

    def get_debug_level(self):
        """
        Set the debug PARI C library variable.
        """
        return int(self.default('debug'))

    cdef GEN toGEN(self, x, int i) except NULL:
        cdef gen _x
        if PY_TYPE_CHECK(x, gen):
            _x = x
            return _x.g

        t0heap[i] = self(x)
        _x = t0heap[i]
        return _x.g

    def set_real_precision(self, long n):
        """
        Sets the PARI default real precision.

        This is used both for creation of new objects from strings and for
        printing. It is the number of digits *IN DECIMAL* in which real
        numbers are printed. It also determines the precision of objects
        created by parsing strings (e.g. pari('1.2')), which is *not* the
        normal way of creating new pari objects in Sage. It has *no*
        effect on the precision of computations within the pari library.

        Returns the previous PARI real precision.
        """
        cdef unsigned long k

        k = GP_DATA.fmt.sigd
        s = str(n)
        pari_catch_sig_on()
        sd_realprecision(s, 2)
        pari_catch_sig_off()
        return int(k)  # Python int

    def get_real_precision(self):
        """
        Returns the current PARI default real precision.

        This is used both for creation of new objects from strings and for
        printing. It is the number of digits *IN DECIMAL* in which real
        numbers are printed. It also determines the precision of objects
        created by parsing strings (e.g. pari('1.2')), which is *not* the
        normal way of creating new pari objects in Sage. It has *no*
        effect on the precision of computations within the pari library.
        """
        return GP_DATA.fmt.sigd

    def set_series_precision(self, long n):
        global precdl
        precdl = n

    def get_series_precision(self):
        return precdl


    ###########################################
    # Create a gen from a GEN object.
    # This *steals* a reference to the GEN, as it
    # frees the memory the GEN occupied.
    ###########################################

    cdef gen new_gen(self, GEN x):
        """
        Create a new gen, then free the \*entire\* stack and call
        pari_catch_sig_off().
        """
        cdef gen g
        g = _new_gen(x)
        global mytop, avma
        avma = mytop
        pari_catch_sig_off()
        return g

    cdef object new_gen_to_string(self, GEN x):
        """
        Converts a gen to a Python string, free the \*entire\* stack and call
        pari_catch_sig_off(). This is meant to be used in place of new_gen().
        """
        cdef char* c
        cdef int n
        c = GENtostr(x)
        s = str(c)
        pari_free(c)
        global mytop, avma
        avma = mytop
        pari_catch_sig_off()
        return s

    cdef void clear_stack(self):
        """
        Clear the entire PARI stack and call pari_catch_sig_off().
        """
        global mytop, avma
        avma = mytop
        pari_catch_sig_off()

    cdef void set_mytop_to_avma(self):
        global mytop, avma
        mytop = avma

    cdef gen new_gen_noclear(self, GEN x):
        """
        Create a new gen, but don't free any memory on the stack and don't
        call pari_catch_sig_off().
        """
        z = _new_gen(x)
        return z

    cdef gen new_gen_from_mpz_t(self, mpz_t value):
        """
        Create a new gen from a given MPIR-integer ``value``.

        EXAMPLES::

            sage: pari(42)       # indirect doctest
            42

        TESTS:

        Check that the hash of an integer does not depend on existing
        garbage on the stack (#11611)::

            sage: foo = pari(2^(32*1024));  # Create large integer to put PARI stack in known state
            sage: a5 = pari(5);
            sage: foo = pari(0xDEADBEEF * (2^(32*1024)-1)//(2^32 - 1));  # Dirty PARI stack
            sage: b5 = pari(5);
            sage: a5.__hash__() == b5.__hash__()
            True
        """
        pari_catch_sig_on()
        return self.new_gen(self._new_GEN_from_mpz_t(value))

    cdef inline GEN _new_GEN_from_mpz_t(self, mpz_t value):
        r"""
        Create a new PARI ``t_INT`` from a ``mpz_t``.

        For internal use only; this directly uses the PARI stack.
        One should call ``pari_catch_sig_on()`` before and ``pari_catch_sig_off()`` after.
        """
        cdef unsigned long limbs = mpz_size(value)

        cdef GEN z = cgeti(limbs + 2)
        # Set sign and "effective length"
        z[1] = evalsigne(mpz_sgn(value)) + evallgefint(limbs + 2)
        mpz_export(int_LSW(z), NULL, -1, sizeof(long), 0, 0, value)

        return z

    cdef gen new_gen_from_int(self, int value):
        pari_catch_sig_on()
        return self.new_gen(stoi(value))

    cdef gen new_gen_from_mpq_t(self, mpq_t value):
        """
        Create a new gen from a given MPIR-rational ``value``.

        EXAMPLES::

            sage: pari(-2/3)
            -2/3
            sage: pari(QQ(42))
            42
            sage: pari(QQ(42)).type()
            't_INT'
            sage: pari(QQ(1/42)).type()
            't_FRAC'

        TESTS:

        Check that the hash of a rational does not depend on existing
        garbage on the stack (#11854)::

            sage: foo = pari(2^(32*1024));  # Create large integer to put PARI stack in known state
            sage: a5 = pari(5/7);
            sage: foo = pari(0xDEADBEEF * (2^(32*1024)-1)//(2^32 - 1));  # Dirty PARI stack
            sage: b5 = pari(5/7);
            sage: a5.__hash__() == b5.__hash__()
            True
        """
        pari_catch_sig_on()
        return self.new_gen(self._new_GEN_from_mpq_t(value))

    cdef inline GEN _new_GEN_from_mpq_t(self, mpq_t value):
        r"""
        Create a new PARI ``t_INT`` or ``t_FRAC`` from a ``mpq_t``.

        For internal use only; this directly uses the PARI stack.
        One should call ``pari_catch_sig_on()`` before and ``pari_catch_sig_off()`` after.
        """
        cdef GEN num = self._new_GEN_from_mpz_t(mpq_numref(value))
        if mpz_cmpabs_ui(mpq_denref(value), 1) == 0:
            # Denominator is 1, return the numerator (an integer)
            return num
        cdef GEN denom = self._new_GEN_from_mpz_t(mpq_denref(value))
        return mkfrac(num, denom)

    cdef gen new_t_POL_from_int_star(self, int *vals, int length, long varnum):
        """
        Note that degree + 1 = length, so that recognizing 0 is easier.

        varnum = 0 is the general choice (creates a variable in x).
        """
        cdef GEN z
        cdef int i

        pari_catch_sig_on()
        z = cgetg(length + 2, t_POL)
        z[1] = evalvarn(varnum)
        if length != 0:
            setsigne(z,1)
            for i from 0 <= i < length:
                set_gel(z,i+2, stoi(vals[i]))
        else:
            ## polynomial is zero
            setsigne(z,0)

        return self.new_gen(z)

    cdef gen new_gen_from_padic(self, long ordp, long relprec,
                                mpz_t prime, mpz_t p_pow, mpz_t unit):
        cdef GEN z
        pari_catch_sig_on()
        z = cgetg(5, t_PADIC)
        z[1] = evalprecp(relprec) + evalvalp(ordp)
        set_gel(z, 2, self._new_GEN_from_mpz_t(prime))
        set_gel(z, 3, self._new_GEN_from_mpz_t(p_pow))
        set_gel(z, 4, self._new_GEN_from_mpz_t(unit))
        return self.new_gen(z)

    def double_to_gen(self, x):
        cdef double dx
        dx = float(x)
        return self.double_to_gen_c(dx)

    cdef gen double_to_gen_c(self, double x):
        """
        Create a new gen with the value of the double x, using Pari's
        dbltor.

        EXAMPLES::

            sage: pari.double_to_gen(1)
            1.00000000000000
            sage: pari.double_to_gen(1e30)
            1.00000000000000 E30
            sage: pari.double_to_gen(0)
            0.E-15
            sage: pari.double_to_gen(-sqrt(RDF(2)))
            -1.41421356237310
        """
        # Pari has an odd concept where it attempts to track the accuracy
        # of floating-point 0; a floating-point zero might be 0.0e-20
        # (meaning roughly that it might represent any number in the
        # range -1e-20 <= x <= 1e20).

        # Pari's dbltor converts a floating-point 0 into the Pari real
        # 0.0e-307; Pari treats this as an extremely precise 0.  This
        # can cause problems; for instance, the Pari incgam() function can
        # be very slow if the first argument is very precise.

        # So we translate 0 into a floating-point 0 with 53 bits
        # of precision (that's the number of mantissa bits in an IEEE
        # double).

        pari_catch_sig_on()
        if x == 0:
            return self.new_gen(real_0_bit(-53))
        else:
            return self.new_gen(dbltor(x))

    cdef GEN double_to_GEN(self, double x):
        if x == 0:
            return real_0_bit(-53)
        else:
            return dbltor(x)

    def complex(self, re, im):
        """
        Create a new complex number, initialized from re and im.
        """
        t0GEN(re)
        t1GEN(im)
        cdef GEN cp
        pari_catch_sig_on()
        cp = cgetg(3, t_COMPLEX)
        set_gel(cp, 1, t0)
        set_gel(cp, 2, t1)
        return self.new_gen(cp)

    cdef GEN deepcopy_to_python_heap(self, GEN x, pari_sp* address):
        return deepcopy_to_python_heap(x, address)

    cdef gen new_ref(self, GEN g, gen parent):
        """
        Create a new gen pointing to the given GEN, which is allocated as a
        part of parent.g.

        .. note::

           As a rule, there should never be more than one sage gen
           pointing to a given Pari GEN. So that means there is only
           one case where this function should be used: when a
           complicated Pari GEN is allocated with a single gen
           pointing to it, and one needs a gen pointing to one of its
           components.

           For example, doing x = pari("[1,2]") allocates a gen pointing to
           the list [1,2], but x[0] has no gen wrapping it, so new_ref
           should be used there. Then parent would be x in this
           case. See __getitem__ for an example of usage.

        EXAMPLES::

            sage: pari("[[1,2],3]")[0][1] ## indirect doctest
            2
        """
        cdef gen p = PY_NEW(gen)

        p.b = 0
        p._parent = self
        p._refers_to = {-1:parent}
        p.g = g

        return p

    def __call__(self, s):
        """
        Create the PARI object obtained by evaluating s using PARI.

        EXAMPLES::

            sage: pari([2,3,5])
            [2, 3, 5]
            sage: pari(Matrix(2,2,range(4)))
            [0, 1; 2, 3]
            sage: pari(x^2-3)
            x^2 - 3

        ::

            sage: a = pari(1); a, a.type()
            (1, 't_INT')
            sage: a = pari(1/2); a, a.type()
            (1/2, 't_FRAC')
            sage: a = pari(1/2); a, a.type()
            (1/2, 't_FRAC')

        See :func:`pari` for more examples.
        """
        cdef int length, i
        cdef gen v

        late_import()

        if isinstance(s, gen):
            return s
        elif isinstance(s, Integer):
            return self.new_gen_from_mpz_t(<void *>s + mpz_t_offset)
        elif PyObject_HasAttrString(s, "_pari_"):
            return s._pari_()

        # Check basic Python types
        if PyInt_Check(s):
            pari_catch_sig_on()
            return self.new_gen(stoi(PyInt_AS_LONG(s)))
        if PyBool_Check(s):
            return self.PARI_ONE if s else self.PARI_ZERO
        cdef mpz_t mpz_int
        cdef GEN g
        if PyLong_Check(s):
            pari_catch_sig_on()
            mpz_init(mpz_int)
            mpz_set_pylong(mpz_int, s)
            g = self._new_GEN_from_mpz_t(mpz_int)
            mpz_clear(mpz_int)
            return self.new_gen(g)
        if PyFloat_Check(s):
            pari_catch_sig_on()
            return self.new_gen(dbltor(PyFloat_AS_DOUBLE(s)))
        if PyComplex_Check(s):
            pari_catch_sig_on()
            g = cgetg(3, t_COMPLEX)
            set_gel(g, 1, dbltor(PyComplex_RealAsDouble(s)))
            set_gel(g, 2, dbltor(PyComplex_ImagAsDouble(s)))
            return self.new_gen(g)

        if isinstance(s, (types.ListType, types.XRangeType,
                            types.TupleType, types.GeneratorType)):
            length = len(s)
            v = self._empty_vector(length)
            for i from 0 <= i < length:
                v[i] = self(s[i])
            return v

        t = str(s)
        pari_catch_sig_str('evaluating PARI string')
        g = gp_read_str(t)
        if g == gnil:
            pari_catch_sig_off()
            return None
        return self.new_gen(g)

    cdef GEN _new_GEN_from_mpz_t_matrix(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc):
        r"""
        Create a new PARI ``t_MAT`` with ``nr`` rows and ``nc`` columns
        from a ``mpz_t**``.

        For internal use only; this directly uses the PARI stack.
        One should call ``pari_catch_sig_on()`` before and ``pari_catch_sig_off()`` after.
        """
        cdef GEN x
        cdef GEN A = zeromatcopy(nr, nc)
        cdef Py_ssize_t i, j
        for i in range(nr):
            for j in range(nc):
                x = self._new_GEN_from_mpz_t(B[i][j])
                set_gcoeff(A, i+1, j+1, x)  # A[i+1, j+1] = x (using 1-based indexing)
        return A

    cdef GEN _new_GEN_from_mpz_t_matrix_rotate90(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc):
        r"""
        Create a new PARI ``t_MAT`` with ``nr`` rows and ``nc`` columns
        from a ``mpz_t**`` and rotate the matrix 90 degrees
        counterclockwise.  So the resulting matrix will have ``nc`` rows
        and ``nr`` columns.  This is useful for computing the Hermite
        Normal Form because Sage and PARI use different definitions.

        For internal use only; this directly uses the PARI stack.
        One should call ``pari_catch_sig_on()`` before and ``pari_catch_sig_off()`` after.
        """
        cdef GEN x
        cdef GEN A = zeromatcopy(nc, nr)
        cdef Py_ssize_t i, j
        for i in range(nr):
            for j in range(nc):
                x = self._new_GEN_from_mpz_t(B[i][nc-j-1])
                set_gcoeff(A, j+1, i+1, x)  # A[j+1, i+1] = x (using 1-based indexing)
        return A

    cdef gen integer_matrix(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc, bint permute_for_hnf):
        """
        EXAMPLES::

            sage: matrix(ZZ,2,[1..6])._pari_()   # indirect doctest
            [1, 2, 3; 4, 5, 6]
        """
        pari_catch_sig_on()
        cdef GEN g
        if permute_for_hnf:
            g = self._new_GEN_from_mpz_t_matrix_rotate90(B, nr, nc)
        else:
            g = self._new_GEN_from_mpz_t_matrix(B, nr, nc)
        return self.new_gen(g)

    cdef GEN _new_GEN_from_mpq_t_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc):
        cdef GEN x
        # Allocate zero matrix
        cdef GEN A = zeromatcopy(nr, nc)
        cdef Py_ssize_t i, j
        for i in range(nr):
            for j in range(nc):
                x = self._new_GEN_from_mpq_t(B[i][j])
                set_gcoeff(A, i+1, j+1, x)  # A[i+1, j+1] = x (using 1-based indexing)
        return A

    cdef gen rational_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc):
        """
        EXAMPLES::

            sage: matrix(QQ,2,[1..6])._pari_()   # indirect doctest
            [1, 2, 3; 4, 5, 6]
        """
        pari_catch_sig_on()
        cdef GEN g = self._new_GEN_from_mpq_t_matrix(B, nr, nc)
        return self.new_gen(g)

    cdef _coerce_c_impl(self, x):
        """
        Implicit canonical coercion into a PARI object.
        """
        try:
            return self(x)
        except (TypeError, AttributeError):
            raise TypeError("no canonical coercion of %s into PARI"%x)

    cdef _an_element_c_impl(self):  # override this in Cython
        return self.PARI_ZERO

    def new_with_bits_prec(self, s, long precision):
        r"""
        pari.new_with_bits_prec(self, s, precision) creates s as a PARI
        gen with (at most) precision *bits* of precision.
        """
        cdef unsigned long old_prec
        old_prec = GP_DATA.fmt.sigd
        precision = prec_bits_to_dec(precision)
        if not precision:
            precision = old_prec
        self.set_real_precision(precision)
        x = self(s)
        self.set_real_precision(old_prec)
        return x



    cdef long get_var(self, v):
        """
        Converts a Python string into a PARI variable reference number. Or
        if v = -1, returns -1.
        """
        if v != -1:
            s = str(v)
            return fetch_user_var(s)
        return -1

    ############################################################
    # Initialization
    ############################################################

    def allocatemem(self, s=0, silent=False):
        r"""
        Double the *PARI* stack.
        """
        if s == 0 and not silent:
            print "Doubling the PARI stack."
        s = int(s)
        cdef size_t a = s
        if int(a) != s:
            raise ValueError, "s must be nonnegative and not too big."
        init_stack(s)

    def pari_version(self):
        return str(PARIVERSION)

    def init_primes(self, _M):
        """
        Recompute the primes table including at least all primes up to M
        (but possibly more).

        EXAMPLES::

            sage: pari.init_primes(200000)

        We make sure that ticket #11741 has been fixed, and double check to
        make sure that diffptr has not been freed::

            sage: pari.init_primes(2^62)
            Traceback (most recent call last):
            ...
            PariError: not enough memory (28)            # 64-bit
            OverflowError: long int too large to convert # 32-bit
            sage: pari.init_primes(200000)
        """
        cdef unsigned long M
        cdef char *tmpptr
        M = _M
        global diffptr, num_primes
        if M <= num_primes:
            return
        pari_catch_sig_on()
        tmpptr = initprimes(M)
        pari_catch_sig_off()
        pari_free(<void*> diffptr)
        num_primes = M
        diffptr = tmpptr


    ##############################################
    ## Support for GP Scripts
    ##############################################

    def read(self, bytes filename):
        r"""
        Read a script from the named filename into the interpreter.  The
        functions defined in the script are then available for use from
        Sage/PARI.  The result of the last command in ``filename`` is
        returned.

        EXAMPLES:

        Create a gp file::

            sage: import tempfile
            sage: gpfile = tempfile.NamedTemporaryFile(mode="w")
            sage: gpfile.file.write("mysquare(n) = {\n")
            sage: gpfile.file.write("    n^2;\n")
            sage: gpfile.file.write("}\n")
            sage: gpfile.file.write("polcyclo(5)\n")
            sage: gpfile.file.flush()

        Read it in Sage, we get the result of the last line::

            sage: pari.read(gpfile.name)
            x^4 + x^3 + x^2 + x + 1

        Call the function defined in the gp file::

            sage: pari('mysquare(12)')
            144
        """
        pari_catch_sig_on()
        return self.new_gen(gp_read_file(filename))


    ##############################################

    def _primelimit(self):
        """
        Return the number of primes already computed
        in this Pari instance.

        EXAMPLES:
            sage: pari._primelimit()
            500519
            sage: pari.init_primes(600000)
            sage: pari._primelimit()
            600000
        """
        global num_primes
        from sage.rings.all import ZZ
        return ZZ(num_primes)

    def prime_list(self, long n):
        """
        prime_list(n): returns list of the first n primes

        To extend the table of primes use pari.init_primes(M).

        INPUT:


        -  ``n`` - C long


        OUTPUT:


        -  ``gen`` - PARI list of first n primes


        EXAMPLES::

            sage: pari.prime_list(0)
            []
            sage: pari.prime_list(-1)
            []
            sage: pari.prime_list(3)
            [2, 3, 5]
            sage: pari.prime_list(10)
            [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
            sage: pari.prime_list(20)
            [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
            sage: len(pari.prime_list(1000))
            1000
        """
        if n >= 2:
            self.nth_prime(n)
        pari_catch_sig_on()
        return self.new_gen(primes(n))

    def primes_up_to_n(self, long n):
        """
        Return the primes <= n as a pari list.

        EXAMPLES::

            sage: pari.primes_up_to_n(1)
            []
            sage: pari.primes_up_to_n(20)
            [2, 3, 5, 7, 11, 13, 17, 19]
        """
        if n <= 1:
            return pari([])
        self.init_primes(n+1)
        return self.prime_list(pari(n).primepi())

##         cdef long k
##         k = (n+10)/math.log(n)
##         p = 2
##         while p <= n:
##             p = self.nth_prime(k)
##             k = 2
##         v = self.prime_list(k)
##         return v[:pari(n).primepi()]

    def __nth_prime(self, long n):
        """
        nth_prime(n): returns the n-th prime, where n is a C-int
        """
        global num_primes

        if n <= 0:
            raise ValueError, "nth prime meaningless for non-positive n (=%s)"%n
        cdef GEN g
        pari_catch_sig_on()
        g = prime(n)
        return self.new_gen(g)


    def nth_prime(self, long n):
        try:
            return self.__nth_prime(n)
        except PariError:
            self.init_primes(max(2*num_primes,20*n))
            return self.nth_prime(n)

    def euler(self, precision=0):
        """
        Return Euler's constant to the requested real precision (in bits).

        EXAMPLES::

            sage: pari.euler()
            0.577215664901533
            sage: pari.euler(precision=100).python()
            0.577215664901532860606512090082...
        """
        pari_catch_sig_on()
        return self.new_gen(mpeuler(pbw(precision)))

    def pi(self, precision=0):
        """
        Return the value of the constant pi to the requested real precision
        (in bits).

        EXAMPLES::

            sage: pari.pi()
            3.14159265358979
            sage: pari.pi(precision=100).python()
            3.1415926535897932384626433832...
        """
        pari_catch_sig_on()
        return self.new_gen(mppi(pbw(precision)))

    def pollegendre(self, long n, v=-1):
        """
        pollegendre(n, v=x): Legendre polynomial of degree n (n C-integer),
        in variable v.

        EXAMPLES::

            sage: pari.pollegendre(7)
            429/16*x^7 - 693/16*x^5 + 315/16*x^3 - 35/16*x
            sage: pari.pollegendre(7, 'z')
            429/16*z^7 - 693/16*z^5 + 315/16*z^3 - 35/16*z
            sage: pari.pollegendre(0)
            1
        """
        pari_catch_sig_on()
        return self.new_gen(pollegendre(n, self.get_var(v)))

    def poltchebi(self, long n, v=-1):
        """
        poltchebi(n, v=x): Chebyshev polynomial of the first kind of degree
        n, in variable v.

        EXAMPLES::

            sage: pari.poltchebi(7)
            64*x^7 - 112*x^5 + 56*x^3 - 7*x
            sage: pari.poltchebi(7, 'z')
            64*z^7 - 112*z^5 + 56*z^3 - 7*z
            sage: pari.poltchebi(0)
            1
        """
        pari_catch_sig_on()
        return self.new_gen(polchebyshev1(n, self.get_var(v)))

    def factorial(self, long n):
        """
        Return the factorial of the integer n as a PARI gen.

        EXAMPLES::

            sage: pari.factorial(0)
            1
            sage: pari.factorial(1)
            1
            sage: pari.factorial(5)
            120
            sage: pari.factorial(25)
            15511210043330985984000000
        """
        pari_catch_sig_on()
        return self.new_gen(mpfact(n))

    def polcyclo(self, long n, v=-1):
        """
        polcyclo(n, v=x): cyclotomic polynomial of degree n, in variable
        v.

        EXAMPLES::

            sage: pari.polcyclo(8)
            x^4 + 1
            sage: pari.polcyclo(7, 'z')
            z^6 + z^5 + z^4 + z^3 + z^2 + z + 1
            sage: pari.polcyclo(1)
            x - 1
        """
        pari_catch_sig_on()
        return self.new_gen(polcyclo(n, self.get_var(v)))

    def polcyclo_eval(self, long n, v):
        """
        polcyclo_eval(n, v): value of the nth cyclotomic polynomial at value v.

        EXAMPLES::

            sage: pari.polcyclo_eval(8, 2)
            17
            sage: cyclotomic_polynomial(8)(2)
            17
        """
        t0GEN(v)
        pari_catch_sig_on()
        return self.new_gen(polcyclo_eval(n, t0))

    def polsubcyclo(self, long n, long d, v=-1):
        """
        polsubcyclo(n, d, v=x): return the pari list of polynomial(s)
        defining the sub-abelian extensions of degree `d` of the
        cyclotomic field `\QQ(\zeta_n)`, where `d`
        divides `\phi(n)`.

        EXAMPLES::

            sage: pari.polsubcyclo(8, 4)
            [x^4 + 1]
            sage: pari.polsubcyclo(8, 2, 'z')
            [z^2 - 2, z^2 + 1, z^2 + 2]
            sage: pari.polsubcyclo(8, 1)
            [x - 1]
            sage: pari.polsubcyclo(8, 3)
            []
        """
        cdef gen plist
        pari_catch_sig_on()
        plist = self.new_gen(polsubcyclo(n, d, self.get_var(v)))
        if typ(plist.g) != t_VEC:
            return pari.vector(1, [plist])
        else:
            return plist
        #return self.new_gen(polsubcyclo(n, d, self.get_var(v)))

    def polzagier(self, long n, long m):
        pari_catch_sig_on()
        return self.new_gen(polzag(n, m))

    def setrand(self, seed):
        """
        Sets PARI's current random number seed.

        INPUT:

        - ``seed`` -- either a strictly positive integer or a GEN of
          type ``t_VECSMALL`` as output by ``getrand()``

        This should not be called directly; instead, use Sage's global
        random number seed handling in ``sage.misc.randstate``
        and call ``current_randstate().set_seed_pari()``.

        EXAMPLES::

            sage: pari.setrand(50)
            sage: a = pari.getrand(); a
            Vecsmall([...])
            sage: pari.setrand(a)
            sage: a == pari.getrand()
            True

        TESTS:

        Check that invalid inputs are handled properly (#11825)::

            sage: pari.setrand(0)
            Traceback (most recent call last):
            ...
            PariError: incorrect type (11)
            sage: pari.setrand("foobar")
            Traceback (most recent call last):
            ...
            PariError: incorrect type (11)
        """
        t0GEN(seed)
        pari_catch_sig_on()
        setrand(t0)
        pari_catch_sig_off()

    def getrand(self):
        """
        Returns PARI's current random number seed.

        OUTPUT:

        GEN of type t_VECSMALL

        EXAMPLES::

            sage: pari.setrand(50)
            sage: a = pari.getrand(); a
            Vecsmall([...])
            sage: pari.setrand(a)
            sage: a == pari.getrand()
            True
        """
        pari_catch_sig_on()
        return self.new_gen(getrand())

    def vector(self, long n, entries=None):
        """
        vector(long n, entries=None): Create and return the length n PARI
        vector with given list of entries.

        EXAMPLES::

            sage: pari.vector(5, [1, 2, 5, 4, 3])
            [1, 2, 5, 4, 3]
            sage: pari.vector(2, [x, 1])
            [x, 1]
            sage: pari.vector(2, [x, 1, 5])
            Traceback (most recent call last):
            ...
            IndexError: length of entries (=3) must equal n (=2)
        """
        cdef gen v = self._empty_vector(n)
        if entries is not None:
            if len(entries) != n:
                raise IndexError, "length of entries (=%s) must equal n (=%s)"%\
                      (len(entries), n)
            for i, x in enumerate(entries):
                v[i] = x
        return v

    cdef gen _empty_vector(self, long n):
        cdef gen v
        pari_catch_sig_on()
        v = self.new_gen(zerovec(n))
        return v

    def matrix(self, long m, long n, entries=None):
        """
        matrix(long m, long n, entries=None): Create and return the m x n
        PARI matrix with given list of entries.
        """
        cdef long i, j, k
        cdef gen A
        cdef gen x

        pari_catch_sig_on()
         # The gtomat is very important!!  Without sage/PARI will segfault.
         # I do not know why. -- William Stein
        A = self.new_gen(gtomat(zeromat(m,n)))
        if entries is not None:
            if len(entries) != m*n:
                raise IndexError, "len of entries (=%s) must be %s*%s=%s"%(len(entries),m,n,m*n)
            k = 0
            for i from 0 <= i < m:
                for j from 0 <= j < n:
                    x = pari(entries[k])
                    A._refers_to[(i,j)] = x
                    (<GEN>(A.g)[j+1])[i+1] = <long>(x.g)
                    k = k + 1
        return A


##############################################
# Used in integer factorization -- must be done
# after the pari_instance creation above:

cdef gen _tmp = P.new_gen_noclear(gp_read_str('1000000000000000'))
cdef GEN ten_to_15 = _tmp.g

##############################################

def init_pari_stack(size=8000000):
    """
    Change the PARI scratch stack space to the given size.

    The main application of this command is that you've done some
    individual PARI computation that used a lot of stack space. As a
    result the PARI stack may have doubled several times and is now
    quite large. That object you computed is copied off to the heap,
    but the PARI stack is never automatically shrunk back down. If you
    call this function you can shrink it back down.

    If you set this too small then it will automatically be increased
    if it is exceeded, which could make some calculations initially
    slower (since they have to be redone until the stack is big
    enough).

    INPUT:


    -  ``size`` - an integer (default: 8000000)


    EXAMPLES::

        sage: get_memory_usage()                       # random output
        '122M+'
        sage: a = pari('2^100000000')
        sage: get_memory_usage()                       # random output
        '157M+'
        sage: del a
        sage: get_memory_usage()                       # random output
        '145M+'

    Hey, I want my memory back!

    ::

        sage: sage.libs.pari.gen.init_pari_stack()
        sage: get_memory_usage()                       # random output
        '114M+'

    Ahh, that's better.
    """
    init_stack(size)

cdef int init_stack(size_t size) except -1:
    cdef size_t s
    cdef pari_sp cur_stack_size

    global top, bot, avma, mytop


    err = False    # whether or not a memory allocation error occurred.


    # delete this if get core dumps and change the 2* to a 1* below.
    if bot:
        sage_free(<void*>bot)

    prev_stack_size = top - bot
    if size == 0:
        size = 2 * prev_stack_size

    # Decide on size
    s = fix_size(size)

    # Allocate memory for new stack using Python's memory allocator.
    # As explained in the python/C API reference, using this instead
    # of malloc is much better (and more platform independent, etc.)
    bot = <pari_sp> sage_malloc(s)

    while not bot:
        err = True
        s = fix_size(prev_stack_size)
        bot = <pari_sp> sage_malloc(s)
        if not bot:
            prev_stack_size /= 2

    top = bot + s
    mytop = top
    avma = top

    if err:
        raise MemoryError, "Unable to allocate %s bytes memory for PARI."%size


cdef size_t fix_size(size_t a):
    cdef size_t b
    b = a - (a & (sizeof(long)-1))     # sizeof(long) | b <= a
    if b < 1024:
        b = 1024
    return b

cdef GEN deepcopy_to_python_heap(GEN x, pari_sp* address):
    cdef size_t s = <size_t> gsizebyte(x)
    cdef pari_sp tmp_bot, tmp_top

    tmp_bot = <pari_sp> sage_malloc(s)
    tmp_top = tmp_bot + s
    address[0] = tmp_bot
    return gcopy_avma(x, &tmp_top)

cdef gen _new_gen (GEN x):
    cdef GEN h
    cdef pari_sp address
    cdef gen y
    h = deepcopy_to_python_heap(x, &address)
    y = PY_NEW(gen)
    y.init(h, address)
    return y

cdef GEN _Vec_append(GEN v, GEN a, long n):
    """
    This implements appending zeros (or another constant GEN ``a``) to
    the result of :meth:`Vec` and similar functions.

    This is a shallow function, copying ``a`` and entries of ``v`` to
    the result.  The result is simply stored on the PARI stack.

    INPUT:

    - ``v`` -- GEN of type ``t_VEC`` or ``t_COL``

    - ``a`` -- GEN which will be used for the added entries.
      Normally, this would be ``gen_0``.

    - ``n`` -- Make the vector of minimal length `|n|`. If `n > 0`,
      append zeros; if `n < 0`, prepend zeros.

    OUTPUT:

    A GEN of the same type as ``v``.
    """
    cdef long lenv = lg(v)-1
    cdef GEN w
    cdef long i
    # Do we need to extend the vector with zeros?
    if n > lenv:
        w = cgetg(n+1, typ(v))
        for i from 1 <= i <= lenv:
            set_gel(w, i, gel(v, i))
        for i from 1 <= i <= n-lenv:
            set_gel(w, i+lenv, a)
        return w
    elif n < -lenv:
        n = -n  # Make n positive
        w = cgetg(n+1, typ(v))
        for i from 1 <= i <= lenv:
            set_gel(w, i+(n-lenv), gel(v, i))
        for i from 1 <= i <= n-lenv:
            set_gel(w, i, a)
        return w
    else:
        return v


#######################
# Base gen class
#######################


cdef extern from "pari/pari.h":
    char *errmessage[]
    int talker2, bugparier, alarmer, openfiler, talker, flagerr, impl, \
        archer, notfuncer, precer, typeer, consister, user, errpile, \
        overflower, matinv1, mattype1, arither1, primer1, invmoder, \
        constpoler, notpoler, redpoler, zeropoler, operi, operf, gdiver, \
        memer, negexper, sqrter5, noer
    int warner, warnprec, warnfile, warnmem

cdef extern from "misc.h":
    int     factorint_withproof_sage(GEN* ans, GEN x, GEN cutoff)
    int     gcmp_sage(GEN x, GEN y)

def __errmessage(d):
    if d <= 0 or d > noer:
        return "unknown"
    return errmessage[d]

# FIXME: we derive PariError from RuntimeError, for backward
# compatibility with code that catches the latter. Once this is
# in production, we should change the base class to StandardError.
from exceptions import RuntimeError

# can we have "cdef class" ?
# because of the inheritance, need to somehow "import" the built-in
# exception class...
class PariError (RuntimeError):

    errmessage = staticmethod(__errmessage)

    def errnum(self):
        r"""
        Return the PARI error number corresponding to this exception.

        EXAMPLES::

            sage: try:
            ...     pari('1/0')
            ... except PariError, err:
            ...     print err.errnum()
            27
        """
        return self.args[0]

    def __repr__(self):
        r"""
        TESTS::

            sage: PariError(11)
            PariError(11)
        """
        return "PariError(%d)"%self.errnum()

    def __str__(self):
        r"""
        EXAMPLES::

            sage: try:
            ...     pari('1/0')
            ... except PariError, err:
            ...     print err
            division by zero (27)
        """
        return "%s (%d)"%(self.errmessage(self.errnum()), self.errnum())


# We expose a trap function to C.
# If this function returns without raising an exception,
# the code is retried.
# This is a proof-of-concept example.
# THE TRY CODE IS NOT REENTRANT -- NO CALLS TO PARI FROM HERE !!!
#              - Gonzalo Tornario

cdef void _pari_trap "_pari_trap"(long errno, long retries) except *:
    if retries > 100:
        pari_catch_sig_off()
        raise RuntimeError("_pari_trap recursion too deep")
    if errno == errpile:
        P.allocatemem(silent=True)
    elif errno == user:
        pari_catch_sig_off()
        raise RuntimeError("PARI user exception")
    else:
        pari_catch_sig_off()
        raise PariError(errno)


def vecsmall_to_intlist(gen v):
    """
    INPUT:


    -  ``v`` - a gen of type Vecsmall


    OUTPUT: a Python list of Python ints
    """
    if typ(v.g) != t_VECSMALL:
        raise TypeError, "input v must be of type vecsmall (use v.Vecsmall())"
    return [v.g[k+1] for k in range(glength(v.g))]



cdef _factor_int_when_pari_factor_failed(x, failed_factorization):
    """
    This is called by factor when PARI's factor tried to factor, got
    the failed_factorization, and it turns out that one of the factors
    in there is not proved prime. At this point, we don't care too much
    about speed (so don't write everything below using the PARI C
    library), since the probability this function ever gets called is
    infinitesimal. (That said, we of course did test this function by
    forcing a fake failure in the code in misc.h.)
    """
    P = failed_factorization[0]  # 'primes'
    E = failed_factorization[1]  # exponents
    if len(P) == 1 and E[0] == 1:
        # Major problem -- factor can't split the integer at all, but it's composite.  We're stuffed.
        print "BIG WARNING: The number %s wasn't split at all by PARI, but it's definitely composite."%(P[0])
        print "This is probably an infinite loop..."
    w = []
    for i in range(len(P)):
        p = P[i]
        e = E[i]
        if not p.isprime():
            # Try to factor further -- assume this works.
            F = p.factor(proof=True)
            for j in range(len(F[0])):
                w.append((F[0][j], F[1][j]))
        else:
            w.append((p, e))
    m = pari.matrix(len(w), 2)
    for i in range(len(w)):
        m[i,0] = w[i][0]
        m[i,1] = w[i][1]
    return m

