r"""
Rational Numbers

AUTHORS:
    -- William Stein (2005): first version
    -- William Stein (2006-02-22): floor and ceil (pure fast GMP versions).
    -- Gonzalo Tornaria and William Stein (2006-03-02): greatly improved
                    python/GMP conversion; hashing
    -- William Stein and Naqi Jaffery (2006-03-06): height, sqrt examples,
          and improve behavior of sqrt.
    -- David Harvey (2006-09-15): added nth_root
    -- Pablo De Napoli (2007-04-01): corrected the implementations of
       multiplicative_order, is_one; optimized __nonzero__ ; documented: lcm,gcd

TESTS:
    sage: a = -2/3
    sage: a == loads(dumps(a))
    True
"""


###########################################################################
#       Copyright (C) 2004, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

include "../ext/interrupt.pxi"  # ctrl-c interrupt block support
include "../ext/gmp.pxi"
include "../ext/stdsage.pxi"


import sys
import operator

from sage.misc.mathml import mathml

import sage.misc.misc as misc
import sage.rings.rational_field
import sage.libs.pari.all

cimport integer
import integer

from integer_ring import ZZ

from sage.structure.element cimport Element, RingElement, ModuleElement
from sage.structure.element import bin_op
from sage.categories.morphism cimport Morphism

import sage.rings.real_mpfr

cimport sage.ext.arith
import  sage.ext.arith

cdef sage.ext.arith.arith_int ai
ai = sage.ext.arith.arith_int()

cdef extern from "mpz_pylong.h":
    cdef mpz_get_pylong(mpz_t src)
    cdef int mpz_set_pylong(mpz_t dst, src) except -1
    cdef long mpz_pythonhash(mpz_t src)

cdef class Rational(sage.structure.element.FieldElement)

cdef public void set_from_mpq(Rational self, mpq_t value):
    mpq_set(self.value, value)

cdef public void set_from_Rational(Rational self, Rational other):
    mpq_set(self.value, other.value)

cdef public void set_from_Integer(Rational self, integer.Integer other):
    mpq_set_z(self.value, other.value)

cdef object Rational_mul_(Rational a, Rational b):
    cdef Rational x
    x = <Rational> PY_NEW(Rational)

    _sig_on
    mpq_mul(x.value, a.value, b.value)
    _sig_off

    return x

cdef object Rational_div_(Rational a, Rational b):
    cdef Rational x
    x = <Rational> PY_NEW(Rational)

    _sig_on
    mpq_div(x.value, a.value, b.value)
    _sig_off

    return x

cdef Rational_add_(Rational self, Rational other):
    cdef Rational x
    x = <Rational> PY_NEW(Rational)
    _sig_on
    mpq_add(x.value, self.value, other.value)
    _sig_off
    return x

cdef Rational_sub_(Rational self, Rational other):
    cdef Rational x
    x = <Rational> PY_NEW(Rational)

    _sig_on
    mpq_sub(x.value, self.value, other.value)
    _sig_off

    return x

cdef object the_rational_ring
the_rational_ring = sage.rings.rational_field.Q

# make sure zero/one elements are set
cdef set_zero_one_elements():
    global the_rational_ring
    the_rational_ring._zero_element = Rational(0)
    the_rational_ring._one_element = Rational(1)

set_zero_one_elements()

cdef class Rational(sage.structure.element.FieldElement):
    """
    A Rational number.

    Rational numbers are implemented using the GMP C library.

    EXAMPLES:
        sage: a = -2/3
        sage: type(a)
        <type 'sage.rings.rational.Rational'>
        sage: parent(a)
        Rational Field
        sage: Rational('1/0')
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 1/0 to a rational
        sage: Rational(1.5)
        3/2
        sage: Rational('9/6')
        3/2
        sage: Rational((2^99,2^100))
        1/2
        sage: Rational(("2", "10"), 16)
        1/8
        sage: Rational(QQbar(125/8).nth_root(3))
        5/2
        sage: Rational(AA(209735/343 - 17910/49*golden_ratio).nth_root(3) + 3*golden_ratio)
        53/7

    Conversion from PARI:
        sage: Rational(pari('-939082/3992923'))
        -939082/3992923
    """
    def __new__(self, x=None, int base=0):
        global the_rational_ring
        mpq_init(self.value)
        self._parent = the_rational_ring

    def __init__(self, x=None, unsigned int base=0):
        """
        Create a new rational number.

        INPUT:
            x -- object (default: None)
            base -- base if x is a string

        EXAMPLES:
            sage: a = Rational()
            sage: a.__init__(7); a
            7
            sage: a.__init__('70', base=8); a
            56
            sage: a.__init__(pari('2/3')); a
            2/3
        """
        if x is not None:
            self.__set_value(x, base)

    def __reduce__(self):
        """
        Used in pickling rational numbers.

        EXAMPLES:
            sage: a = 3/5
            sage: a.__reduce__()
            (<built-in function make_rational>, ('3/5',))
        """
        return sage.rings.rational.make_rational, (self.str(32),)

    def __index__(self):
        """
        Needed so integers can be used as list indices.

        EXAMPLES:
            sage: v = [1,2,3,4,5]
            sage: v[3/1]
            4
            sage: v[3/2]
            Traceback (most recent call last):
            ...
            TypeError: rational is not an integer
        """
        if self.denominator() == 1:
            return int(self)
        raise TypeError, "rational is not an integer"

    def _reduce_set(self, s):
        """
        Used in setting a rational number when unpickling.  Do not
        call this from external code since it violates immutability.

        INPUT:
            s -- string representation of rational in base 32

        EXAMPLES:
            sage: a = -17/3730; _, (s,) = a.__reduce__(); s
            '-h/3ki'
            sage: b = 2/3; b._reduce_set('-h/3ki'); b
            -17/3730
        """
        mpq_set_str(self.value, s, 32)

    def __set_value(self, x, unsigned int base):
        """
        Set the value of this rational number.  This function
        handles numerous cases.

        INPUT:
            x -- object
            base -- positive integer

        EXAMPLES:
            sage: a = 3/5
            sage: a.__set_value('-h/3ki', 32); a
            -17/3730
        """
        cdef int n
        cdef Rational temp_rational
        cdef integer.Integer a, b

        if isinstance(x, Rational):
            set_from_Rational(self, x)

        elif isinstance(x, int):
            i = x
            mpq_set_si(self.value, i, 1)

        elif isinstance(x, long):
            mpz_set_pylong(mpq_numref(self.value), x)

        elif isinstance(x, integer.Integer):
            set_from_Integer(self, x)

        elif isinstance(x, sage.rings.real_mpfr.RealNumber):

            if x == 0:
                mpq_set_si(self.value, 0, 1)
                return
            if not base:
                set_from_Rational(self, x.simplest_rational())
            else:
                xstr = x.str(base)
                if '.' in xstr:
                    exp = (len(xstr) - (xstr.index('.') +1))
                    p = base**exp
                    pstr = '1'+'0'*exp
                    s = xstr.replace('.','') +'/'+pstr
                    n = mpq_set_str( self.value, s, base)
                    if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                        raise TypeError, "unable to convert %s to a rational"%x
                    mpq_canonicalize(self.value)
                else:
                    n = mpq_set_str(self.value, xstr, base)
                    if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                        raise TypeError, "unable to convert %s to a rational"%x
                    mpq_canonicalize(self.value)

        elif isinstance(x, str):
            n = mpq_set_str(self.value, x, base)
            if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                raise TypeError, "unable to convert %s to a rational"%x
            mpq_canonicalize(self.value)

        elif hasattr(x, "_rational_"):
            set_from_Rational(self, x._rational_())

        elif isinstance(x, tuple) and len(x) == 2:
            num = x[0]
            denom = x[1]
            if isinstance(num, int) and isinstance(denom, int):
                mpq_set_si(self.value, num, denom)
            else:
                if not isinstance(num, integer.Integer):
                    num = integer.Integer(num, base)
                if not isinstance(denom, integer.Integer):
                    denom = integer.Integer(denom, base)
                mpz_set(mpq_numref(self.value), (<integer.Integer>num).value)
                mpz_set(mpq_denref(self.value), (<integer.Integer>denom).value)
            if mpz_sgn(mpq_denref(self.value)) == 0:
                raise ValueError, "denominator must not be 0"
            mpq_canonicalize(self.value)

        elif isinstance(x, list) and len(x) == 1:
            self.__set_value(x[0], base)

        elif isinstance(x, sage.libs.pari.all.pari_gen):
            a = integer.Integer(x.numerator())
            b = integer.Integer(x.denominator())
            mpz_set(mpq_numref(self.value), a.value)
            mpz_set(mpq_denref(self.value), b.value)

        elif hasattr(x, 'rational_reconstruction'):
            temp_rational = x.rational_reconstruction()
            mpq_set(self.value, temp_rational.value)

        else:

            raise TypeError, "Unable to coerce %s (%s) to Rational"%(x,type(x))

    cdef void set_from_mpq(Rational self, mpq_t value):
        mpq_set(self.value, value)

    def list(self):
        """
        Return a list with the rational element in it, to be
        compatible with the method for number fields.

        OUTPUT:
            list -- the list [self]

        EXAMPLES:
            sage: m = 5/3
            sage: m.list()
            [5/3]
        """
        return [ self ]


    def __richcmp__(left, right, int op):
        """
        INPUT:
            left, right -- objects
            op -- integer

        EXAMPLES:
            sage: 1/3 < 2/3
            True
            sage: 2/3 < 1/3
            False
            sage: 4/5 < 2.0
            True
            sage: 4/5 < 0.8
            False
        """
        return (<sage.structure.element.Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, sage.structure.element.Element right) except -2:
        cdef int i
        i = mpq_cmp((<Rational>left).value, (<Rational>right).value)
        if i < 0: return -1
        elif i == 0: return 0
        else: return 1

    def copy(self):
        """
        Return a copy of self.

        OUTPUT:
            Rational

        EXAMPLES:
            sage: a = -17/37
            sage: a.copy() is a
            False

        Coercion does not make a new copy:
            sage: QQ(a) is a
            True

        The constructor also makes a new copy:
            sage: Rational(a) is a
            False
        """
        cdef Rational z
        z = <Rational> PY_NEW(Rational)
        mpq_set(z.value, self.value)
        return z

    def  __dealloc__(self):
        """
        Free memory occupied by this rational number.

        EXAMPLES:
            sage: a = -17/37
            sage: del a          # indirect test
        """
        mpq_clear(self.value)

    def __repr__(self):
        """
        Return string representation of this rational number.

        EXAMPLES:
            sage: a = -17/37; a.__repr__()
            '-17/37'
        """
        return self.str()

    def _latex_(self):
        """
        Return Latex representation of this rational number.

        EXAMPLES:
            sage: a = -17/37
            sage: a._latex_()
            '-\\frac{17}{37}'
        """
        if self.denom() == 1:
            return str(self.numer())
        else:
            if self < 0:
                return "-\\frac{%s}{%s}"%(-self.numer(), self.denom())
            else:
               return "\\frac{%s}{%s}"%(self.numer(), self.denom())

    def _sympy_(self):
        """
        Convert Sage Rational to SymPy Rational.

        EXAMPLES:
            sage: n = 1/2; n._sympy_()
            1/2
            sage: n = -1/5; n._sympy_()
            -1/5
            sage: from sympy import Symbol
            sage: QQ(1)+Symbol('x')*QQ(2)
            2*x + 1

        """
        import sympy
        return sympy.Rational(int(self.numerator()), int(self.denominator()))


    def _mathml_(self):
        """
        Return mathml representation of this rational number.

        EXAMPLES:
            sage: a = -17/37; a._mathml_()
            '<mo>-</mo><mfrac><mrow><mn>17</mn></mrow><mrow><mn>37</mn></mrow></mfrac>'
        """
        if self.denom() == 1:
            return '<mn>%s</mn>'%(self.numer())
        else:
            t = ''
            if self < 0:
                t = t + '<mo>-</mo>'
            t = t + '<mfrac><mrow>%s</mrow><mrow>%s</mrow></mfrac>'%(
                mathml(abs(self.numer())), mathml(self.denom()))
            return t

    def _im_gens_(self, codomain, im_gens):
        """
        Return the image of self under the homomorphism from the
        rational field to codomain.  This always just returns self
        coerced into the codomain.

        INPUT:
            codomain -- object (usually a ring)
            im_gens -- list of elements of codomain

        EXAMPLES:
            sage: a = -17/37
            sage: a._im_gens_(QQ, [1/1])
            -17/37
        """
        return codomain._coerce_(self)

    def lcm(self, Rational other):
        r"""
        Return the least common multiple of self and other.

        One way to define this notion is the following:

        Note that each rational positive rational number can be written
        as a product of primes with integer (positive or negative)
        powers in a unique way.

        Then, the LCM of two rational numbers x,y can be defined by
        specifying that the  exponent of every prime p in lcm(x,y)
        is the supremum of the exponents of p in x,
        and the exponent of p  in y
        (The primes that does not appear in the decomposition of x
        or y are considered to have exponent zero).

        This definition  is consistent with the definition of the LCM
        in the rational integers. Our hopefully interesting notion of LCM
        for rational numbers is illustrated in the examples below.

        EXAMPLES:
            sage: lcm(2/3,1/5)
            2

        This is consistent with the definition above, since:
        $$
           2/3 = 2^1 * 3^{-1}*5^0
        $$
        $$
           1/5 = 2^0 * 3^0   *5^{-1}
        $$
        and hence,
        $$
           lcm(2/3,1/5)= 2^1*3^0*5^0 = 2.
        $$

            sage: lcm(2/3,7/5)
            14

        In this example:
        $$
           2/3 = 2^1*3^{-1}*5^0    * 7^0
        $$
        $$
           7/5 = 2^0*3^0   *5^{-1} * 7^1
        $$
        $$
           lcm(2/3,7/5) = 2^1*3^0*5^0*7^1 = 14
        $$

            sage: lcm(1/3,1/5)
            1

        In this example:
        $$
           1/3 = 3^{-1}*5^0
        $$
        $$
           1/5 = 3^0 * 5^{-1}
        $$
        $$
           lcm(1/3,1/5)=3^0*5^0=1
        $$

            sage: lcm(1/3,1/6)
            1/3

        In this example:
        $$
           1/3 = 2^0*3^{-1}
        $$
        $$
           1/6 = 2^{-1}*3^{-1}
        $$
        $$
           lcm(1/3,1/6)=2^0*3^{-1}=1/3
        $$
        """
        d = self.denom()*other.denom()
        self_d = self.numer()*other.denom()
        other_d = other.numer()*self.denom()
        return self_d.lcm(other_d) / d

    def gcd(self, Rational other):
        """
        Return the least common multiple of self and other.

        One way to define this notion is the following:

        Note that each rational positive rational number can be written
        as a product of primes  with integer (positive or negative)
        powers in a unique way.

        Then, the GCD of two rational numbers x,y can be defined by
        specifying that the exponent of every prime p in gcd(x,y) is
        the infimum of the exponents of p in x,
        and  the exponent of p  in y
        (The primes that does not appear in the decomposition of x or y
        are considered to have exponent zero).

        This definition is consistent with the definition of the GCD
        in the rational integers.  Our hopefully interesting notion of GCD
        for rational numbers is illustrated in the examples below.

        EXAMPLES:

            sage: gcd(2/3,1/5)
            1/15

        This is consistent with the definition above, since:
        $$
           2/3 = 2^1 * 3^{-1}*5^0
        $$
        $$
           1/5 = 2^0 * 3^0   *5^{-1}
        $$
        and hence,
        $$
           gcd(2/3,1/5)= 2^0*3^{-1}*5^{-1} = 1/15
        $$

            sage: gcd(2/3,7/5)
            1/15

        In this example:
        $$
          2/3 = 2^1*3^{-1}*5^0    * 7^0
        $$
        $$
          7/5 = 2^0*3^0   *5^{-1} * 7^1
        $$
        $$
          gcd(2/3,7/5) = 2^0*3^{-1}*5^{-1}*7^0 = 1/15
        $$

            sage: gcd(1/3,1/6)
            1/6

        In this example:
        $$
          1/3 = 2^0*3^{-1}
        $$
        $$
          1/6 = 2^{-1}*3^{-1}
        $$
        $$
          gcd(1/3,1/6)=2^{-1}*3^{-1}=1/6
        $$

            sage: gcd(6/7,9/7)
            3/7

        In this example:
        $$
          6/7 = 2^1*3^1*7^{-1}
        $$
        $$
          9/7 = 2^0*3^2*7^{-1}
        $$
        $$
          gcd(6/7,9/7)=2^0*3^1*7^{-1}=3/7
        $$
        """
        d = self.denom()*other.denom()
        self_d = self.numer()*other.denom()
        other_d = other.numer()*self.denom()
        return self_d.gcd(other_d) / d

    def valuation(self, p):
        """
        Return the largest power of p that divides self.

        INPUT:
            p -- a prime number

        EXAMPLES:
            sage: x = -5/9
            sage: x.valuation(5)
            1
            sage: x.valuation(3)
            -2
            sage: x.valuation(2)
            0

        Some edge cases:
            sage: (0/1).valuation(4)
            +Infinity
            sage: (7/16).valuation(4)
            -2
        """
        return self.numerator().valuation(p) - self.denominator().valuation(p)

    def is_square(self):
        """
        Return whether or not this rational number is a square.

        OUTPUT:
            bool

        EXAMPLES:
            sage: x = 9/4
            sage: x.is_square()
            True
            sage: x = (7/53)^100
            sage: x.is_square()
            True
            sage: x = 4/3
            sage: x.is_square()
            False
            sage: x = -1/4
            sage: x.is_square()
            False
        """
        return mpq_sgn(self.value) >= 0 and mpz_perfect_square_p(mpq_numref(self.value)) and mpz_perfect_square_p(mpq_denref(self.value))

    def squarefree_part(self):
        """
        Return the square free part of $x$, i.e., an integer z such that $x = z y^2$,
        for a perfect square $y^2$.

        EXAMPLES:
            sage: a = 1/2
            sage: a.squarefree_part()
            2
            sage: b = a/a.squarefree_part()
            sage: b, b.is_square()
            (1/4, True)
            sage: a = 24/5
            sage: a.squarefree_part()
            30
        """
        return self.numer().squarefree_part() * self.denom().squarefree_part()

    def sqrt_approx(self, prec=None, all=False):
        """
        Return numerical approximation with given number of bits
        of precision to this rational number.  If all is given,
        return both approximations.

        INPUT:
            prec -- integer
            all -- bool

        EXAMPLES:
            sage: (5/3).sqrt_approx()
            1.29099444873581
            sage: (990829038092384908234098239048230984/4).sqrt_approx()
            4.9770197862083713747374920870362581922510725585130996993055116540856385e17
            sage: (5/3).sqrt_approx(prec=200)
            1.2909944487358056283930884665941332036109739017638636088625
            sage: (9/4).sqrt_approx()
            3/2
        """
        try:
            return self.sqrt(extend=False,all=all)
        except ValueError:
            pass
        if prec is None:
            prec = max(max(53, 2*(mpz_sizeinbase(mpq_numref(self.value), 2)+2)),
                   2*(mpz_sizeinbase(mpq_denref(self.value), 2)+2))
        return self.sqrt(prec=prec, all=all)

    def val_unit(self, p):
        r"""
        Returns a pair: the p-adic valuation of self, and the p-adic
        unit of self, as a Rational.

        We do not require the p be prime, but it must be at least 2.
        For more documentation see \code{Integer.val_unit}

        INPUT:
            p -- a prime
        OUTPUT:
            int -- the p-adic valuation of this rational
            Rational -- p-adic unit part of self

        EXAMPLES:
            sage: (-4/17).val_unit(2)
            (2, -1/17)
            sage: (-4/17).val_unit(17)
            (-1, -4)
            sage: (0/1).val_unit(17)
            (+Infinity, 1)

        AUTHOR:
            -- David Roe (4/12/07)
        """
        return self._val_unit(p)

    # TODO -- change to use cpdef?  If so, must fix
    # code in padics, etc.  Do search_src('_val_unit').
    cdef _val_unit(Rational self, integer.Integer p):
        """
        This is called by self.val_unit.
        """
        cdef integer.Integer v
        cdef Rational u
        if mpz_cmp_ui(p.value, 2) < 0:
            raise ValueError, "p must be at least 2."
        if mpq_sgn(self.value) == 0:
            import sage.rings.infinity
            u = PY_NEW(Rational)
            mpq_set_ui(u.value, 1, 1)
            return (sage.rings.infinity.infinity, u)
        v = PY_NEW(integer.Integer)
        u = PY_NEW(Rational)
        _sig_on
        mpz_set_ui(v.value, mpz_remove(mpq_numref(u.value), mpq_numref(self.value), p.value))
        _sig_off
        if mpz_sgn(v.value) != 0:
            mpz_set(mpq_denref(u.value), mpq_denref(self.value))
        else:
            _sig_on
            mpz_set_ui(v.value, mpz_remove(mpq_denref(u.value), mpq_denref(self.value), p.value))
            _sig_off
            mpz_neg(v.value, v.value)
        return (v, u)

    def sqrt(self, prec=None, extend=True, all=False):
        r"""
        The square root function.

        INPUT:
            prec -- integer (default: None): if None, returns an exact
                 square root; otherwise returns a numerical square
                 root if necessary, to the given bits of precision.
            extend -- bool (default: True); if True, return a square
                 root in an extension ring, if necessary. Otherwise,
                 raise a ValueError if the square is not in the base
                 ring.
            all -- bool (default: False); if True, return all square
                   roots of self, instead of just one.

        EXAMPLES:
            sage: x = 25/9
            sage: x.sqrt()
            5/3
            sage: x = 64/4
            sage: x.sqrt()
            4
            sage: x = 100/1
            sage: x.sqrt()
            10
            sage: x.sqrt(all=True)
            [10, -10]
            sage: x = 81/5
            sage: x.sqrt()
            9/sqrt(5)
            sage: x = -81/3
            sage: x.sqrt()
            3*sqrt(3)*I

            sage: n = 2/3
            sage: n.sqrt()
            sqrt(2)/sqrt(3)
            sage: n.sqrt(prec=10)
            0.82
            sage: n.sqrt(prec=100)
            0.81649658092772603273242802490
            sage: n.sqrt(prec=100)^2
            0.66666666666666666666666666667
            sage: n.sqrt(prec=53, all=True)
            [0.816496580927726, -0.816496580927726]
            sage: n.sqrt(extend=False, all=True)
            Traceback (most recent call last):
            ...
            ValueError: square root of 2/3 not a rational number
            sage: sqrt(-2/3, all=True)
            [sqrt(2)*I/sqrt(3), -sqrt(2)*I/sqrt(3)]
            sage: sqrt(-2/3, prec=53)
            0.816496580927726*I
            sage: sqrt(-2/3, prec=53, all=True)
            [0.816496580927726*I, -0.816496580927726*I]

        AUTHOR:
            -- Naqi Jaffery (2006-03-05): some examples
        """
        if mpq_sgn(self.value) == 0:
            return [self] if all else self

        if mpq_sgn(self.value) < 0:
            if not extend:
                raise ValueError, "square root of negative number not rational"
            if prec:
                from sage.rings.complex_field import ComplexField
                K = ComplexField(prec)
                return K(self).sqrt(all=all)
            from sage.calculus.calculus import sqrt
            return sqrt(self, all=all)

        cdef Rational z = <Rational> PY_NEW(Rational)
        cdef mpz_t tmp
        cdef int non_square = 0

        _sig_on
        mpz_init(tmp)
        mpz_sqrtrem(mpq_numref(z.value), tmp, mpq_numref(self.value))
        if mpz_sgn(tmp) != 0:
            non_square = 1
        else:
            mpz_sqrtrem(mpq_denref(z.value), tmp, mpq_denref(self.value))
            if mpz_sgn(tmp) != 0:
                non_square = 1
        mpz_clear(tmp)
        _sig_off

        if non_square:
            if not extend:
                raise ValueError, "square root of %s not a rational number"%self
            if prec:
                from sage.rings.real_mpfr import RealField
                K = RealField(prec)
                return K(self).sqrt(all=all)
            from sage.calculus.calculus import sqrt
            return sqrt(self, all=all)
        if all:
            return [z, -z]
        return z

    def period(self):
        r"""
        Return the period of the repeating part of the decimal
        expansion of this rational number.

        ALGORITHM: When a rational number $n/d$ with $(n,d)==1$ is
        expanded, the period begins after $s$ terms and has length
        $t$, where $s$ and $t$ are the smallest numbers satisfying
        $10^s=10^(s+t) (mod d)$.  When $d$ is coprime to 10, this
        becomes a purely periodic decimal with $10^t=1 (mod d)$.
        (Lehmer 1941 and Mathworld).

        EXAMPLES:
            sage: (1/7).period()
            6
            sage: RR(1/7)
            0.142857142857143
            sage: (1/8).period()
            1
            sage: RR(1/8)
            0.125000000000000
            sage: RR(1/6)
            0.166666666666667
            sage: (1/6).period()
            1
            sage: x = 333/106
            sage: x.period()
            13
            sage: RealField(200)(x)
            3.1415094339622641509433962264150943396226415094339622641509
        """
        cdef unsigned int alpha, beta
        d = self.denominator()
        alpha = d.valuation(2)
        beta = d.valuation(5)
        P = d.parent()
        if alpha > 0 or beta > 0:
            d = d//(P(2)**alpha * P(5)**beta)
        from sage.rings.integer_mod import Mod
        a = Mod(P(10),d)
        return a.multiplicative_order()

    def nth_root(self, int n):
        r"""
        Computes the nth root of self, or raises a \exception{ValueError}
        if self is not a perfect nth power.

        INPUT:
            n -- integer (must fit in C int type)

        AUTHOR:
           -- David Harvey (2006-09-15)

        EXAMPLES:
          sage: (25/4).nth_root(2)
          5/2
          sage: (125/8).nth_root(3)
          5/2
          sage: (-125/8).nth_root(3)
          -5/2
          sage: (25/4).nth_root(-2)
          2/5

          sage: (9/2).nth_root(2)
          Traceback (most recent call last):
          ...
          ValueError: not a perfect nth power

          sage: (-25/4).nth_root(2)
          Traceback (most recent call last):
          ...
          ValueError: cannot take even root of negative number

        """
        # TODO -- this could be quicker, by using GMP directly.
        cdef integer.Integer num
        cdef integer.Integer den
        cdef int negative

        if n > 0:
            negative = 0
        elif n < 0:
            n = -n
            negative = 1
        else:
            raise ValueError, "n cannot be zero"

        num, exact = self.numerator().nth_root(n, 1)
        if not exact:
            raise ValueError, "not a perfect nth power"

        den, exact = self.denominator().nth_root(n, 1)
        if not exact:
            raise ValueError, "not a perfect nth power"

        if negative:
            return den / num
        else:
            return num / den


    def str(self, int base=10):
        """
        INPUT:
            base -- integer (default: 10); base must be between 2 and 36.

        OUTPUT:
            string

        EXAMPLES:
            sage: (-4/17).str()
            '-4/17'
            sage: (-4/17).str(2)
            '-100/10001'

        Note that the base must be at most 36.
            sage: (-4/17).str(40)
            Traceback (most recent call last):
            ...
            ValueError: base (=40) must be between 2 and 36
            sage: (-4/17).str(1)
            Traceback (most recent call last):
            ...
            ValueError: base (=1) must be between 2 and 36
        """
        if base < 2 or base > 36:
            raise ValueError, "base (=%s) must be between 2 and 36"%base
        cdef size_t n
        cdef char *s

        n = mpz_sizeinbase (mpq_numref(self.value), base) \
            + mpz_sizeinbase (mpq_denref(self.value), base) + 3
        s = <char *>PyMem_Malloc(n)
        if s == NULL:
            raise MemoryError, "Unable to allocate enough memory for the string representation of an integer."

        _sig_on
        mpq_get_str(s, base, self.value)
        _sig_off
        k = <object> PyString_FromString(s)
        PyMem_Free(s)
        return k

    def __float__(self):
        """
        Return floating point approximation to self as a Python float.

        OUTPUT:
            float

        EXAMPLES:
            sage: (-4/17).__float__()
            -0.23529411764705882
            sage: float(-4/17)
            -0.23529411764705882
        """
        return mpq_get_d(self.value)

    def __hash__(self):
        """
        Return hash of self.

        OUTPUT:
            integer

        EXAMPLES:
            sage: (-4/17).__hash__()
            -19
            sage: (-5/1).__hash__()
            -5
        """
        cdef long n, d
        n = mpz_pythonhash(mpq_numref(self.value))
        d = mpz_pythonhash(mpq_denref(self.value))
        if d == 1:
            return n
        n = n ^ d
        if n == -1:
            return -2
        return n

    def __getitem__(self, int n):
        """
        Return n-th element of self, viewed as a list.  This is for
        consistency with how number field elements work.

        INPUT:
            n -- an integer  (error if not 0 or -1)

        OUTPUT:
            Rational

        EXAMPLES:
            sage: (-4/17)[0]
            -4/17
            sage: (-4/17)[1]
            Traceback (most recent call last):
            ...
            IndexError: index n (=1) out of range; it must be 0
            sage: (-4/17)[-1]   # indexing from the right
            -4/17
        """
        if n == 0 or n == -1:
            return self
        raise IndexError, "index n (=%s) out of range; it must be 0"%n

    ################################################################
    # Optimized arithmetic
    ################################################################
    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_add(x.value, self.value, (<Rational>right).value)
        return x

    cdef ModuleElement _iadd_c_impl(self, ModuleElement right):
        mpq_add(self.value, self.value, (<Rational>right).value)
        return self

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        # self and right are guaranteed to be Integers
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_sub(x.value, self.value, (<Rational>right).value)
        return x

    cdef ModuleElement _isub_c_impl(self, ModuleElement right):
        mpq_sub(self.value, self.value, (<Rational>right).value)
        return self

    cdef ModuleElement _neg_c_impl(self):
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_neg(x.value, self.value)
        return x

    cdef RingElement _mul_c_impl(self, RingElement right):
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        if mpz_sizeinbase (mpq_numref(self.value), 2)  > 100000 or \
             mpz_sizeinbase (mpq_denref(self.value), 2) > 100000:
            # We only use the signal handler (to enable ctrl-c out) in case
            # self is huge, so the product might actually take a while to compute.
            _sig_on
            mpq_mul(x.value, self.value, (<Rational>right).value)
            _sig_off
        else:
            mpq_mul(x.value, self.value, (<Rational>right).value)
        return x

    cdef RingElement _imul_c_impl(self, RingElement right):
        if mpz_sizeinbase (mpq_numref(self.value), 2)  > 100000 or \
             mpz_sizeinbase (mpq_denref(self.value), 2) > 100000:
            # We only use the signal handler (to enable ctrl-c out) in case
            # self is huge, so the product might actually take a while to compute.
            _sig_on
            mpq_mul(self.value, self.value, (<Rational>right).value)
            _sig_off
        else:
            mpq_mul(self.value, self.value, (<Rational>right).value)
        return self

    cdef RingElement _div_c_impl(self, RingElement right):
        """
        EXAMPLES:
            sage: 2/3
            2/3
            sage: 3/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Rational division by zero
        """
        if mpq_cmp_si((<Rational> right).value, 0, 1) == 0:
            raise ZeroDivisionError, "Rational division by zero"
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_div(x.value, self.value, (<Rational>right).value)
        return x

    cdef RingElement _idiv_c_impl(self, RingElement right):
        if mpq_cmp_si((<Rational> right).value, 0, 1) == 0:
            raise ZeroDivisionError, "Rational division by zero"
        mpq_div(self.value, self.value, (<Rational>right).value)
        return self


    ################################################################
    # Other arithmetic operations.
    ################################################################

    def __invert__(self):
        """
        Return the multiplicative inverse of this rational number.

        OUTPUT:
            Rational

        EXAMPLES:
            sage: (-4/17).__invert__()
            -17/4
            sage: ~(-4/17)
            -17/4
        """
        if self.is_zero():
            raise ZeroDivisionError, "rational division by zero"
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_inv(x.value, self.value)
        return x

    def __pow__(self, n, dummy):
        """
        Raise self to the integer power n.

        EXAMPLES:
            sage: (2/3)^5
            32/243
            sage: (-1/1)^(1/3)
            (-1)^(1/3)

        We raise to some interesting powers:
            sage: (2/3)^I
            2^I/3^I
            sage: (2/3)^sqrt(2)
            2^sqrt(2)/3^sqrt(2)
            sage: x,y,z,n = var('x,y,z,n')
            sage: (2/3)^(x^n + y^n + z^n)
            3^(-z^n - y^n - x^n)*2^(z^n + y^n + x^n)
            sage: (-7/11)^(tan(x)+exp(x))
            11^(-tan(x) - e^x)*(-7)^(tan(x) + e^x)
            sage: (2/3)^(3/4)
            2^(3/4)/3^(3/4)
            sage: (-1/3)^0
            1
            sage: a = (0/1)^(0/1); a
            1
            sage: type(a)
            <type 'sage.rings.rational.Rational'>

        The exponent must fit in a long unless the base is -1, 0, or 1.
            sage: s = (1/2)^(2^100)
            Traceback (most recent call last):
            ...
            RuntimeError: exponent must be at most 2147483647  # 32-bit
            RuntimeError: exponent must be at most 9223372036854775807 # 64-bit
            sage: s = (1/2)^(-2^100)
            Traceback (most recent call last):
            ...
            RuntimeError: exponent must be at most 2147483647  # 32-bit
            RuntimeError: exponent must be at most 9223372036854775807 # 64-bit
            sage: (-3/3)^(2^100)
            1

        This works even if the base is a float or Python complex or other type:
            sage: float(1.2)**(1/2)
            1.0954451150103321
            sage: complex(1,2)**(1/2)
            (1.272019649514069+0.786151377757423...j)
            sage: int(2)^(1/2)
            sqrt(2)
            sage: a = int(2)^(3/1); a
            8
            sage: type(a)
            <type 'sage.rings.rational.Rational'>
        """
        if dummy is not None:
            raise ValueError, "__pow__ dummy variable not used"

        if not PY_TYPE_CHECK(self, Rational):
            # If the base is not a rational, e.g., it is an int, complex, float, user-defined type, etc.
            try:
                self_coerced = Rational(self)
            except TypeError:
                n_coerced = type(self)(n)
                if n != n_coerced:
                    # dangerous coercion -- don't use -- try symbolic result
                    from sage.calculus.calculus import SR
                    return SR(self)**SR(n)
                return self.__pow__(n_coerced)
            return self_coerced.__pow__(n)

        cdef Rational _self = <Rational>self
        cdef long nn

        try:
            nn = PyNumber_Index(n)
        except TypeError:
            if PY_TYPE_CHECK(n, Rational):
                # this is the only sensible answer that avoids rounding and
                # an infinite recursion.
                from sage.calculus.calculus import SR
                return SR(self)**SR(n)
            if PY_TYPE_CHECK(n, Element):
                return (<Element>n)._parent(self)**n
            try:
                return n.parent()(self)**n
            except AttributeError:
                try:
                    return type(n)(self)**n
                except:
                    raise TypeError, "exponent (=%s) must be an integer.\nCoerce your numbers to real or complex numbers first."%n

        except OverflowError:
            if mpz_cmp_si(mpq_denref(_self.value), 1) == 0:
                if mpz_cmp_si(mpq_numref(_self.value), 1) == 0:
                    return self
                elif mpz_cmp_si(mpq_numref(_self.value), 0) == 0:
                    return self
                elif mpz_cmp_si(mpq_numref(_self.value), -1) == 0:
                    return self if n % 2 else -self
            raise RuntimeError, "exponent must be at most %s" % sys.maxint

        cdef Rational x = <Rational> PY_NEW(Rational)

        if nn == 0:
            mpq_set_si(x.value, 1, 1)
            return x

        if nn < 0:
            _sig_on
            # mpz_pow_ui(mpq_denref(x.value), mpq_numref(_self.value), <unsigned long int>(-nn))
            # mpz_pow_ui(mpq_numref(x.value), mpq_denref(_self.value), <unsigned long int>(-nn))
            # The above causes segfaults, so swap after instead...
            mpz_pow_ui(mpq_numref(x.value), mpq_numref(_self.value), -nn)
            mpz_pow_ui(mpq_denref(x.value), mpq_denref(_self.value), -nn)
            # mpz_swap(mpq_numref(x.value), mpq_denref(x.value)) # still a segfault
            mpq_inv(x.value, x.value)
            _sig_off
            return x
        elif nn > 0:
            _sig_on
            mpz_pow_ui(mpq_numref(x.value), mpq_numref(_self.value), nn)
            mpz_pow_ui(mpq_denref(x.value), mpq_denref(_self.value), nn)
            _sig_off
            return x


        if n < 0:  # this doesn't make sense unless n is an integer.
            x = _self**(-n)
            return x.__invert__()

        cdef mpz_t num, den

        _sig_on
        mpz_init(num)
        mpz_init(den)
        mpz_pow_ui(num, mpq_numref(_self.value), nn)
        mpz_pow_ui(den, mpq_denref(_self.value), nn)
        mpq_set_num(x.value, num)
        mpq_set_den(x.value, den)
        mpz_clear(num)
        mpz_clear(den)
        _sig_off

        return x














        _sig_on
        if nn < 0:  # we used to call (self**(-n)).__invert__()) -- this should be loads faster
            mpz_pow_ui(mpq_denref(x.value), mpq_numref(_self.value), -nn) #we switch den and num to invert
            mpz_pow_ui(mpq_numref(x.value), mpq_denref(_self.value), -nn)
        else:
            mpz_pow_ui(mpq_numref(x.value), mpq_numref(_self.value), nn)
            mpz_pow_ui(mpq_denref(x.value), mpq_denref(_self.value), nn)
        _sig_off

        #print "returning"

        return x

    def __pos__(self):
        """
        Return this rational number.

        OUTPUT:
            Rational

        EXAMPLES:
            sage: (-4/17).__pos__()
            -4/17
            sage: +(-4/17)
            -4/17
        """
        return self

    def __neg__(self):
        """
        Return the negative of this rational number.

        OUTPUT:
            Rational

        EXAMPLES:
            sage: (-4/17).__neg__()
            4/17
            sage: - (-4/17)
            4/17
        """
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_neg(x.value, self.value)
        return x

    def __nonzero__(self):
        """
        Return True if this rational number is nonzero.

        OUTPUT:
            bool

        EXAMPLES:
            sage: (-4/17).__nonzero__()
            True
            sage: (0/5).__nonzero__()
            False
            sage: bool(-4/17)
            True
        """
        # A rational number is zero iff its numerator is zero.
        return mpq_sgn(self.value) != 0

    def __abs__(self):
        """
        Return the absolute value of this rational number.

        OUTPUT:
            Rational

        EXAMPLES:
            sage: (-4/17).__abs__()
            4/17
            sage: abs(-4/17)
            4/17
        """
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_abs(x.value, self.value)
        return x

    def mod_ui(Rational self, unsigned long int n):
        """
        Return the remainder upon division of self by the unsigned long integer n.

        INPUT:
            n -- an unsigned long integer
        OUTPUT:
            integer

        EXAMPLES:
            sage: (-4/17).mod_ui(3)
            1
            sage: (-4/17).mod_ui(17)
            Traceback (most recent call last):
            ...
            ArithmeticError: The inverse of 0 modulo 17 is not defined.
        """
        cdef unsigned int num, den, a

        # Documentation from GMP manual:
        # "For the ui variants the return value is the remainder, and
        # in fact returning the remainder is all the div_ui functions do."
        _sig_on
        num = mpz_fdiv_ui(mpq_numref(self.value), n)
        den = mpz_fdiv_ui(mpq_denref(self.value), n)
        _sig_off
        return int((num * ai.inverse_mod_int(den, n)) % n)

    def __mod__(Rational self, other):
        """
        Return the remainder of division of self by other, where other
        is coerced to an integer

        INPUT:
            other -- object that coerces to an integer.
        OUTPUT:
            integer

        EXAMPLES:
            sage: (-4/17).__mod__(3/1)
            1
        """
        other = integer.Integer(other)
        if not other:
            raise ZeroDivisionError, "Rational modulo by zero"
        n = self.numer() % other
        d = self.denom() % other
        _sig_on
        d = d.inverse_mod(other)
        _sig_off
        return (n*d)%other

    def norm(self):
        """
        Returns the norm from Q to Q of x (which is just x). This was
        added for compatibility with NumberFields.

        OUTPUT:
            Rational -- reference to self

        EXAMPLES:
            sage: (1/3).norm()
             1/3

        AUTHOR:
          -- Craig Citro
        """
        return self

    def trace(self):
        """
        Returns the trace from Q to Q of x (which is just x). This was
        added for compatibility with NumberFields.

        OUTPUT:
            Rational -- reference to self

        EXAMPLES:
            sage: (1/3).trace()
             1/3

        AUTHOR:
          -- Craig Citro
        """
        return self

    def charpoly(self, var):
        """
        Return the characteristic polynomial of this rational number.
        This will always be just var - self; this is really here
        so that code written for number fields won't crash when
        applied to rational numbers.

        INPUT:
             var -- a string
        OUTPUT:
             Polynomial

        EXAMPLES:
            sage: (1/3).charpoly('x')
             x - 1/3

        AUTHOR:
          -- Craig Citro
        """
        QQ = self.parent()
        return QQ[var]([-self,1])

    def minpoly(self, var):
        """
        Return the minimal polynomial of this rational number.
        This will always be just x - self; this is really here
        so that code written for number fields won't crash when
        applied to rational numbers.

        INPUT:
             var -- a string
        OUTPUT:
             Polynomial

        EXAMPLES:
            sage: (1/3).minpoly('x')
             x - 1/3

        AUTHOR:
          -- Craig Citro
        """
        QQ = self.parent()
        return QQ[var]([-self,1])

    cdef integer.Integer _integer_c(self):
        if not mpz_cmp_si(mpq_denref(self.value), 1) == 0:
            raise TypeError, "no coercion of this rational to integer"
        cdef integer.Integer n
        n = PY_NEW(integer.Integer)
        n.set_from_mpz(mpq_numref(self.value))
        return n

    def _integer_(self):
        """
        Return self coerced to an integer.  Of course this rational number
        have a denominator of 1.

        OUTPUT:
            Integer

        EXAMPLES:
            sage: (-4/17)._integer_()
            Traceback (most recent call last):
            ...
            TypeError: no coercion of this rational to integer
            sage: (-4/1)._integer_()
            -4
        """
        return self._integer_c()

    # TODO -- this should be deprecated
    def numer(self):
        """
        Return the numerator of this rational number.

        EXAMPLE:
            sage: x = -5/11
            sage: x.numer()
            -5
        """
        cdef integer.Integer n
        n = PY_NEW(integer.Integer)
        n.set_from_mpz(mpq_numref(self.value))
        return n

    def numerator(self):
        """
        Return the numerator of this rational number.

        EXAMPLE:
            sage: x = 5/11
            sage: x.numerator()
            5

            sage: x = 9/3
            sage: x.numerator()
            3
        """
        cdef integer.Integer n
        n = PY_NEW(integer.Integer)
        n.set_from_mpz(mpq_numref(self.value))
        return n

    def __int__(self):
        """
        Return coercion of self to Python int.

        This takes the floor of self if self has a denominator (which
        is consistent with Python's long(floats)).

        EXAMPLES:
            sage: int(7/3)
            2
            sage: int(-7/3)
            -3
        """
        return int(self.__long__())

    def __long__(self):
        """
        Return coercion of self to Python long.

        This takes the floor of self if self has a denominator (which
        is consistent with Python's long(floats)).

        EXAMPLES:
            sage: long(7/3)
            2L
            sage: long(-7/3)
            -3L
        """
        cdef mpz_t x
        if mpz_cmp_si(mpq_denref(self.value),1) != 0:
            mpz_init(x)
            mpz_fdiv_q(x, mpq_numref(self.value), mpq_denref(self.value))
            n = mpz_get_pylong(x)
            mpz_clear(x)
            return n
        else:
            return mpz_get_pylong(mpq_numref(self.value))

    def denom(self):
        """
        self.denom(): Return the denominator of this rational number.

        EXAMPLES:
            sage: x = 5/13
            sage: x.denom()
            13
            sage: x = -9/3
            sage: x.denom()
            1
        """
        cdef integer.Integer n
        n = PY_NEW(integer.Integer)
        n.set_from_mpz(mpq_denref(self.value))
        return n

    def denominator(self):
        """
        self.denominator(): Return the denominator of this rational number.

        EXAMPLES:
            sage: x = -5/11
            sage: x.denominator()
            11
            sage: x = 9/3
            sage: x.denominator()
            1
        """
        cdef integer.Integer n
        n = PY_NEW(integer.Integer)
        n.set_from_mpz(mpq_denref(self.value))
        return n

    def factor(self):
        """
        Return the factorization of this rational number.

        OUTPUT:
            Factorization

        EXAMPLES:
            sage: (-4/17).factor()
            -1 * 2^2 * 17^-1

        Trying to factor 0 gives an arithmetic error:
            sage: (0/1).factor()
            Traceback (most recent call last):
            ...
            ArithmeticError: Prime factorization of 0 not defined.
        """
        return sage.rings.rational_field.factor(self)

    def support(self):
        """
        Return a sorted list of the primes where this rational number
        has non-zero valuation.

        OUTPUT:
            The set of primes appearing in the factorization of
            this rational with nonzero exponent, as a sorted list.

        EXAMPLES:
            sage: (-4/17).support()
            [2, 17]

        Trying to find the support of 0 gives an arithmetic error:
            sage: (0/1).support()
            Traceback (most recent call last):
            ...
            ArithmeticError: Support of 0 not defined.
        """
        if self.is_zero():
            raise ArithmeticError, "Support of 0 not defined."
        return sage.rings.arith.prime_factors(self)

    def gamma(self, prec=None):
        """
        Return the gamma function evaluated at self. This value is exact for
        integers and half-integers, otherwise a numerical approximation is
        returned.

        EXAMPLES:
            sage: gamma(1/2)
            sqrt(pi)
            sage: gamma(7/2)
            15*sqrt(pi)/8
            sage: gamma(-3/2)
            4*sqrt(pi)/3
            sage: gamma(6/1)
            120
            sage: gamma(1/3)
            2.67893853470775

        This function accepts an optional precision argument:
            sage: (1/3).gamma(prec=100)
            2.6789385347077476336556929410
            sage: (1/2).gamma(prec=100)
            1.7724538509055160272981674833
        """
        if prec is None:
            if mpz_cmp_ui(mpq_denref(self.value), 1) == 0:
                return integer.Integer(self).gamma()
            elif mpz_cmp_ui(mpq_denref(self.value), 2) == 0:
                numer = self.numer()
                rat_part = Rational((numer-2).multifactorial(2)) >> ((numer-1)//2)
                from sage.functions.constants import pi
                from sage.calculus.calculus import sqrt
                return sqrt(pi) * rat_part
            else:
                prec = 53
        return self.n(prec).gamma()

    def floor(self):
        """
        Return the floor of this rational number as an integer.

        OUTPUT:
            Integer

        EXAMPLES:
            sage: n = 5/3; n.floor()
            1
            sage: n = -17/19; n.floor()
            -1
            sage: n = -7/2; n.floor()
            -4
            sage: n = 7/2; n.floor()
            3
            sage: n = 10/2; n.floor()
            5
        """
        cdef integer.Integer n
        n = integer.Integer()
        mpz_fdiv_q(n.value, mpq_numref(self.value), mpq_denref(self.value))
        return n

    def ceil(self):
        """
        Return the ceiling of this rational number.

        OUTPUT:
            Integer

        If this rational number is an integer, this returns this
        number, otherwise it returns the floor of this number +1.

        EXAMPLES:
            sage: n = 5/3; n.ceil()
            2
            sage: n = -17/19; n.ceil()
            0
            sage: n = -7/2; n.ceil()
            -3
            sage: n = 7/2; n.ceil()
            4
            sage: n = 10/2; n.ceil()
            5
        """
        cdef integer.Integer n
        n = integer.Integer()
        mpz_cdiv_q(n.value, mpq_numref(self.value), mpq_denref(self.value))
        return n

    def round(Rational self, mode = "toward"):
        """
        Returns the nearest integer to self.

        INPUT:
            self -- a rational number
            mode -- a rounding mode for half integers:
                    'toward' (default) rounds toward zero
                    'away' rounds away from zero
                    'up' rounds up
                    'down' rounds down
                    'even' rounds toward the even integer
                    'odd' rounds toward the odd integer
        OUTPUT:
            Integer

        EXAMPLES:
            sage: n = 4/3; n.round()
            1
            sage: n = -17/4; n.round()
            -4
            sage: n = -5/2; n.round()
            -2
            sage: n.round("away")
            -3
            sage: n.round("up")
            -2
            sage: n.round("down")
            -3
            sage: n.round("even")
            -2
            sage: n.round("odd")
            -3
        """
        if not (mode in ['toward', 'away', 'up', 'down', 'even', 'odd']):
            raise ValueError, "rounding mode must be one of 'toward', 'away', 'up', 'down', 'even', or 'odd'"
        if self.denominator() == 1:
            from sage.rings.integer import Integer
            return Integer(self)
        if self.denominator() == 2:
            # round down:
            if (mode == "down") or \
                   (mode == "toward" and self > 0) or \
                   (mode == "away" and self < 0) or \
                   (mode == "even" and self.numerator() % 4 == 1) or \
                   (mode == "odd" and self.numerator() % 4 == 3):
                return self.numerator() // self.denominator()
            else:
                return self.numerator() // self.denominator() + 1
        else:
            q, r = self.numerator().quo_rem(self.denominator())
            if r < self.denominator() / 2:
                return q
            else:
                return q+1

    def height(self):
        """
        The max absolute value of the numerator and denominator of self,
        as an Integer.

        OUTPUT:
            Integer

        EXAMPLES:
            sage: a = 2/3
            sage: a.height()
            3
            sage: a = 34/3
            sage: a.height()
            34
            sage: a = -97/4
            sage: a.height()
            97

        AUTHOR:
            -- Naqi Jaffery (2006-03-05): examples
        """
        x = abs(self.numer())
        if x > self.denom():
            return x
        return self.denom()

    def _lcm(self, Rational other):
        """
        Returns the least common multiple, in the rational numbers,
        of self and other.  This function returns either 0 or 1 (as
        a rational number).

        INPUT:
            other -- Rational
        OUTPUT:
            Rational -- 0 or 1

        EXAMPLES:
            sage: (2/3)._lcm(3/5)
            1
            sage: (0/1)._lcm(0/1)
            0
            sage: type((2/3)._lcm(3/5))
            <type 'sage.rings.rational.Rational'>
        """
        if mpz_cmp_si(mpq_numref(self.value), 0) == 0 and \
               mpz_cmp_si(mpq_numref(other.value), 0) == 0:
            return Rational(0)
        return Rational(1)

    def _gcd(self, Rational other):
        """
        Returns the least common multiple, in the rational numbers,
        of self and other.  This function returns either 0 or 1 (as
        a rational number).

        INPUT:
            other -- Rational
        OUTPUT:
            Rational -- 0 or 1

        EXAMPLES:
            sage: (2/3)._gcd(3/5)
            1
            sage: (0/1)._gcd(0/1)
            0
        """
        if mpz_cmp_si(mpq_numref(self.value), 0) == 0 and \
               mpz_cmp_si(mpq_numref(other.value), 0) == 0:
            return Rational(0)
        return Rational(1)


    def additive_order(self):
        """
        Return the additive order of self.

        OUTPUT:
            integer or infinity

        EXAMPLES:
            sage: QQ(0).additive_order()
            1
            sage: QQ(1).additive_order()
            +Infinity
        """
        import sage.rings.infinity
        if self.is_zero():
            return integer.Integer(1)
        else:
            return sage.rings.infinity.infinity


    def multiplicative_order(self):
        """
        Return the multiplicative order of self.

        OUTPUT:
            Integer of infinity

        EXAMPLES:
            sage: QQ(1).multiplicative_order()
            1
            sage: QQ('1/-1').multiplicative_order()
            2
            sage: QQ(0).multiplicative_order()
            +Infinity
            sage: QQ('2/3').multiplicative_order()
            +Infinity
            sage: QQ('1/2').multiplicative_order()
            +Infinity
        """
        import sage.rings.infinity
        if self.is_one():
            return integer.Integer(1)
        elif mpz_cmpabs(mpq_numref(self.value),mpq_denref(self.value))==0:
	    # if the numerator and the denominator are equal in absolute value,
	    # then the rational number is -1
            return integer.Integer(2)
        else:
            return sage.rings.infinity.infinity

    def is_one(self):
        r"""
        Determine if a rational number is one.

        OUTPUT:
            bool

        EXAMPLES:
            sage: QQ(1/2).is_one()
            False
            sage: QQ(4/4).is_one()
            True
        """
        # A rational number is equal to 1 iff its numerator and denominator are equal
        return mpz_cmp(mpq_numref(self.value),mpq_denref(self.value))==0
        r"""Test if a rational number is zero

        EXAMPLES:

        sage: QQ(1/2).is_zero()
        False
        sage: QQ(0/4).is_zero()
        True
        """

    def is_integral(self):
        r"""
        Determine if a rational number is integral (i.e is in $\Z$).

        OUTPUT:
            bool

        EXAMPLES:
            sage: QQ(1/2).is_integral()
            False
            sage: QQ(4/4).is_integral()
            True
        """
        return bool(self in ZZ)

    cdef _lshift(self, long int exp):
        r"""
        Return $self*2^exp$
        """
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        _sig_on
        if exp < 0:
            mpq_div_2exp(x.value,self.value,-exp)
        else:
            mpq_mul_2exp(x.value,self.value,exp)
        _sig_off
        return x

    def __lshift__(x,y):
        """
        Left shift operator x << y.

        INPUT:
            x, y -- integer or rational
        OUTPUT:
            Rational

        EXAMPLES:
            sage: (2/3).__lshift__(4/1)
            32/3
            sage: (2/3).__lshift__(4/7)
            Traceback (most recent call last):
            ...
            ValueError: denominator must be 1
            sage: (2).__lshift__(4/1)
            32
            sage: (2/3).__lshift__(4)
            32/3
            sage: (2/3) << (4/1)
            32/3
        """
        if PY_TYPE_CHECK(x, Rational):
            if isinstance(y, (int, long, integer.Integer)):
                return (<Rational>x)._lshift(y)
            if PY_TYPE_CHECK(y, Rational):
                if mpz_cmp_si(mpq_denref((<Rational>y).value), 1) != 0:
                    raise ValueError, "denominator must be 1"
                return (<Rational>x)._lshift(y)
        return bin_op(x, y, operator.lshift)

    cdef _rshift(self, long int exp):
        r"""
        Return $self/2^exp$
        """
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        _sig_on
        if exp < 0:
            mpq_mul_2exp(x.value,self.value,-exp)
        else:
            mpq_div_2exp(x.value,self.value,exp)
        _sig_off
        return x

    def __rshift__(x,y):
        """
        Right shift operator x << y.

        INPUT:
            x, y -- integer or rational
        OUTPUT:
            Rational

        EXAMPLES:
            sage: (2/3).__rshift__(4/1)
            1/24
            sage: (2/3).__rshift__(4/7)
            Traceback (most recent call last):
            ...
            ValueError: denominator must be 1
            sage: (2).__rshift__(4/1)
            0
            sage: (2/1).__rshift__(4)
            1/8
            sage: (2/1) >>(4/1)
            1/8
        """
        if PY_TYPE_CHECK(x, Rational):
            if isinstance(y, (int, long, integer.Integer)):
                return (<Rational>x)._rshift(y)
            if PY_TYPE_CHECK(y, Rational):
                if mpz_cmp_si(mpq_denref((<Rational>y).value), 1) != 0:
                    raise ValueError, "denominator must be 1"
                return (<Rational>x)._rshift(y)
        return bin_op(x, y, operator.rshift)

    ##################################################
    # Support for interfaces
    ##################################################

    def _pari_(self):
        """
        Return this rational coerced into the PARI C library.

        OUTPUT:
            pari gen

        EXAMPLES:
            sage: (-9/17)._pari_()
            -9/17
        """
        return self.numerator()._pari_()/self.denominator()._pari_()

    def _interface_init_(self):
        """
        Return representation of this rational suitable for coercing
        into most any computer algebra system.

        OUTPUT:
            string

        EXAMPLES:
            sage: (2/3)._interface_init_()
            '2/3'
            sage: kash(3/1).Type()              # optional
            elt-fld^rat
            sage: magma(3/1).Type()             # optional
            FldRatElt
        """
        return '%s/%s'%(self.numerator(), self.denominator())



def pyrex_rational_reconstruction(integer.Integer a, integer.Integer m):
    """
    Find the rational reconstruction of a mod m, if it exists.

    INPUT:
        a -- Integer
        m -- Integer
    OUTPUT:
        x -- rings.rational.Rational

    EXAMPLES:
        sage: Integers(100)(2/3)
        34
        sage: sage.rings.rational.pyrex_rational_reconstruction(34, 100)
        2/3
    """
    cdef Rational x
    x = <Rational> PY_NEW(Rational)
    mpq_rational_reconstruction(x.value, a.get_value()[0], m.get_value()[0])
    return x

def make_rational(s):
    """
    Make a rational number from s (a string in base 32)

    INPUT:
        s -- string in base 32
    OUTPUT:
        Rational

    EXAMPLES:
        sage: (-7/15).str(32)
        '-7/f'
        sage: sage.rings.rational.make_rational('-7/f')
        -7/15
    """
    r = Rational()
    r._reduce_set(s)
    return r

cdef class Z_to_Q(Morphism):

    def __init__(self):
        """
        Create morphism from integers to rationals.

        EXAMPLES:
            sage: sage.rings.rational.Z_to_Q()
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
        """
        import integer_ring
        import rational_field
        import sage.categories.homset
        Morphism.__init__(self, sage.categories.homset.Hom(integer_ring.ZZ, rational_field.QQ))

    cpdef Element _call_(self, x):
        cdef Rational rat
        rat = <Rational> PY_NEW(Rational)
        mpq_set_z(rat.value, (<integer.Integer>x).value)
        return rat

    def _repr_type(self):
        """
        Return string that describes the type of morphism.

        EXAMPLES:
            sage: sage.rings.rational.Z_to_Q()._repr_type()
            'Natural'
        """
        return "Natural"

cdef class int_to_Q(Morphism):
    def __init__(self):
        """
        EXAMPLES:
            sage: sage.rings.rational.int_to_Q()
            Native morphism:
              From: Set of Python objects of type 'int'
              To:   Rational Field
        """
        import rational_field
        import sage.categories.homset
        from sage.structure.parent import Set_PythonType
        Morphism.__init__(self, sage.categories.homset.Hom(Set_PythonType(int), rational_field.QQ))

    cpdef Element _call_(self, a):
        cdef Rational rat
        rat = <Rational> PY_NEW(Rational)
        mpq_set_si(rat.value, PyInt_AS_LONG(a), 1)
        return rat

    def _repr_type(self):
        """
        Return string that describes the type of morphism.

        EXAMPLES:
            sage: sage.rings.rational.int_to_Q()._repr_type()
            'Native'
        """
        return "Native"

