"""
Field of real quad double numbers. These are deprecated.

sage: RQDF(1)
doctest:...: DeprecationWarning: RQDF is deprecated; use RealField(212) instead.
1.000000000000000000000000000000000000000000000000000000000000000

Quad double numbers allow us to represent real numbers with 212 bits
(or 64 decimal digits) of accuracy.  Computation of special functions
is extremely fast in this field, and arithmetic is slightly faster
than arithmetic in the MPFR real field.
Bloh
EXAMPLES:
RQDF numbers can be mixed with RDF and MPFR reals, Integers,
Rationals, etc.:
    sage: RQDF( 123.2) + RR (1.0)
    124.200000000000
    sage: RQDF( 12.2) + RDF (0.56)
    12.76
    sage: RQDF( 12.2) + (9)
    21.19999999999999928945726423989981412887573242187500000000000000
    sage: RQDF( 12.2) + (9/3)
    15.19999999999999928945726423989981412887573242187500000000000000

Note that the result will always be coerced to the field
with the lowest precision:
    sage: RR = RealField (300)
    sage: RQDF(123.2) * RR (.543)
    66.89760000000000154329882207093760371208190917968750000000000000

Mixing of symbolic an quad double elements:
    sage: a = RQDF(2) / log(10); a
    2.0/log(10)
    sage: parent(a)
    Symbolic Ring

Note that the following numerical imprecision is caused by passing via maxima:
    sage: RQDF(a)
    0.86858896380650365530225783783321016458879401160733313222890756...
"""


#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2007 didier deshommes <dfdeshom@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include '../ext/cdefs.pxi'
include '../ext/interrupt.pxi'


from sage.libs.mpfr cimport *

from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.real_mpfr cimport RealNumber, RealField
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage.structure.parent_base cimport ParentWithBase
from sage.structure.coerce_maps import NamedConvertMap

import operator

cdef extern from "solaris_fix.h": pass

from sage.misc.sage_eval import sage_eval

# if we don't import anything, cython complains with
# "undeclared name not builtin: sage"
import sage.rings

from sage.rings.real_double import RealDoubleElement

_R = None
cpdef RR():
    global _R
    if _R is None:
        from real_mpfr import RealField
        _R = RealField(212)
    return _R

cdef qd *qd_from_mpz(mpz_t z):
    cdef double d[4]
    cdef int i
    cdef mpz_t cur, res
    # The most significant double
    d[0] = mpz_get_d(z)
    if isinf(d[0]):
        d[0] = INFINITY
        return qd_from_double(d[0])

    if not finite(d[0]):
        d[0] = NAN
        return qd_from_double(d[0])
    mpz_init(cur)
    mpz_init(res)
    mpz_set_d(cur, d[0])
    mpz_sub(res, z, cur)
    # now we repeatedly find the most significant part of the remainder
    for i from 1 <= i < 4:
        d[i] = mpz_get_d(res)
        mpz_set_d(cur, d[i])
        mpz_sub(res, res, cur)
    mpz_clear(cur)
    mpz_clear(res)
    return qd_from_qd(d[0], d[1], d[2], d[3])

cdef qd *qd_from_mpfr(mpfr_t rr):
    cdef double d[4]
    cdef int i
    cdef mpfr_t cur, res
    cdef int isnan
    isnan = 0
    # The most significant double
    # We use GMP_RNDN here, which guarantees an exact
    # conversion if prec(rr) <= 4*53+3=215.
    d[0] = mpfr_get_d(rr, GMP_RNDN)
    mpfr_init2(cur, 53)
    mpfr_init2(res, mpfr_get_prec(rr))
    mpfr_set_d(cur, d[0], GMP_RNDN)
    mpfr_sub(res, rr, cur, GMP_RNDN)
    # now we repeatedly find the most significant part of the remainder
    for i from 1 <= i < 4:
        d[i] = mpfr_get_d(res, GMP_RNDN)
        mpfr_set_d(cur, d[i], GMP_RNDN)
        mpfr_sub(res, res, cur, GMP_RNDN)
    mpfr_clear(cur)
    mpfr_clear(res)
    # check if result is nan
    if 0 == d[0]: return qd_from_qd(0.0, 0.0, 0.0, 0.0)
    for i from 0 <= i < 3:
        if d[i] == d[i+1]: isnan += 1

    if 3 == isnan:
        return qd_from_double(NAN)
    else:
        return qd_from_qd(d[0], d[1], d[2], d[3])

cdef class RealQuadDoubleField_class(Field):
    """
    An approximation to a real number using quad double precision
    floating point numbers. Answers derived from calculations with
    such approximations may differ from what they would be if those
    calculations were performed with true real numbers. This is due
    to the rounding errors inherent to finite precision calculations.
    """
    def __init__(self):
        self._populate_coercion_lists_(coerce_list=[RealField(213)], embedding=NamedConvertMap(self, RealField(212), '_mpfr_'))

    def __dealloc__(self):
        pass

    cpdef bint is_exact(self) except -2:
        return False

    cdef _an_element_c_impl(self):  # override this in SageX
        return self(1.23)

    def _latex_(self):
        return "\\mathbf{R}"

    def __repr__(self):
        """
        Print out this real quad double field.

        EXAMPLES:
            sage: RQDF
            Real Quad Double Field

        """
        return "Real Quad Double Field"

    def __cmp__(self, x):
        """
        EXAMPLES:
            sage: RQDF == 1
            False
            sage: RQDF == RealQuadDoubleField()
            True

        """
        if PY_TYPE_CHECK(x, RealQuadDoubleField_class):
            return 0
        return cmp(type(self), type(x))

    def _element_constructor_(self, x):
        """
        Create a real quad double using x.

        EXAMPLES:
            sage: RQDF (-1)
            -1.000000000000000000000000000000000000000000000000000000000000000

            sage: RQDF (-1/3)
            -0.333333333333333333333333333333333333333333333333333333333333333

            sage: RQDF('-0.33333333333333333333')
            -0.333333333333333333330000000000000000000000000000000000000000000

        """
        if PY_TYPE_CHECK(x, QuadDoubleElement):
            return x
        elif hasattr(x,'_real_rqdf_'):
            return x._real_rqdf_(self)
        return QuadDoubleElement(x)

    def _magma_init_(self, magma):
        r"""
        Return a string representation of self in the Magma language.

        EXAMPLES:
            sage: magma(RQDF) # optional
            Real field of precision 63
            sage: floor(RR(log(2**212, 10)))
            63
        """
        return "RealField(%s : Bits := true)" % self.prec()

    def prec(self):
        """
        Return the precision of this real quad double field (to be more
        similar to RealField).  Always returns 212.

        EXAMPLES:
            sage: RQDF.prec()
            212
        """
        return 212

    def gen(self, i=0):
        """
        Return the generator of the field.

        EXAMPLES:
            sage: RQDF.gen()
            1.000000000000000000000000000000000000000000000000000000000000000

        """
        if i == 0:
            return self(int(1))
        else:
            raise IndexError

    def ngens(self):
        """
        Return the number of generators of the field.

        EXAMPLES:
            sage: RQDF.ngens()
            1
        """
        return 1

    def gens(self):
        """
        Returns a list of the generators of this field

        EXAMPLES:
            sage: RQDF.gens()
            [1.000000000000000000000000000000000000000000000000000000000000000]

        """
        return [self.gen()]

    def construction(self):
        """
        Returns the functorial construction of self, namely, completion of
        the rational numbers with respect to the prime at $\infinity$.

        Also preserves other information that makes this field unique
        (i.e. the Real Quad Double Field).

        EXAMPLES:
            sage: c, S = RQDF.construction(); S
            Rational Field
            sage: RQDF == c(S)
            True
        """
        from sage.categories.pushout import CompletionFunctor
        return (CompletionFunctor(sage.rings.infinity.Infinity,
                                  212,
                                  {'type': 'RQDF'}),
               sage.rings.rational_field.QQ)

    cpdef _coerce_map_from_(self, R):
        """
        Canonical coercion of x to the real quad double field.

        The rings that canonically coerce to the real quad double field are:
             * the real quad double field itself
             * int, long, integer, and rational rings
             * the mpfr real field, if its precision is at least 212 bits
             * anything that canonically coerces to the mpfr real field.

        EXAMPLES:
            sage: RQDF._coerce_(1/2)
            0.500000000000000000000000000000000000000000000000000000000000000
            sage: RQDF._coerce_(26)
            26.00000000000000000000000000000000000000000000000000000000000000
            sage: RQDF._coerce_(26r)
            26.00000000000000000000000000000000000000000000000000000000000000

            sage: RR = RealField (300)
            sage: RQDF._coerce_(RR('0.245465643523545656345356677563'))
            0.245465643523545656345356677563000000000000000000000000000000000
        """
        if R in (int, long, ZZ, QQ):
            return True

        if isinstance(R, RealField):
            return R.prec() > 212

    def name(self):
        return "QuadDoubleField"

    def __hash__(self):
        return hash(self.name())

    def pi(self):
        """
        Returns pi

        EXAMPLES:
            sage: RQDF.pi()
            3.141592653589793238462643383279502884197169399375105820974944592
            sage: RQDF.pi().sqrt()/2
            0.886226925452758013649083741670572591398774728061193564106903895

        """
        cdef qd z
        return QuadDoubleElement((z._pi.x[0], z._pi.x[1], z._pi.x[2], z._pi.x[3]))

    def log2(self):
        """
        Returns log(2)

        EXAMPLES:
            sage: RQDF.log2()
            0.693147180559945309417232121458176568075500134360255254120680009
        """
        cdef qd z
        return QuadDoubleElement((z._log2.x[0], z._log2.x[1], z._log2.x[2], z._log2.x[3]))

    def e(self):
        """
        Returns the natural log constant E

        EXAMPLES:
            sage: RQDF.e()
            2.718281828459045235360287471352662497757247093699959574966967628
        """
        cdef qd z
        return QuadDoubleElement((z._e.x[0], z._e.x[1], z._e.x[2], z._e.x[3]))

    def NaN(self):
        """
        Returns NaN

        EXAMPLES:
            sage: RQDF.NaN()
            'NaN'
        """
        return "NaN"

    def random_element(self,x=0,y=1):
        """
        Generate a random quad double between x and y
        RQDF.random_element() -- random real number between 0 and 1
        RQDF.random_element(n) -- return an real number between 0 and n-1, inclusive.
        RQDF.random_element(min, max) -- return a real number between min and max-1, inclusive.

        EXAMPLES:
            sage: RQDF.random_element(-10,10)
            2.070062766573917070526967401926937175104037963758214054970605794

            sage: RQDF.random_element(10)
            1.876986321449735379165332456424771674656237391205272916375069342

            sage: [RQDF.random_element(10) for _ in range(5)]
            [8.084281420673943787229609237458327203934426138730686265930328254,
            8.806298363947585826181224013949012875074955779827726577990392512,
            9.475217302151261462676501400096032242004946644547518035714857940,
            2.629993393648438854368552598518628332583301315453134053887857396,
            1.946435375426485270381085392826467583700841913266187372985425003]

        """
        # Switch to generating random numbers through MPFR (to make
        # the numbers work with Sage's global random number seeds).
        # Unfortunately, this takes about 3 times as long
        # as the c_qd_rand() code that's commented out below, but
        # I think it's worth it.  (15 microseconds vs. 5 microseconds,
        # on my laptop.)
        return RQDF(RR().random_element(y, x))
#         cdef QuadDoubleElement res, upper, lower
#         res = QuadDoubleElement(0)
#         upper = self(x)
#         lower = self(y)
#         c_qd_rand(res.initptr.x)
#         return (upper-lower)*res + lower

cdef class QuadDoubleElement(FieldElement):
    """
    A floating point approximation to a real number using quad
   double precision. Answers derived from calculations with such
   approximations may differ from what they would be if those
   calculations were performed with true real numbers. This is due
   to the rounding errors inherent to finite precision calculations.
    """
    cdef _new(self):
        cdef QuadDoubleElement q
        q = PY_NEW(QuadDoubleElement)
        q.initptr = new_qd_real()
        return q

    cdef _new_c(self, qd a):
        cdef QuadDoubleElement q
        q = PY_NEW(QuadDoubleElement)
        q.initptr =  qd_from_qd(a.x[0],a.x[1],a.x[2],a.x[3])
        return q

    def __dealloc__(self):
        cdef unsigned int cw
        fpu_fix_start(&cw)
        delete(self.initptr)
        fpu_fix_end(&cw)
        pass

    def __new__(self, x=None):
        # explicit cast required for C++
        self._parent = <ParentWithBase> _RQDF
        from sage.misc.misc import deprecation
        deprecation('RQDF is deprecated; use RealField(212) instead.')

    def __init__(self, x):
        """
        Create a quad double real number. A quad double real number is
        an unevaluated sum of 4 IEEE double numbers, capable of representing
        at least 212 bits of significand.The quad double number
        $(a_0,a_1,a_2,a_3)$ represents the exact sum $a=a_0+a_1+a_2+a_3$,
        where $a_0$ is the most significant component.

        EXAMPLES:
            sage: RQDF('1.1111111111100001111111011')
            1.111111111110000111111101100000000000000000000000000000000000000
            sage: RQDF(124/34)
            3.647058823529411764705882352941176470588235294117647058823529412
            sage: RQDF(12434)
            12434.00000000000000000000000000000000000000000000000000000000000
            sage: RQDF(2^60 + 9 )
            1.15292150460684698500000000000000000000000000000000000000000000e18

            You can also create a quad double number from a 4-tuple:
            sage: w= (2432323.0r,2.323e-12r,4.3423e-34r,-2.323e-52r)
            sage: RQDF(w)
            2432323.000000000002323000000000000098074113547887953405079397297

            or from RDF and MPFR reals:
            sage: RQDF(RDF(3434.34342))
            3434.343420000000151048880070447921752929687500000000000000000000
            sage: RQDF(RR(1091.34342))
            1091.343419999999923675204627215862274169921875000000000000000000
            sage: RR = RealField (300)
            sage: RQDF(RR('1091.34342'))
            1091.343420000000000000000000000000000000000000000000000000000000
        """
        self._set(x)

    cdef _set(self, x):
        cdef unsigned int cw
        fpu_fix_start(&cw)
        cdef qd *n,*d

        if PY_TYPE_CHECK(x, RealNumber):
            self.initptr = qd_from_mpfr((<RealNumber>x).value)

        elif PY_TYPE_CHECK(x, int):
            self.initptr = qd_from_int(x)

        elif PY_TYPE_CHECK(x, Integer):
            self.initptr = qd_from_mpz((<Integer>x).value)

        elif PY_TYPE_CHECK(x, Rational):
            n = qd_from_mpz(mpq_numref((<Rational>x).value))
            d = qd_from_mpz(mpq_denref((<Rational>x).value))
            self.initptr = new_qd_real()
            c_qd_div(n.x,d.x,self.initptr.x)

        elif PY_TYPE_CHECK(x, RealDoubleElement) or PY_TYPE_CHECK(x, float):
            self.initptr = qd_from_double(x)

        elif PY_TYPE_CHECK(x, QuadDoubleElement):
            n = (<QuadDoubleElement>x).initptr
            self.initptr = qd_from_qd(n.x[0],n.x[1],n.x[2],n.x[3])

        elif PY_TYPE_CHECK(x, long) or PY_TYPE_CHECK(x, str):
            s = str(x)
            _sig_on
            self.initptr = qd_from_str(s)
            _sig_off

        elif PY_TYPE_CHECK(x, tuple):
            _sig_on
            self.initptr = qd_from_qd(float(x[0]),float(x[1]),
                                      float(x[2]),float(x[3]))
            _sig_off

        else:
            self._set(RR()(x))

        fpu_fix_end(&cw)

    def get_doubles(self):
        """
        Return the 4 doubles that constitute this quad double real number

        EXAMPLES:
            sage: w= RQDF('1.34435343435344446376457677898097745635222')
            sage: t=w.get_doubles ()
            sage: RQDF (t)
            1.344353434353444463764576778980977456352220000000000000000000000
            sage: RQDF (t) == w
            True
        """
        return (self.initptr.x[0], self.initptr.x[1],
                self.initptr.x[2], self.initptr.x[3])

    def prec(self):
        """
        Return the precision of this real quad double element.  Always returns 212.

        EXAMPLES:
            sage: RQDF(0).prec()
            212
        """
        return 212

    def real(self):
        """
        Returns itself

        EXAMPLES:
            sage: w=RQDF(2) ; w.real()
            2.000000000000000000000000000000000000000000000000000000000000000
        """
        return self

    def imag(self):
        """
        Returns the imaginary part of this number (ie 0).

        EXAMPLES:
            sage: w=RQDF(2)
            sage: w.imag() == 0
            True
        """
        return QuadDoubleElement(0)

    def __complex__(self):
        """
        Returns self as a complex number
        EXAMPLES:
           sage: w=RQDF(2)
           sage: complex(w)
           (2+0j)

        """
        return complex(float(self))

    def __reduce__(self):
        """
        EXAMPLES:
            sage: s = dumps(RQDF('7.123456789'))
            sage: loads(s)
            7.123456789000000000000000000000000000000000000000000000000000000
        """
        doubles = self.get_doubles()
        return QuadDoubleElement,(doubles,)



    def __str__(self):
        return self.str()

    def __repr__(self):
        return self.str()

    def __str_no_scientific(self):
        """
        Returns a string representation of this number
        """
        cdef int MAX_DIGITS
        cdef char *s
        cdef int point_index
        result = ''
        MAX_DIGITS = 64

        # A negative number
        if self.initptr.x[0] < 0.0:
            result += '-'

        s = <char*>PyMem_Malloc(MAX_DIGITS+1)
        s[MAX_DIGITS]=0
        _sig_on
        self.initptr.to_digits(s,point_index,MAX_DIGITS)
        _sig_off

        t = str(s)
        PyMem_Free(s)

        # If this number is < 0
        if point_index < 0:
            point_index = - point_index
            result += '0.'

            # The value of point_index gives how many
            # additional zeroes there might be after the dot
            for i in range(1,point_index): result += '0'
            result +=  t[:63]
            return result

        # we always want to print the number with exactly 63 digits after
        # the dot
        result += t[:point_index+1] + '.' + t[point_index+1:63+point_index+1]

        return result

    def __str_scientific(self,precision=63):
        """
        Returns this number in scientific notation.
        """
        cdef char *s
        s = <char*>PyMem_Malloc(precision+8) # See docs for write()
        _sig_on
        self.initptr.write(s,precision,0,0)
        _sig_off
        t = str(s)
        PyMem_Free(s)
        return t

    def str(self):
        """
        Returns the string representation of self.

        If this number is less than $10^12$ and greater than $10^-12$,
        it is printed normally.

        Otherwise, it is printed in scientific notation

        EXAMPLES:
            sage: r= RQDF('222222222224.00000000001111111111') ; r
            222222222224.0000000000111111111100000000000000000000000000000000
            sage: r= RQDF('2222222222245.00000000001111111111') ; r
            2.22222222224500000000001111111111000000000000000000000000000000e12
            sage: r= RQDF('0.00000000001111111111') ; r
            0.0000000000111111111100000000000000000000000000000000000000000000000000000
            sage: r= RQDF('-0.00000000001111111111') ; r
            -0.0000000000111111111100000000000000000000000000000000000000000000000000000
            sage: RQDF(-10)/RQDF(0)
            -inf
            sage: RQDF(10.1)/RQDF(0)
            inf
            sage: RQDF(0)/RQDF(0)
            NaN
        """
        # 10**12 is used as the cutoff  because that is how RDF handles
        # large numbers
        cdef double x = self.initptr.x[0]
        if self.is_infinity():
            if x < 0: return "-inf"
            return "inf"
        if self.is_NaN(): return "NaN"
        cdef unsigned int cw
        fpu_fix_start(&cw)
        if 0.0 == x:
            result = self.__str_no_scientific()
        elif 1e-12 < x and x < 1e12:
            result = self.__str_no_scientific()
        elif -1e-12 > x and x > -1e12:
            result =  self.__str_no_scientific()
        else:
            result = self.__str_scientific()
        fpu_fix_end(&cw)
        return result

    def parent(self):
        """
        Returns the parent of this number

        EXAMPLES:
           sage: w=RQDF(-21.2) ; w.parent()
           Real Quad Double Field
        """
        return self._parent

    def __copy__(self):
        """
        Return copy of self, which since self is immutable, is just self.

        EXAMPLES:
            sage: w=RQDF('-21.2') ; copy(w)
            -21.20000000000000000000000000000000000000000000000000000000000000
        """
        return self

    def integer_part(self):
        """
        If in decimal this number is written n.defg, returns n.

        EXAMPLES:
            sage: test = 253536646425647436353675786864535364746
            sage: RQDF(test).integer_part() == test
            True

        """
        if 0 == self: return Integer("0")
        cdef unsigned int cw
        fpu_fix_start(&cw)
        s = self.__str_no_scientific()
        num = s.split('.')[0]
        result = Integer(num)
        fpu_fix_end(&cw)
        return result

    ########################
    #   Basic Arithmetic
    ########################
    def __invert__(self):
        """
        Compute the multiplicative inverse of self.

        EXAMPLES:
            sage: w=RQDF(1/3)
            sage: 1/w
            3.000000000000000000000000000000000000000000000000000000000000000
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_npwr(self.initptr.x,-1,res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Add two quad double numbers

        EXAMPLES:
            sage: RQDF(1/3) + RQDF(1)
            1.333333333333333333333333333333333333333333333333333333333333333
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        c_qd_add(self.initptr.x,(<QuadDoubleElement>right).initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        return res

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Substract two quad double numbers

        EXAMPLES:
            sage: RQDF(1/3) - RQDF(1)
            -0.666666666666666666666666666666666666666666666666666666666666666
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        c_qd_sub(self.initptr.x,(<QuadDoubleElement>right).initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        return res

    cpdef RingElement _mul_(self, RingElement right):
        """
        Multiply two quad double numbers

        EXAMPLES:
            sage: RQDF('1.3') * RQDF(10)
            13.00000000000000000000000000000000000000000000000000000000000000
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        c_qd_mul(self.initptr.x,(<QuadDoubleElement>right).initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        return res

    cpdef RingElement _div_(self, RingElement right):
        """
        Divide two quad double numbers

        EXAMPLES:
            sage: RQDF(1/3) / RQDF(100)
            0.00333333333333333333333333333333333333333333333333333333333333333
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        c_qd_div(self.initptr.x,(<QuadDoubleElement>right).initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        return res

    cpdef ModuleElement _neg_(self):
        """
        Negates a quad double number.

        EXAMPLES:
            sage: -RQDF('0.056')
            -0.0560000000000000000000000000000000000000000000000000000000000000
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)          # This might not be needed here.
        res = self._new()
        c_qd_neg(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        return res

    def __abs__(self):
        """
        Returns the absolute value of a quad double number.

        EXAMPLES:
            sage: abs(RQDF('-0.45'))
            0.450000000000000000000000000000000000000000000000000000000000000
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        c_qd_abs(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        return res

    def __lshift__(self, n):
        """
        LShifting a quad double is not supported
        """
        raise TypeError, "unsupported operand type(s) for <<: '%s' and '%s'"%(type(self), type(n))

    def __rshift__(self, n):
        """
        RShifting a quad double is not supported
        """
        raise TypeError, "unsupported operand type(s) for >>: '%s' and '%s'"%(type(self), type(n))

    def multiplicative_order(self):
        """
        Returns the multiplicative order of self

        EXAMPLES:
            sage: w=RQDF(-1)
            sage: w.multiplicative_order()
            -1
            sage: w=RQDF(1)
            sage: w.multiplicative_order()
            1
            sage: w=RQDF(0)
            sage: w.multiplicative_order()
            +Infinity

        """
        if self == 1: return 1
        if self == -1: return -1
        return sage.rings.infinity.infinity


    ###################
    # Rounding etc
    ###################

    def round(self):
        """
        Given real number x, rounds up if fractional part is greater than .5,
        rounds down if fractional part is lesser than .5.
        EXAMPLES:
            sage: RQDF(0.49).round()
            0.000000000000000000000000000000000000000000000000000000000000000
            sage: RQDF(0.51).round()
            1.000000000000000000000000000000000000000000000000000000000000000
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        if self.frac() < 0.5:
            _sig_on
            c_qd_floor(self.initptr.x,res.initptr.x)
            fpu_fix_end(&cw)
            _sig_off
            return res

        _sig_on
        c_qd_ceil(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res


    def floor(self):
        """
        Returns the floor of this number

        EXAMPLES:
            sage: RQDF(2.99).floor()
            2
            sage: RQDF(2.00).floor()
            2
            sage: RQDF(-5/2).floor()
            -3
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        c_qd_floor(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        return res.integer_part()

    def ceil(self):
        """
        Returns the ceiling of this number, as an integer.

        EXAMPLES:
            sage: RQDF(2.99).ceil()
            3
            sage: RQDF(2.00).ceil()
            2
            sage: RQDF(-5/2).ceil()
            -2
            sage: type(RQDF(-5/2).ceil())
            <type 'sage.rings.integer.Integer'>
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        c_qd_ceil(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        return res.integer_part()

    def ceiling(self):
        """
        EXAMPLES:
            sage: RQDF('2.000000000001').ceiling()
            3
        """
        return self.ceil()

    def trunc(self):
        """
        Truncates this number (returns integer part).

        EXAMPLES:
            sage: RQDF(2.99).trunc()
            2.000000000000000000000000000000000000000000000000000000000000000
            sage: RQDF(-2.00).trunc()
            -2.000000000000000000000000000000000000000000000000000000000000000
            sage: RQDF(0.00).trunc()
            0.000000000000000000000000000000000000000000000000000000000000000
        """
        return QuadDoubleElement(self.floor())


    def frac(self):
        """
        frac returns a real number > -1 and < 1. that satisfies the
        relation:
            x = x.trunc() + x.frac()

        EXAMPLES:
            sage: RQDF('2.99').frac()
            0.990000000000000000000000000000000000000000000000000000000000000
            sage: RQDF(2.50).frac()
            0.500000000000000000000000000000000000000000000000000000000000000
            sage: RQDF('-2.79').frac()
            -0.790000000000000000000000000000000000000000000000000000000000000
        """
        return self - self.integer_part()


    ###########################################
    # Conversions
    ###########################################

    def __float__(self):
        """
        Returns the floating-point value of this number

        EXAMPLES:
            sage: w = RQDF('-23.79'); w
            -23.79000000000000000000000000000000000000000000000000000000000000
            sage: float(w)
            -23.789999999999999
        """
        cdef double d
        cdef unsigned int cw
        fpu_fix_start(&cw)
        d = qd_to_double(qd_deref(self.initptr))
        fpu_fix_end(&cw)
        return d

    def _rpy_(self):
        """
        Returns self.__float__() for rpy to convert into the
        appropriate R object.

        EXAMPLES:
            sage: n = RQDF(2.0)
            sage: n._rpy_()
            2.0
            sage: type(n._rpy_())
            <type 'float'>
        """
        return self.__float__()

    def __int__(self):
        """
        Returns integer truncation of this real number.

        EXAMPLES:
            sage: w = RQDF('-23.79'); w
            -23.79000000000000000000000000000000000000000000000000000000000000
            sage: int(w)
            -23
        """
        cdef int i
        cdef unsigned int cw
        fpu_fix_start(&cw)
        i = qd_to_int(qd_deref(self.initptr))
        fpu_fix_end(&cw)
        return i

    def __long__(self):
        """
        Returns long integer truncation of this real number.

        EXAMPLES:
            sage: w = RQDF.pi(); w
            3.141592653589793238462643383279502884197169399375105820974944592
            sage: long(w)
            3L
        """
        return long(self.__int__())


    def _real_double_(self, R):
        """
        EXAMPLES:
            sage: w = RQDF.e(); w
            2.718281828459045235360287471352662497757247093699959574966967628
            sage: RDF(w)
            2.71828182846
            sage: w._real_double_(RDF)
            2.71828182846
        """
        return  R(float(self))

    def _mpfr_(self, RealField R):
        """
        TESTS:
            sage: w = RQDF('2.345001').sqrt(); w
            1.531339609622894852128128425884749978483262262653204338472911277
            sage: RealField(212)(w)
            1.53133960962289485212812842588474997848326226265320433847291128

            sage: RR(RQDF('324324.0098736633445565765349760000276353865'))
            324324.009873663
            sage: R200 = RealField (200)
            sage: R200(RQDF('324324.0098736633445565765349760000276353865'))
            324324.00987366334455657653497600002763538650000000000000000

            sage: w = RQDF('2.345001').sqrt(); w
            1.531339609622894852128128425884749978483262262653204338472911277
            sage: RealField(212)(w)
            1.53133960962289485212812842588474997848326226265320433847291128
            sage: RealField(250)(w)
            1.5313396096228948521281284258847499784832622626532043384729112767048720070
        """
        cdef int i
        cdef mpfr_t curr
        cdef mpfr_rnd_t rnd = (<RealField>R).rnd
        cdef RealNumber rr = R()
        mpfr_set_d(rr.value, self.initptr.x[0], rnd)
        mpfr_init2(curr, 53)
        for i from 1 <= i < 4:
            mpfr_set_d(curr, self.initptr.x[i], rnd)
            mpfr_add(rr.value, rr.value, curr, rnd)
        mpfr_clear(curr)
        return rr

    def _complex_double_(self, C):
        """
        EXAMPLES:
            sage: RQDF(-1)._complex_double_(CDF)
            -1.0
            sage: CDF(RQDF(-1))
            -1.0
        """
        return  C(float(self))

    def _complex_mpfr_field_(self, K):
        """
        EXAMPLES:
            sage: w = RQDF('2.345001').sqrt(); w
            1.531339609622894852128128425884749978483262262653204338472911277
            sage: w._complex_mpfr_field_(ComplexField(212))
            1.53133960962289485212812842588474997848326226265320433847291128
            sage: ComplexField(212)(w)
            1.53133960962289485212812842588474997848326226265320433847291128
        """
        return K(repr(self))

    def _pari_(self):
        """
        Return the PARI real number corresponding to self.

        EXAMPLES:
            sage: w = RQDF('2').sqrt(); w
            1.414213562373095048801688724209698078569671875376948073176679738
            sage: x = w._pari_()

        Note that x appears to have lower precision:
            sage: x
            1.41421356237310

        In fact it doesn't.
            sage: pari.set_real_precision(64)
            15
            sage: x
            1.414213562373095048801688724209698078569671875376948073176679738

        Set the precision in PARI back.
            sage: pari.set_real_precision(15)
            64
        """
        return  sage.libs.pari.all.pari.new_with_bits_prec(repr(self), 212)


    ###########################################
    # Comparisons: ==, !=, <, <=, >, >=
    ###########################################

    def is_NaN(self):
        """
        Returns True if self is NaN

        EXAMPLES:
            sage: w=RQDF(9)
            sage: w.is_NaN()
            False

            sage: w=RQDF(-10)/RQDF(0)
            sage: (0*w).is_NaN()
            True
            sage: (RQDF(0)/0).is_NaN()
            True
            sage: w=RQDF(-10)/RQDF(0)
            sage: (w/w).is_NaN()
            True
            sage: (w-w).is_NaN()
            True
        """
        return qd_is_nan(qd_deref(self.initptr))

    def is_infinity(self):
        """
        Returns True if self is infinifty

        EXAMPLES:
            sage: w=RQDF(32)/RQDF(0) ; w.is_infinity()
            True
            sage: w=RQDF(-10)/RQDF(0) ; w.is_infinity()
            True
            sage: w=RQDF(1<<2900000)
            sage: w.is_infinity ()
            True
        """
        return qd_is_inf(qd_deref(self.initptr))

    def is_positive_infinity(self):
        """
        Returns True if this number is positive infinity

        EXAMPLES:
            sage: r=RQDF(13)/RQDF(0)
            sage: r.is_positive_infinity ()
            True
        """
        cdef double x = self.initptr.x[0]
        if self.is_infinity() and x > 0 :
            return True
        return False

    def is_negative_infinity(self):
        """
        Returns True if this number is negative infinity

        EXAMPLES:
            sage: r=RQDF(-13)/RQDF(0)
            sage: r.is_negative_infinity ()
            True
            sage: r=RQDF(0)
            sage: r.is_negative_infinity ()
            False
        """
        cdef double x = self.initptr.x[0]
        if self.is_infinity() and x < 0 :
            return True
        return False

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Compares 2 quad double numbers

        Returns True if self is NaN
        EXAMPLES:
            sage: RQDF('3233') > RQDF('323')
            True
            sage: RQDF('-3233') > RQDF('323')
            False
            sage: RQDF('-3233') == RQDF('-3233')
            True
        """
        cdef int i
        c_qd_comp(left.initptr.x,(<QuadDoubleElement>right).initptr.x,&i)
        return i

    ############################
    # Special Functions
    ############################

    def sqrt(self, extend=True, all=False):
        """
        The square root function.

        INPUT:
            extend -- bool (default: True); if True, return
                a square root in a complex field if necessary
                if self is negative; otherwise raise a ValueError
            all -- bool (default: False); if True, return a list
                of all square roots.

        EXAMPLES:
            sage: RQDF('-3233').sqrt()
            56.8594759033179974315326998383777083570361710615905085284278958*I
            sage: RQDF('-4').sqrt()
            2.00000000000000000000000000000000000000000000000000000000000000*I
            sage: RQDF('4.1234').sqrt()
            2.030615670184784130018521630263322574805816770689719741871597292
            sage: w = RQDF(2).sqrt(); w
            1.414213562373095048801688724209698078569671875376948073176679738
            sage: RealField(212)(2).sqrt()
            1.41421356237309504880168872420969807856967187537694807317667974
            sage: w = RQDF(-2).sqrt(); w
            1.41421356237309504880168872420969807856967187537694807317667974*I
            sage: w = RQDF('2.345001').sqrt(); w
            1.531339609622894852128128425884749978483262262653204338472911277
            sage: w = RQDF(2).sqrt(all=True); w
            [1.414213562373095048801688724209698078569671875376948073176679738,
            -1.414213562373095048801688724209698078569671875376948073176679738]

            sage: RQDF(-2).sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: negative number -2.000000000000000000000000000000000000000000000000000000000000000 does not have a square root in the real quad double field
        """
        if self.initptr.x[0] ==0:
            if all:
                return [self]
            else:
                return self
        elif self.initptr.x[0] > 0:
            z = self._square_root()
            if all:
                return [z, -z]
            return z


        if not extend:
            raise ValueError, "negative number %s does not have a square root in the real quad double field"%self

        from sage.rings.complex_field import ComplexField
        return ComplexField(212)(self).sqrt(all=all)


    def _square_root(self):
        """
        Return a square root of self.  A real number will always be
        returned (though it will be NaN if self is negative).

        Use self.sqrt() to get a complex number if self is negative.

        EXAMPLES:
            sage: RQDF('4')._square_root()
            2.000000000000000000000000000000000000000000000000000000000000000
            sage: RQDF('-4')._square_root()
            NaN
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        if self.initptr.x[0] > 0:
            c_qd_sqrt(self.initptr.x, res.initptr.x)
            fpu_fix_end(&cw)
            return res
        else:
            res = self._new_c(self.initptr._nan)
            fpu_fix_end(&cw)
            return res

    def cube_root(self):
        """
        Return the cubic root (defined over the real numbers) of self.

        EXAMPLES:
            sage: r = RQDF(125.0); r.cube_root()
            5.000000000000000000000000000000000000000000000000000000000000000
            sage: RQDF('-8').cube_root()
            -2.000000000000000000000000000000000000000000000000000000000000000
        """
        return self.nth_root(3)


    def nth_root(self, int n):
        """
        Returns the $n^{th}$ root of self.
        Returns NaN if self is negative and $n$ is even.

        EXAMPLES:
            sage: r = RQDF(125.0); r.nth_root(3)
            5.000000000000000000000000000000000000000000000000000000000000000
            sage: r.nth_root(5)
            2.626527804403767236455131266496479582115662802810898530034436330
            sage: RQDF('-4987').nth_root(2)
            NaN
        """

        cdef QuadDoubleElement res
        cdef double neg[4]
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        if self.initptr.is_negative() and n % 2 == 1:
            c_qd_neg(self.initptr.x, neg)
            c_qd_nroot(neg, n, res.initptr.x)
            c_qd_neg(res.initptr.x, res.initptr.x)
        else:
            c_qd_nroot(self.initptr.x, n, res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    def __pow(self, n, modulus):
        cdef QuadDoubleElement res, n2
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_npwr(self.initptr.x, n, res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    def __pow__(self,n,d):
        """
        Compute self raised to the power of exponent.

        EXAMPLES:
            sage: a = RQDF('1.23456')
            sage: a^20
            67.64629770385396909318545781437142167431212492061816351749107831
            sage: RQDF (2)^3.2
            9.18958683997628
        """
        if not isinstance(n,(Integer, int)):
            res = n * self.log()
            res = res.exp()
        else:
            res = self.__pow(n,d)
        return res

    def log(self):
        """
        Returns the log of this number.

        Returns NaN is self is negative,a nd minus infinity if self is 0.

        EXAMPLES:
            sage: RQDF(2).log()
            0.693147180559945309417232121458176568075500134360255254120680009
            sage: RQDF(0).log()
            -inf
            sage: RQDF(-1).log()
            NaN
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        if self.initptr.x[0] < 0.0:
            res = self._new_c(self.initptr._nan)
        elif self.initptr.x[0] == 0.0:
            res = self._new()
            res.initptr = qd_from_double(-INFINITY)
        else:
            res = self._new()
            _sig_on
            c_qd_log(self.initptr.x,res.initptr.x)
            _sig_off
        fpu_fix_end(&cw)

        return res

    def log10(self):
        """
        Returns log to the base 10 of self

        Returns NaN is self is negative,a nd minus infinity if self is 0.

        EXAMPLES:
            sage: r = RQDF(16); r.log10()
            1.204119982655924780854955578897972107072759525848434165241709845
            sage: r.log() / RQDF(log(10))
            1.204119982655924780854955578897972107072759525848434165241709844
            sage: r = RQDF('-16.0'); r.log10()
            NaN

        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        if self.initptr.x[0] < 0.0:
            res = self._new_c(self.initptr._nan)
        elif self.initptr.x[0] == 0.0:
            res = self._new()
            res.initptr = qd_from_double(-INFINITY)
        else:
            res = self._new()
            _sig_on
            c_qd_log10(self.initptr.x,res.initptr.x)
            _sig_off
        fpu_fix_end(&cw)

        return res

    def exp(self):
        r"""
        Returns $e^\code{self}$

        EXAMPLES:
            sage: r = RQDF(0.0) ; r.exp()
            1.000000000000000000000000000000000000000000000000000000000000000
            sage: RQDF('-32.3').exp()
            9.38184458849865778521574138078207776624385053667970865778536257e-15
            sage: r = RQDF('16.0');
            sage: r.log().exp() == r
            True
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_exp(self.initptr.x,res.initptr.x)
        _sig_off
        fpu_fix_end(&cw)
        return res

    def cos(self):
        """
        Returns the cosine of this number

        EXAMPLES:
            sage: t=RQDF(pi/2)
            sage: t.cos()
            -2.84867032372793962424197268086858714626628975008598348988393400e-65
            sage: t.cos()^2 + t.sin()^2
            1.000000000000000000000000000000000000000000000000000000000000000
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_cos(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    def sin(self):
        """
        Returns the sine of this number

        EXAMPLES:
            sage: RQDF(pi).sin()
            -5.69734064745587924848394536173717429253257950017196697976786799e-65
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_sin(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    def tan(self):
        """
        Returns the tangent of this number

        EXAMPLES:
            sage: q = RQDF(pi/3)
            sage: q.tan()
            1.732050807568877293527446341505872366942805253810380628055806980
            sage: q = RQDF(pi/6)
            sage: q.tan()
            0.577350269189625764509148780501957455647601751270126876018602326
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_tan(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    def sincos(self):
        """
        Returns a pair consisting of the sine and cosine.

        EXAMPLES:
            sage: t = RQDF(pi/6)
            sage: t.sincos()
            (0.500000000000000000000000000000000000000000000000000000000000000, 0.866025403784438646763723170752936183471402626905190314027903489)
        """
        return self.sin(), self.cos()

    def arccos(self):
        """
        Returns the inverse cosine of this number

        EXAMPLES:
            sage: q = RQDF(pi/3)
            sage: i = q.cos()
            sage: q
            1.047197551196597746154214461093167628065723133125035273658314864
            sage: i.arccos()
            1.047197551196597746154214461093167628065723133125035273658314864
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_acos(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    def arcsin(self):
        """
        Returns the inverse sine of this number

        EXAMPLES:
            sage: q = RQDF(pi/3)
            sage: i = q.sin()
            sage: q
            1.047197551196597746154214461093167628065723133125035273658314864
            sage: i.arcsin()
            1.047197551196597746154214461093167628065723133125035273658314864
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_asin(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    def arctan(self):
        """
        Returns the inverse tangent of this number

        EXAMPLES:
            sage: q = RQDF(pi/3)
            sage: i = q.tan()
            sage: q
            1.047197551196597746154214461093167628065723133125035273658314864
            sage: i.arctan()
            1.047197551196597746154214461093167628065723133125035273658314864
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_atan(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    def cosh(self):
        """
        Returns the hyperbolic cosine of this number

        EXAMPLES:
            sage: q = RQDF(pi/12)
            sage: q.cosh()
            1.034465640095510565271865251179886560959831568117717546138668562
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_cosh(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    def sinh(self):
        """
        Returns the hyperbolic sine of this number

        EXAMPLES:
            sage: q = -RQDF(pi/12)
            sage: q.sinh()
            -0.264800227602270757698096543949405541727737186661923151601337995
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_sinh(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    def tanh(self):
        """
        Returns the hyperbolic tangent of this number

        EXAMPLES:
            sage: q = RQDF(pi/12)
            sage: q.tanh()
            0.255977789245684539459617840766661476446446454939881314446181622
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_tanh(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    def arccosh(self):
        """
        Returns the hyperbolic inverse cosine of this number

        EXAMPLES:
            sage: q = RQDF(pi/2)
            sage: i = q.cosh() ; i
            2.509178478658056782009995643269405948212024358148152274047975682
            sage: i.arccosh()
            1.570796326794896619231321691639751442098584699687552910487472296
            sage: q
            1.570796326794896619231321691639751442098584699687552910487472296
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_acosh(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    def arcsinh(self):
        """
        Returns the hyperbolic inverse sine of this number

        EXAMPLES:
            sage: q = RQDF(pi/2)
            sage: i = q.sinh() ; i
            2.301298902307294873463040023434427178178146516516382665972839798
            sage: i.arcsinh() ; q
            1.570796326794896619231321691639751442098584699687552910487472296
            1.570796326794896619231321691639751442098584699687552910487472296
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_asinh(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    def arctanh(self):
        """
        Returns the hyperbolic inverse tangent of this number

        EXAMPLES:
            sage: q = RQDF(pi/2)
            sage: i = q.tanh() ; i
            0.917152335667274346373092921442618775367927148601088945343574124
            sage: i.arctanh() ; q
            1.570796326794896619231321691639751442098584699687552910487472291
            1.570796326794896619231321691639751442098584699687552910487472296
        """
        cdef QuadDoubleElement res
        cdef unsigned int cw
        fpu_fix_start(&cw)
        res = self._new()
        _sig_on
        c_qd_atanh(self.initptr.x,res.initptr.x)
        fpu_fix_end(&cw)
        _sig_off
        return res

    def agm(self, other):
        """
        Return the arithmetic-geometric mean of self and other. The
        arithmetic-geometric mean is the common limit of the sequences
        $u_n$ and $v_n$, where $u_0$ is self, $v_0$ is other,
        $u_{n+1}$ is the arithmetic mean of $u_n$ and $v_n$, and
        $v_{n+1}$ is the geometric mean of u_n and v_n. If any operand
        is negative, the return value is \code{NaN}.

        EXAMPLES:
            sage: RQDF(2).sqrt().agm(RQDF(3).sqrt())
            1.569105802869322326985195456078256167313945200090173796316846190
        """
        R = RR()
        return QuadDoubleElement(R(self).agm(R(other)))

cdef RealQuadDoubleField_class _RQDF
_RQDF = RealQuadDoubleField_class()

RQDF = _RQDF   # external interface

def RealQuadDoubleField():
    global _RQDF
    return _RQDF

def is_QuadDoubleElement(x):
    return PY_TYPE_CHECK(x, QuadDoubleElement)
