include '../ext/cdefs.pxi'
include '../ext/interrupt.pxi'
include '../gsl/gsl_complex.pxi'

import operator

from sage.misc.sage_eval import sage_eval

cimport sage.ext.element
import sage.ext.element

cimport sage.ext.ring
import sage.ext.ring

import complex_number


cdef class ComplexDoubleField_class(sage.ext.ring.Field):
    """
    The field of complex double precision numbers.

    ALGORITHM: Arithmetic is done using GSL (the GNU Scientific Library).
    """

    def __cmp__(self, other):
        """
        Returns True if and only if other is the unique complex double field.

        EXAMPLES:
            sage: CC == CDF
            False
            sage: CDF == ComplexDoubleField     # CDF is the shorthand
            True
        """
        if other is ComplexDoubleField:
            return 0
        return -1

    def __repr__(self):
        """
        Print out this complex double field.

        EXAMPLES:
            sage: ComplexDoubleField
            Complex Double Field
            sage: CDF
            Complex Double Field
        """
        return "Complex Double Field"

    def __call__(self, x, im=None):
        """
        Create a complex double using x and optionally an imaginary part im.

        EXAMPLES:
            sage: CDF(0,1)
            1.0*I
            sage: CDF(2/3)
            0.666666666667
            sage: CDF(5)
            5.0
            sage: CDF('i')
            1.0*I
            sage: CDF(complex(2,-3))
            2.0 - 3.0*I
            sage: CDF(4.5)
            4.5
            sage: CDF(1+I)
            1.0 + 1.0*I

        A TypeError is raised if the coercion doesn't make sense:
            sage: CDF(QQ['x'].0)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to float

        One can convert back and forth between double precision complex
        numbers and higher-precision ones, though of course there may
        be loss of precision:
            sage: a = ComplexField(200)(-2).sqrt(); a
            1.4142135623730950488016887242096980785696718753769480731766796*I
            sage: b = CDF(a); b
            1.41421353817*I
            sage: a.parent()(b)
            1.4142135623700000000000000000000000000000000000000000000000002*I
        """
        if im is None:
            if isinstance(x, ComplexDoubleElement):
                return x
            elif isinstance(x, (float, int, long)):
                return ComplexDoubleElement(x, 0)
            elif isinstance(x, complex):
                return ComplexDoubleElement(x.real, x.imag)
            elif isinstance(x, complex_number.ComplexNumber):
                return ComplexDoubleElement(x.real(), x.imag())
            elif isinstance(x, str):
                return sage_eval(x.replace(' ',''), locals={"I":self.gen(),"i":self.gen()})
            else:
                return ComplexDoubleElement(x, 0)
        else:
            return ComplexDoubleElement(x, im)

    def gen(self, n=0):
        """
        Return the generator of the complex double field.
        EXAMPLES:
            sage: CDF.0
            1.0*I
            sage: CDF.gens()
            (1.0*I,)
        """
        if n != 0:
            raise ValueError, "only 1 generator"
        return I

    def ngens(self):
        return 1

cdef class ComplexDoubleElement(sage.ext.element.FieldElement):
    cdef gsl_complex _complex
    def __init__(self, real, imag):
        self._complex = gsl_complex_rect(real, imag)

    # TODO: real and imag should be elements of RealDoubleField, when that exists.
    def real(self):
        """
        Return the real part of this complex double.

        EXAMPLES:
            sage: a = CDF(3,-2)
            sage: a.real()
            3.0
        """
        return self._complex.dat[0]

    def imag(self):
        """
        Return the imaginary part of this complex double.

        EXAMPLES:
            sage: a = CDF(3,-2)
            sage: a.imag()
            -2.0
        """
        return self._complex.dat[1]

    def __complex__(self):
        """
        EXAMPLES:
            sage: a = complex(2303,-3939)
            sage: CDF(a)
            2303.0 - 3939.0*I
            sage: complex(CDF(a))
            (2303-3939j)
        """
        return complex(self._complex.dat[0], self._complex.dat[1])

    def parent(self):
        """
        Return the complex double field, which is the parent of self.

        EXAMPLES:
            sage: a = CDF(2,3)
            sage: a.parent()
            Complex Double Field
            sage: parent(a)
            Complex Double Field
        """
        return CDF

    def __repr__(self):
        """
        Return print version of self.

        EXAMPLES:
            sage: a = CDF(2,-3); a
            2.0 - 3.0*I
            sage: a^2
            -5.0 - 12.0*I
        """
        # todo -- redo completely in C
        cdef float y
        s = ""
        if self._complex.dat[0] != 0:
            s = str(self._complex.dat[0])
        y  =  self._complex.dat[1]
        if y != 0:
            if s != "":
                if y < 0:
                    s = s+" - "
                    y = -y
                else:
                    s = s+" + "
            s = s+"%s*I"%y
        if len(s) == 0:
            s = "0"
        return s

    def _add_(ComplexDoubleElement self, ComplexDoubleElement right):
        """
        Add self and right.

        EXAMPLES:
            sage: CDF(2,-3)._add_(CDF(1,-2))
            3.0 - 5.0*I
        """
        cdef ComplexDoubleElement z
        z = ComplexDoubleElement.__new__(ComplexDoubleElement)
        z._complex = gsl_complex_add(self._complex, right._complex)
        return z

    def __add__(x, y):
        try:
            return x._add_(y)
        except (TypeError, AttributeError):
            return sage.rings.coerce.bin_op(x, y, operator.add)
        #if isinstance(x, ComplexDoubleElement) and isinstance(y, ComplexDoubleElement):
        #    return x._add_(y)
        #return sage.rings.coerce.bin_op(x, y, operator.add)

    def _sub_(ComplexDoubleElement self, ComplexDoubleElement right):
        """
        Subtract self and right.

        EXAMPLES:
            sage: CDF(2,-3)._sub_(CDF(1,-2))
            1.0 - 1.0*I
        """
        cdef ComplexDoubleElement z
        z = ComplexDoubleElement.__new__(ComplexDoubleElement)
        z._complex = gsl_complex_sub(self._complex, right._complex)
        return z

    def __sub__(x, y):
        try:
            return x._sub_(y)
        except (TypeError, AttributeError):
            return sage.rings.coerce.bin_op(x, y, operator.sub)

    def _mul_(ComplexDoubleElement self, ComplexDoubleElement right):
        """
        Multiply self and right.

        EXAMPLES:
            sage: CDF(2,-3)._mul_(CDF(1,-2))
            -4.0 - 7.0*I
        """
        cdef ComplexDoubleElement z
        z = ComplexDoubleElement.__new__(ComplexDoubleElement)
        z._complex = gsl_complex_mul(self._complex, right._complex)
        return z

    def __mul__(x, y):
        try:
            return x._mul_(y)
        except (TypeError, AttributeError):
            return sage.rings.coerce.bin_op(x, y, operator.mul)

    def _div_(ComplexDoubleElement self, ComplexDoubleElement right):
        """
        Divide self by right.

        EXAMPLES:
            sage: CDF(2,-3)._div_(CDF(1,-2))
            1.6 + 0.20000000298*I
        """
        cdef ComplexDoubleElement z
        z = ComplexDoubleElement.__new__(ComplexDoubleElement)
        z._complex = gsl_complex_div(self._complex, right._complex)
        return z

    def __div__(x, y):
        try:
            return x._div_(y)
        except (TypeError, AttributeError):
            return sage.rings.coerce.bin_op(x, y, operator.div)

    def arg(self):
        r"""
        This function returns the argument of the complex number $z$, $\arg(z)$,
        where $-\pi < \arg(z) <= \pi$.

        EXAMPLES:
            sage: CDF(1,0).arg()
            0.0
            sage: CDF(0,1).arg()
            1.5707963267948966
            sage: CDF(0,-1).arg()
            -1.5707963267948966
            sage: CDF(-1,0).arg()
            3.1415926535897931
        """
        return gsl_complex_arg(self._complex)

    def __abs__(self):
        """
        This function returns the magnitude of the complex number $z$, $|z|$.

        EXAMPLES:
            sage: abs(CDF(1,2))
            2.2360679774997898
            sage: abs(CDF(1,0))
            1.0
            sage: abs(CDF(-2,3))
            3.6055512754639891
        """
        return gsl_complex_abs(self._complex)

    def abs(self):
        """
        This function returns the magnitude of the complex number $z$, $|z|$.

        EXAMPLES:
            sage: CDF(2,3).abs()
            3.6055512754639891
        """
        return gsl_complex_abs(self._complex)

    def abs2(self):
        """
        This function returns the squared magnitude of the complex number $z$, $|z|^2$.

        EXAMPLES:
            sage: CDF(2,3).abs2()
            13.0
        """
        return gsl_complex_abs2(self._complex)

    def logabs(self):
        r"""
        This function returns the natural logarithm of the magnitude of the complex
        number $z$, $\log|z|$.

        This allows for an accurate evaluation of $\log|z| when $|z|$
        is close to $1$.  The direct evaluation of \code{log(abs(z))}
        would lead to a loss of precision in this case.

        EXAMPLES:
        We try it out.
            sage: CDF(1.1,0.1).logabs()
            0.099425429372582669
            sage: log(abs(CDF(1.1,0.1)))
            0.099425429373735899

        Which is better?
            sage: log(abs(ComplexField(200)(1.1,0.1)))
            0.099425429372582675602989386713555936556752871164033127857197658

        Indeed, the logabs, wins.
        """
        return gsl_complex_logabs(self._complex)



#####################################################
# unique objects
#####################################################
ComplexDoubleField = ComplexDoubleField_class()
CDF = ComplexDoubleField
I = ComplexDoubleElement(0,1)

