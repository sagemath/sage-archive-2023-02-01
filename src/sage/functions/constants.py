r"""
Mathematical constants

The following standard mathematical constants are defined in \sage,
along with support for coercing them into GAP, GP/PARI, KASH, Maxima,
Mathematica, Maple, Octave, and Singular:

    sage: pi
    pi
    sage: e             # base of the natural logarithm
    e
    sage: NaN           # Not a number
    NaN
    sage: golden_ratio
    golden_ratio
    sage: log2          # natural logarithm of the real number 2
    log2
    sage: euler_gamma   # Euler's gamma constant
    euler_gamma
    sage: catalan       # the Catalon constant
    catalan
    sage: khinchin      # Khinchin's constant
    khinchin
    sage: twinprime
    twinprime
    sage: merten
    merten
    sage: brun
    brun

Support for coercion into the various systems means that if, e.g.,
you want to create $\pi$ in Maxima and Singular, you don't have
to figure out the special notation for each system.  You just
type the following:

    sage: pi.str()
    '3.14159265358979'
    sage: maxima(pi)
    %pi
    sage: singular(pi)
    3.14159265358979
    sage: gap(pi)
    3.14159265358979
    sage: gp(pi)
    3.141592653589793238462643383     # 32-bit
    3.1415926535897932384626433832795028842   # 64-bit
    sage: pari(pi)
    3.141592653589793238462643383     # 32-bit
    3.1415926535897932384626433832795028842   # 64-bit
    sage: kash(pi)                    # optional
    3.14159265358979323846264338328
    sage: mathematica(pi)             # optional
    Pi
    sage: maple(pi)                   # optional
    Pi
    sage: octave(pi)                  # optional
    3.14159

Arithmetic operations with constants also yield constants, which
can be coerced into other systems or evaluated.
    sage: a = pi + e*4/5; a
    pi + 4*e/5
    sage: maxima(a)
    %pi+4*%e/5
    sage: RealField(15)(a)           # 15 *bits* of precision
    5.316
    sage: gp(a)
    5.316218116357029426750873360              # 32-bit
    5.3162181163570294267508733603616328824    # 64-bit
    sage: print mathematica(a)                     # optional
     4 E
     --- + Pi
      5

EXAMPLES: Decimal expansions of constants

We can obtain floating point approximations to each of these constants
by coercing into the real field with given precision.  For example, to
200 decimal places we have the following:

    sage: R = RealField(200); R
    Real Field with 200 bits of precision

    sage: R(pi)
    3.1415926535897932384626433832795028841971693993751058209749

    sage: R(e)
    2.7182818284590452353602874713526624977572470936999595749670

    sage: R(NaN)
    NaN

    sage: R(golden_ratio)
    1.6180339887498948482045868343656381177203091798057628621354

    sage: R(log2)
    0.69314718055994530941723212145817656807550013436025525412068

    sage: R(euler_gamma)
    0.57721566490153286060651209008240243104215933593992359880577

    sage: R(catalan)
    0.91596559417721901505460351493238411077414937428167213426650

    sage: R(khinchin)
    2.6854520010653064453097148354817956938203822939944629530512


EXAMPLES: Arithmetic with constants
    sage: f = I*(e+1); f
    (e + 1)*I
    sage: f^2
    -(e + 1)^2

    sage: pp = pi+pi; pp
    2*pi
    sage: R(pp)
    6.2831853071795864769252867665590057683943387987502116419499

    sage: s = (1 + e^pi); s
    e^pi + 1
    sage: R(s)
    24.140692632779269005729086367948547380266106242600211993445
    sage: R(s-1)
    23.140692632779269005729086367948547380266106242600211993445

    sage: l = (1-log2)/(1+log2); l
    (1 - log(2))/(log(2) + 1)
    sage: R(l)
    0.18123221829928249948761381864650311423330609774776013488056

    sage: pim = maxima(pi)
    sage: maxima.eval('fpprec : 100')
    '100'
    sage: pim.bfloat()
    3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068b0



AUTHORS:
    -- Alex Clemesha  <aclemesh@ucsd.edu>, 2006-01-15
    -- William Stein
    -- Alex Clemesha & William Stein (2006-02-20): added new constants; removed todos
    -- didier deshommes <dfdeshom@gmail.com> (2007-03-27): added constants from RQDF

TESTS:
    Coercion of each constant to the RQDF:
        sage: RQDF(e)
        2.718281828459045235360287471352662497757247093699959574966967630
        sage: RQDF(pi)
        3.141592653589793238462643383279502884197169399375105820974944590
        sage: RQDF(e)
        2.718281828459045235360287471352662497757247093699959574966967630
        sage: RQDF(I)
        Traceback (most recent call last):
        ...
        TypeError
        sage: RQDF(golden_ratio)
        1.618033988749894848204586834365638117720309179805762862135448623
        sage: RQDF(log2)
        0.693147180559945309417232121458176568075500134360255254120680009
        sage: RQDF(euler_gamma)
        0.577215664901532860606512090082402431042159335939923598805767234
        sage: RQDF(catalan)
        0.915965594177219015054603514932384110774149374281672134266498119
        sage: RQDF(khinchin)
        2.685452001065306445309714835481795693820382293994462953051152345
        sage: RQDF(twinprime)
        0.660161815846869573927812110014555778432623360284733413319448422
        sage: RQDF(merten)
        0.261497212847642783755426838608695859051566648261199206192064212
        sage: RQDF(brun)
        Traceback (most recent call last):
        ...
        NotImplementedError: Brun's constant only available up to 41 bits


Coercing the sum of a bunch of the constants to many different
floating point rings:
    sage: a = pi + e + golden_ratio + log2 + euler_gamma + catalan + khinchin + twinprime + merten; a
    twinprime + merten + khinchin + euler_gamma + catalan + log(2) + pi + e + (sqrt(5) + 1)/2
    sage: parent(a)
    Symbolic Ring
    sage: RQDF(a)
    13.27134794019724931009881919957581394087110682000307481783297119
    sage: RR(a)
    13.2713479401972
    sage: RealField(212)(a)
    13.2713479401972493100988191995758139408711068200030748178329712
    sage: RealField(230)(a)
    13.271347940197249310098819199575813940871106820003074817832971189555
    sage: CC(a)
    13.2713479401972
    sage: CDF(a)
    13.2713479402
    sage: ComplexField(230)(a)
    13.271347940197249310098819199575813940871106820003074817832971189555
    sage: RDF(a)
    13.2713479402
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 Alex Clemesha <aclemesh@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import math, operator

from sage.rings.all import CommutativeRing, RealField, Integer, RingElement, QQ
import sage.interfaces.all
import sage.rings.all
from sage.libs.pari.all import pari
from sage.misc.latex import latex


# Default real field used for coercion.

from functions import Function_gen, Function_arith, Function, FunctionRing_class


######################
# Ring of Constants
######################

## class ConstantRing_class(FunctionRing_class):
##     def _repr_(self):
##         return "Ring of Real Mathematical Constants"

##     def __cmp__(self, right):
##         if isinstance(right, ConstantRing_class):
##             return 0
##         return -1

##     def __call__(self, x):
##         try:
##             return self._coerce_(x)
##         except TypeError:
##             return Constant_gen(x)

##     def _coerce_impl(self, x):
##         if isinstance(x, (sage.rings.integer.Integer,
##                           sage.rings.rational.Rational, int, long)):
##             return Constant_gen(x)
##         raise TypeError, 'no canonical coercion of element into self.'

## ConstantRing = ConstantRing_class()


######################
# Constant functions
######################
import sage.calculus.calculus
SR = sage.calculus.calculus.SR

class Constant(Function):
    def __init__(self, conversions={}, parent=sage.calculus.calculus.SR):
        self._conversions = conversions
        RingElement.__init__(self, parent)

    def number_of_arguments(self):
        """
        Returns the number of arguments of this constant, viewed as a function.
        This is of course always 0.

        EXAMPLES:
            sage: pi.number_of_arguments()
            0
        """
        return 0

    # The maxima one is special:
    def _maxima_(self, session=None):
        """
        Returns self as a maxima object.

        EXAMPLES:
            sage: pi._maxima_()
            %pi
        """
        if session is None:
            from sage.calculus.calculus import maxima
            return RingElement._maxima_(self, maxima)
        else:
            return RingElement._maxima_(self, session)

    def _has_op(self, x):
        """
        Check whether or not self contains the operation x.  Since
        self is a constant, this is always False.

        EXAMPLES:
            sage: pi._has_op(operator.add)
            False
        """
        return False

    def number_of_arguments(self):
        """
        Returns the number of arguments of self.  For constants,
        this is just zero.

        EXAMPLES:
            sage: type(pi)
            <class 'sage.functions.constants.Pi'>
            sage: pi.number_of_arguments()
            0
            sage: e.number_of_arguments()
            0
        """
        return 0

    def substitute(self, *args, **kwds):
        """
        Substitute values into self.  For constants, this just returns
        self.

        EXAMPLES:
            sage: pi.substitute(x=3)
            pi
            sage: pi.substitute(3)
            pi
            sage: pi.substitute(4, x=4)
            pi

        """
        return self

    def _recursive_sub(self, kwds):
        """
        Recursively substitute values into self.  For constants, this just
        returns self.

        EXAMPLES:
            sage: pi._recursive_sub({x:3})
            pi
        """
        return self

    def _recursive_sub_over_ring(self, kwds, ring):
        """
        Recursively substitute values into self over a ring.
        For constants, this just returns ring(self).

        EXAMPLES:
            sage: pi._recursive_sub_over_ring({x:3}, RDF)
            3.14159265359
        """
        return ring(self)

    def variables(self):
        """
        Return a list of the variables of self.  For constants, this
        is the empty list.

        EXAMPLES:
            sage: pi.variables()
            []
            sage: e.variables()
            []
        """
        return []

    def _ser(self):
        """
        Returns self as an element of SymbolicRing.

        EXAMPLES:
            sage: s = pi._ser(); s
            pi
            sage: type(s)
            <class 'sage.calculus.calculus.SymbolicConstant'>
            sage: s.parent()
            Symbolic Ring
        """
        try:
            return self.__ser
        except AttributeError:
            self.__ser = sage.calculus.calculus.SR._coerce_impl(self)
            return self.__ser

    def __abs__(self):
        """
        Returns the absolute value of self.

        EXAMPLES:
           sage: abs(pi)
           pi
        """
        if self.str()[0] != '-':
            return self
        return -self

    def _neg_(self):
        """
        Returns the negation of self.

        EXAMPLES:
            sage: -pi
            -1*pi
            sage: -e
            -1*e
        """
        return -Integer(1)*self

    def __call__(self, *args, **kwargs):
        """
        Call self as a function.  Since self is a constant function,
        this just returns self.

        EXAMPLES:
            sage: pi(2,3,4)
            pi
            sage: pi(pi)
            pi
        """
        return self

    def _fast_float_(self, *vars):
        from sage.ext.fast_eval import fast_float_constant
        return fast_float_constant(self)

    def floor(self):
        """
        Returns the floor of self.

        EXAMPLES:
            sage: pi.floor()
            3
            sage: e.floor()
            2
            sage: golden_ratio.floor()
            1
            sage: log2.floor()
            0
        """
        return Integer(int(float(self)))

    def _latex_(self):
        r"""
        Return the \LaTeX representation of self.

        EXAMPLES:
            sage: catalan._latex_()
            '\\text{catalan}'
        """
        return '\\text{%s}'%self

    def _complex_mpfr_field_(self, R):
        """
        EXAMPLES:
            sage: pi._complex_mpfr_field_(ComplexField(53))
            3.14159265358979
            sage: pi._complex_mpfr_field_(ComplexField(200))
            3.1415926535897932384626433832795028841971693993751058209749
        """
        return R(self._mpfr_(R._real_field()))

    def _real_double_(self, R):
        """
        EXAMPLES:
           sage: pi._real_double_(RDF)
           3.14159265359
        """
        return R(float(self))

    def _complex_double_(self, R):
        """
        EXAMPLES:
            sage: pi._complex_double_(CDF)
            3.14159265359
        """
        return R(float(self))

    # The following adds formal arithmetic support for generic constant
    def _add_(self, right):
        """
        EXAMPLES:
            sage: I + 2
            I + 2
            sage: a = I+I
            sage: map(type, a._operands)
            [<class 'sage.calculus.calculus.SymbolicConstant'>,
             <class 'sage.calculus.calculus.SymbolicConstant'>]
        """
        return self._ser() + SR(right)

    def _sub_(self, right):
        """
        EXAMPLES:
            sage: a = I - pi; a
            I - pi
            sage: map(type, a._operands)
            [<class 'sage.calculus.calculus.SymbolicConstant'>,
             <class 'sage.calculus.calculus.SymbolicConstant'>]
        """
        return self._ser() - SR(right)

    def _mul_(self, right):
        """
        EXAMPLES:
            sage: a = I * pi; a
            I*pi
            sage: map(type, a._operands)
            [<class 'sage.calculus.calculus.SymbolicConstant'>,
             <class 'sage.calculus.calculus.SymbolicConstant'>]
        """
        return self._ser() * SR(right)

    def _div_(self, right):
        """
        EXAMPLES:
            sage: a = I / pi; a
            I/pi
            sage: map(type, a._operands)
            [<class 'sage.calculus.calculus.SymbolicConstant'>,
             <class 'sage.calculus.calculus.SymbolicConstant'>]

        """
        return self._ser() / SR(right)

    def __pow__(self, right):
        """
        EXAMPLES:
            sage: a = pi^pi; a
            pi^pi
            sage: map(type, a._operands)
            [<class 'sage.calculus.calculus.SymbolicConstant'>,
             <class 'sage.calculus.calculus.SymbolicConstant'>]
        """
        return self._ser() ** SR(right)

    def _interface_is_cached_(self):
        """
        Return False, since coercion of functions to interfaces
        is not cached.

        We do not cache coercions of functions to interfaces, since
        the precision of the interface may change.

        EXAMPLES:
            sage: gp(pi)
            3.141592653589793238462643383              # 32-bit
            3.1415926535897932384626433832795028842    # 64-bit
            sage: old_prec = gp.set_precision(100)
            sage: gp(pi)
            3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
            sage: _ = gp.set_precision(old_prec)
            sage: gp(pi)
            3.141592653589793238462643383              # 32-bit
            3.1415926535897932384626433832795028842    # 64-bit
        """
        return False

    def __lt__(self, right):
        """
        EXAMPLES:
            sage: pi < 3
            pi < 3
            sage: type(pi<3)
            <class 'sage.calculus.equations.SymbolicEquation'>
            sage: bool(pi<3)
            False
        """
        return self._ser().__lt__(right)

    def __le__(self, right):
        """
        EXAMPLES:
            sage: pi <= 3
            pi <= 3
            sage: type(pi<=3)
            <class 'sage.calculus.equations.SymbolicEquation'>
            sage: bool(pi<=3)
            False
        """
        return self._ser().__le__(right)

    def __eq__(self, right):
        """
        EXAMPLES:
            sage: solve(pi == 2*x)
            [x == pi/2]
            sage: solve(cos(x^2) == pi)
            [x == -sqrt(arccos(pi)), x == sqrt(arccos(pi))]
        """
        return self._ser().__eq__(right)

    def __ne__(self, right):
        """
        EXAMPLES:
            sage: type(pi != 3)
            <class 'sage.calculus.equations.SymbolicEquation'>
            sage: bool(pi != 3)
            True
        """
        return self._ser().__ne__(right)

    def __ge__(self, right):
        """
        EXAMPLES:
            sage: type(pi>=3)
            <class 'sage.calculus.equations.SymbolicEquation'>
            sage: bool(pi>=3)
            True
        """
        return self._ser().__ge__(right)

    def __gt__(self, right):
        """
        EXAMPLES:
            sage: type(pi>3)
            <class 'sage.calculus.equations.SymbolicEquation'>
            sage: bool(pi>3)
            True
        """
        return self._ser().__gt__(right)



class Constant_gen(Constant, Function_gen):
    def __init__(self, x):
        Function_gen.__init__(self, x)
        Constant.__init__(self)

    def __call__(self, x):
        return self.obj()

    def __repr__(self):
        return Function_gen._repr_(self)

class Constant_arith(Constant, Function_arith):
    def __init__(self, x, y, op):
        Function_arith.__init__(self, x, y, op)
        Constant.__init__(self)

class Pi(Constant):
    """
    The ratio of a circle's circumference to its diameter.

    EXAMPLES:
        sage: pi
        pi
        sage: float(pi)
        3.1415926535897931
        sage: gp(pi)
        3.141592653589793238462643383            # 32-bit
        3.1415926535897932384626433832795028842  # 64-bit
        sage: RR(pi)
        3.14159265358979
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(pi)
        3.1415926535897932384626433832795028841971693993751058209749
        sage: pp = pi+pi; pp
        2*pi
        sage: R(pp)
        6.2831853071795864769252867665590057683943387987502116419499
        sage: maxima(pi)
        %pi
        sage: maxima(pi).float()
        3.141592653589793
    """
    def __init__(self):
        """
        EXAMPLES:
            sage: bool(pi == loads(dumps(pi)))
            True
        """
        Constant.__init__(self,
            {'axiom':'%pi',
             'maxima':'%pi','gp':'Pi','kash':'PI','mathematica':'Pi',
             'matlab':'pi','maple':'Pi','octave':'pi','pari':'Pi'})

    def _repr_(self, simplify=True):
        """
        EXAMPLES:
            sage: repr(pi)
            'pi'
        """
        return "pi"

    def _latex_(self):
        """
        EXAMPLES:
            sage: pi._latex_()
            '\\pi'
            sage: latex(pi)
            \pi
        """
        return "\\pi"

    def _mathml_(self):
        """
        EXAMPLES:
            sage: pi._mathml_()
            '<mi>&pi;</mi>'
            sage: mathml(pi)
            <mi>&pi;</mi>
        """
        return "<mi>&pi;</mi>"

    def __float__(self):
        """
        EXAMPLES:
            sage: float(pi)
            3.1415926535897931
        """
        return math.pi

    def _mpfr_(self, R):
        """
        EXAMPLES:
            sage: pi._mpfr_(RealField(100))
            3.1415926535897932384626433833
        """
        return R.pi()

    def _real_double_(self, R):
        """
        EXAMPLES:
            sage: pi._real_double_(RDF)
            3.14159265359
         """
        return R.pi()

    def _real_rqdf_(self, R):
        """
        EXAMPLES:
            sage: pi._real_rqdf_(RQDF)
            3.141592653589793238462643383279502884197169399375105820974944590
        """
        return R.pi()


    def floor(self):
        """
        Returns the floor of self.

        EXAMPLES:
            sage: pi.floor()
            3
        """
        return Integer(3)

    # This just gives a string in singular anyways, and it's
    # *REALLY* slow!
    #def _singular_(self, singular):
    #    singular.lib('general')
    #    return singular('number_pi(%s)'%int((ConstantRing._default_precision/3)+1))

pi = Pi()

python_complex_i = complex(0,1)

class I_class(Constant):
    """
    The formal square root of -1.

    EXAMPLES:
        sage: I
        I
        sage: I^2
        -1
        sage: float(I)
        Traceback (most recent call last):
        ...
        TypeError
        sage: gp(I)
        I
        sage: RR(I)
        Traceback (most recent call last):
        ...
        TypeError
        sage: C = ComplexField(200); C
        Complex Field with 200 bits of precision
        sage: C(I)
        1.0000000000000000000000000000000000000000000000000000000000*I
        sage: z = I + I; z
        2*I
        sage: C(z)
        2.0000000000000000000000000000000000000000000000000000000000*I
        sage: maxima(2*I)
        2*%i
        sage: 1e8*I
        1.00000000000000e8*I
    """
    def __init__(self):
        """
        EXAMPLES:
            sage: bool(I == loads(dumps(I)))
            True
        """
        Constant.__init__(self,
            {'axiom':'%i',
             'maxima':'%i','gp':'I','mathematica':'I',
             'matlab':'i','maple':'I','octave':'i','pari':'I'})

    def _repr_(self, simplify=True):
        """
        EXAMPLES:
            sage: repr(I)
            'I'
        """
        return "I"

    def _latex_(self):
        """
        EXAMPLES:
            sage: I._latex_()
            'i'
            sage: latex(I)
            i
        """
        return "i"

    def minpoly(self, bits=None, degree=None, epsilon=0):
        """
        EXAMPLES:
            sage: I.minpoly()
            x^2 + 1
        """
        return QQ['x'].gen(0)**2 + 1

    def _mathml_(self):
        """
        EXAMPLES:
            sage: I._mathml_()
            '<mi>&i;</mi>'
            sage: mathml(I)
            <mi>&i;</mi>
        """
        return "<mi>&i;</mi>"

    def __float__(self):
        """
        EXAMPLES:
            sage: float(I)
            Traceback (most recent call last):
            ...
            TypeError
        """
        raise TypeError

    def _mpfr_(self, R):
        """
        EXAMPLES:
            sage: I._mpfr_(RealField(53))
            Traceback (most recent call last):
            ...
            TypeError
        """
        raise TypeError

    def _real_rqdf_(self, R):
        """
        EXAMPLES:
            sage: I._real_rqdf_(RQDF)
            Traceback (most recent call last):
            ...
            TypeError
        """
        raise TypeError

    def _complex_mpfr_field_(self, R):
        """
        EXAMPLES:
            sage: I._complex_mpfr_field_(ComplexField(53))
            1.00000000000000*I
        """
        return R.gen()

    def _complex_double_(self, C):
        """
        EXAMPLES:
            sage: I._complex_double_(CDF)
            1.0*I
        """
        return C.gen()

    def __complex__(self):
        """
        EXAMPLES:
            sage: complex(I)
            1j
        """
        return python_complex_i

    def _real_double_(self, R):
        """
        EXAMPLES:
            sage: I._mpfr_(RealField(53))
            Traceback (most recent call last):
            ...
            TypeError
        """
        raise TypeError

    def _algebraic_(self, field):
        """
        EXAMPLES:
            sage: QQbar(I)
            1*I
        """
        import sage.rings.qqbar
        return field(sage.rings.qqbar.QQbar_I)

    def __abs__(self):
        """
        EXAMPLES:
            sage: abs(I)
            1
            sage: I.__abs__()
            1
        """
        return Integer(1)

    def floor(self):
        """
        EXAMPLES:
            sage: I.floor()
            Traceback (most recent call last):
            ...
            TypeError
        """
        raise TypeError

I = I_class()

class E(Constant):
    """
    The base of the natural logarithm.

    EXAMPLES:
        sage: RR(e)
        2.71828182845905
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(e)
        2.7182818284590452353602874713526624977572470936999595749670
        sage: em = 1 + e^(1-e); em
        e^(1 - e) + 1
        sage: R(em)
        1.1793740787340171819619895873183164984596816017589156131574
        sage: maxima(e).float()
        2.718281828459045
        sage: t = mathematica(e)               # optional
        sage: t                                # optional
        E
        sage: float(t)                         # optional
        2.7182818284590451
    """
    def __init__(self):
        """
        EXAMPLES:
            sage: bool( e == loads(dumps(e)) )
            True
        """
        Constant.__init__(self,
            {'axiom':'%e',
             'maxima':'%e',
             'gp':'exp(1)',
             'kash':'E',
             'pari':'exp(1)',
             'mathematica':'E',
             'maple':'exp(1)',
             'octave':'e'})

    def _repr_(self, simplify=True):
        """
        EXAMPLES:
            sage: repr(e)
            'e'
            sage: e._repr_()
            'e'
        """
        return 'e'

    def _latex_(self):
        """
        EXAMPLES:
            sage: e._latex_()
            'e'
            sage: latex(e)
            e
        """
        return 'e'

    def __float__(self):
        """
        EXAMPLES:
            sage: float(e)
            2.7182818284590451
            sage: e.__float__()
            2.7182818284590451
        """
        return math.e

    def _mpfr_(self, R):
        """
        EXAMPLES:
            sage: e._mpfr_(RealField(100))
            2.7182818284590452353602874714
        """
        return R(1).exp()

    def floor(self):
        """
        Returns the floor of self.

        EXAMPLES:
            sage: e.floor()
            2
        """
        return Integer(2)

    def _real_double_(self, R):
        """
        EXAMPLES:
            sage: e._real_double_(RDF)
            2.71828182846
        """
        return R(1).exp()

    def _real_rqdf_(self, R):
        """
        EXAMPLES:
            sage: RQDF = RealQuadDoubleField ()
            sage: RQDF.e()
            2.718281828459045235360287471352662497757247093699959574966967630
        """
        return R.e()

    # This just gives a string in singular anyways, and it's
    # *REALLY* slow!
    #def _singular_(self, singular):
    #    singular.lib('general')
    #    return singular('number_e(%s)'%int((ConstantRing._default_precision/3)+1))

e = E()
ee = e

class NotANumber(Constant):
    """
    Not a Number
    """
    def __init__(self):
        """
        EXAMPLES:
            sage: loads(dumps(NaN))
            NaN
        """
        Constant.__init__(self,
	    {'matlab':'NaN'})

    def _repr_(self, simplify=True):
        """
        EXAMPLES:
            sage: repr(NaN)
            'NaN'
            sage: NaN._repr_()
            'NaN'
        """
        return 'NaN'

    def _mpfr_(self,R):
        """
        EXAMPLES:
            sage: NaN._mpfr_(RealField(53))
            NaN
            sage: type(_)
            <type 'sage.rings.real_mpfr.RealNumber'>
        """
        return R('NaN') #??? nan in mpfr: void mpfr_set_nan (mpfr_t x)

    def _real_double_(self, R):
        """
        EXAMPLES:
            sage: RDF(NaN)
            nan
        """
        return R.NaN()

    def _real_rqdf_(self, R):
        """
        EXAMPLES:
            sage: RQDF(NaN)
            'NaN'
        """
        return R.NaN()

NaN = NotANumber()

class GoldenRatio(Constant):
    """
    The number (1+sqrt(5))/2

    EXAMPLES:
        sage: gr = golden_ratio
        sage: RR(gr)
        1.61803398874989
        sage: R = RealField(200)
        sage: R(gr)
        1.6180339887498948482045868343656381177203091798057628621354
        sage: grm = maxima(golden_ratio);grm
        (sqrt(5)+1)/2
        sage: grm + grm
        sqrt(5)+1
        sage: float(grm + grm)
        3.2360679774997898
    """
    def __init__(self):
        """
        EXAMPLES:
            sage: loads(dumps(golden_ratio))
            golden_ratio
        """
        Constant.__init__(self,{'mathematica':'N[(1+Sqrt[5])/2]','gp':'(1+sqrt(5))/2',
				'maple':'(1+sqrt(5))/2','maxima':'(1+sqrt(5))/2',
				'pari':'(1+sqrt(5))/2','octave':'(1+sqrt(5))/2',
				'kash':'(1+Sqrt(5))/2'})
    def _repr_(self, simplify=True):
        """
        EXAMPLES:
            sage: repr(golden_ratio)
            'golden_ratio'
            sage: golden_ratio._repr_()
            'golden_ratio'
        """
        return 'golden_ratio'

    def _latex_(self):
        """
        EXAMPLES:
            sage: latex(golden_ratio)
            \phi
            sage: golden_ratio._latex_()
            '\\phi'
        """
        return '\\phi'

    def minpoly(self, bits=None, degree=None, epsilon=0):
        """
        EXAMPLES:
            sage: golden_ratio.minpoly()
            x^2 - x - 1
        """
        x = QQ['x'].gen(0)
        return x**2 - x - 1

    def __float__(self):
        """
        EXAMPLES:
            sage: float(golden_ratio)
            1.6180339887498949
            sage: golden_ratio.__float__()
            1.6180339887498949
        """
        return float(0.5)*(float(1.0)+math.sqrt(float(5.0)))

    def _real_double_(self, R):
        """
        EXAMPLES:
            sage: RDF(golden_ratio)
            1.61803398875
        """
        return R('1.61803398874989484820458')

    def _real_rqdf_(self, R):
        """
        EXAMPLES:
            sage: golden_ratio._real_rqdf_(RQDF)
            1.618033988749894848204586834365638117720309179805762862135448623
            sage: RQDF(golden_ratio)
            1.618033988749894848204586834365638117720309179805762862135448623
        """
        return (R(1)+R(5).sqrt())/R(2)

    def _mpfr_(self,R):
        """
        EXAMPLES:
            sage: golden_ratio._mpfr_(RealField(100))
            1.6180339887498948482045868344
            sage: RealField(100)(golden_ratio)
            1.6180339887498948482045868344
        """
	return (R(1)+R(5).sqrt())/R(2)

    def _algebraic_(self, field):
        """
        EXAMPLES:
            sage: golden_ratio._algebraic_(QQbar)
            1.618033988749895?
            sage: QQbar(golden_ratio)
            1.618033988749895?
        """
        import sage.rings.qqbar
        return field(sage.rings.qqbar.get_AA_golden_ratio())

golden_ratio = GoldenRatio()

class Log2(Constant):
    """
    The natural logarithm of the real number 2.

    EXAMPLES:
        sage: log2
        log2
        sage: float(log2)
        0.69314718055994529
        sage: RR(log2)
        0.693147180559945
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(log2)
        0.69314718055994530941723212145817656807550013436025525412068
        sage: l = (1-log2)/(1+log2); l
        (1 - log(2))/(log(2) + 1)
        sage: R(l)
        0.18123221829928249948761381864650311423330609774776013488056
        sage: maxima(log2)
        log(2)
        sage: maxima(log2).float()
        .6931471805599453
        sage: gp(log2)
        0.6931471805599453094172321215             # 32-bit
        0.69314718055994530941723212145817656807   # 64-bit
    """
    def __init__(self):
        """
        EXAMPLES:
            sage: loads(dumps(log2))
            log2
        """
        Constant.__init__(self,{'mathematica':'N[Log[2]]','kash':'Log(2)',
				'maple':'log(2)','maxima':'log(2)','gp':'log(2)',
				'pari':'log(2)','octave':'log(2)'})

    def _repr_(self, simplify=True):
        """
        EXAMPLES:
            sage: repr(log2)
            'log2'
            sage: log2._repr_()
            'log2'
        """
        return 'log2'

    def _latex_(self):
        """
        EXAMPLES:
            sage: log2._latex_()
            '\\log(2)'
            sage: latex(log2)
            \log(2)
        """
        return '\\log(2)'

    def __float__(self):
        """
        EXAMPLES:
            sage: float(log2)
            0.69314718055994529
            sage: log2.__float__()
            0.69314718055994529
        """
        return math.log(2)

    def _real_double_(self, R):
        """
        EXAMPLES:
            sage: RDF(log2)
            0.69314718056
        """
        return R.log2()

    def _real_rqdf_(self, R):
        """
        EXAMPLES:
            sage: RQDF(log2)
            0.693147180559945309417232121458176568075500134360255254120680009
        """
        return R.log2()

    def _mpfr_(self,R):
        """
        EXAMPLES:
            sage: RealField(100)(log2)
            0.69314718055994530941723212146
            sage: log2._mpfr_(RealField(100))
            0.69314718055994530941723212146
        """
        return R.log2()

    def floor(self):
        """
        Returns the floor of self.

        EXAMPLES:
            sage: log2.floor()
            0
        """
        return Integer(0)

log2 = Log2()

class EulerGamma(Constant):
    """
    The limiting difference between the harmonic series and the natural logarithm.

    EXAMPLES:
        sage: R = RealField()
        sage: R(euler_gamma)
        0.577215664901533
        sage: R = RealField(200); R
        Real Field with 200 bits of precision
        sage: R(euler_gamma)
        0.57721566490153286060651209008240243104215933593992359880577
        sage: eg = euler_gamma + euler_gamma; eg
        2*euler_gamma
        sage: R(eg)
        1.1544313298030657212130241801648048620843186718798471976115
    """
    def __init__(self):
        """
        EXAMPLES:
            sage: loads(dumps(euler_gamma))
            euler_gamma
        """
        Constant.__init__(self,
	    {'kash':'EulerGamma(R)','maple':'gamma',
             'mathematica':'EulerGamma','pari':'Euler',
             'maxima':'%gamma', 'maxima':'euler_gamma'})

    def _repr_(self, simplify=True):
        """
        EXAMPLES:
            sage: repr(euler_gamma)
            'euler_gamma'
            sage: euler_gamma._repr_()
            'euler_gamma'
        """
        return 'euler_gamma'

    def _latex_(self):
        """
        EXAMPLES:
            sage: euler_gamma._latex_()
            '\\gamma'
            sage: latex(euler_gamma)
            \gamma
        """
        return '\\gamma'

    def _mpfr_(self,R):
        """
        EXAMPLES:
            sage: RealField(100)(euler_gamma)
            0.57721566490153286060651209008
            sage: euler_gamma._mpfr_(RealField(100))
            0.57721566490153286060651209008
        """
        return R.euler_constant()

    def _real_double_(self, R):
        """
        EXAMPLES:
            sage: RDF(euler_gamma)
            0.577215664902
        """
        return R.euler_constant()

    def _real_rqdf_(self, R):
        """
        EXAMPLES:
            sage: RQDF(euler_gamma)
            0.577215664901532860606512090082402431042159335939923598805767234
            sage: euler_gamma._real_rqdf_(RQDF)
            0.577215664901532860606512090082402431042159335939923598805767234
        """
        return R('0.577215664901532860606512090082402431042159335939923598805767235')

    def floor(self):
        """
        Return the floor of self.

        EXAMPLES:
            sage: euler_gamma.floor()
            0
        """
        return Integer(0)

euler_gamma = EulerGamma()

class Catalan(Constant):
    """
    A number appaering in combinatorics defined as the Dirichlet beta
    function evaluated at the number 2.

    EXAMPLES:
        sage: catalan^2 + merten
        merten + catalan^2
    """
    def __init__(self):
        """
        EXAMPLES:
            sage: loads(dumps(catalan))
            catalan
        """
        Constant.__init__(self,
             {'mathematica':'Catalan','kash':'Catalan(R)', #kash: R is default prec
              'maple':'Catalan', 'maxima':'catalan'})

    def _repr_(self, simplify=True):
        """
        EXAMPLES:
            sage: repr(catalan)
            'catalan'
            sage: catalan._repr_(catalan)
            'catalan'
        """
        return 'catalan'

    def _mpfr_(self, R):
        """
        EXAMPLES:
            sage: RealField(100)(catalan)
            0.91596559417721901505460351493
            sage: catalan._mpfr_(RealField(100))
            0.91596559417721901505460351493
        """
        return R.catalan_constant()

    def _real_double_(self, R):
        """
        EXAMPLES:
        We coerce to the real double field:
            sage: RDF(catalan)
            0.915965594177
        """
        return R('0.91596559417721901505460351493252')

    def _real_rqdf_(self, R):
        """
        EXAMPLES:
            sage: RQDF(catalan)
            0.915965594177219015054603514932384110774149374281672134266498119
        """
        return R('0.915965594177219015054603514932384110774149374281672134266498120')

    def __float__(self):
        """
        EXAMPLES:
            sage: float(catalan)
            0.91596559417721901
        """
        return 0.91596559417721901505460351493252

    def floor(self):
        """
        Return the floor of self.

        EXAMPLES:
            sage: catalan.floor()
            0
        """
        return Integer(0)

catalan = Catalan()

class Khinchin(Constant):
    """
    The geometric mean of the continued fraction expansion
    of any (almost any) real number.

    EXAMPLES:
        sage: float(khinchin)
        2.6854520010653062
        sage: khinchin.str(100)
        '2.6854520010653064453097148355'
        sage: m = mathematica(khinchin); m             # optional
        Khinchin
        sage: m.N(200)                                 # optional
        2.68545200106530644530971483548179569382038229399446295305115234555721885953715200280114117493184769799515346590528809008289767771641096305179253348325966838185231542133211949962603932852204481940961807          # 32-bit
        2.6854520010653064453097148354817956938203822939944629530511523455572188595371520028011411749318476979951534659052880900828976777164109630517925334832596683818523154213321194996260393285220448194096181                # 64-bit
    """
    def __init__(self):
        """
        EXAMPLES:
            sage: loads(dumps(khinchin))
            khinchin
        """
        Constant.__init__(self,
             {'maxima':'khinchin', 'mathematica':'Khinchin'}) #Khinchin is only implemented in Mathematica

        # digits come from http://pi.lacim.uqam.ca/piDATA/khintchine.txt
        self.__value = "2.6854520010653064453097148354817956938203822939944629530511523455572188595371520028011411749318476979951534659052880900828976777164109630517925334832596683818523154213321194996260393285220448194096180686641664289308477880620360737053501033672633577289049904270702723451702625237023545810686318501032374655803775026442524852869468234189949157306618987207994137235500057935736698933950879021244642075289741459147693018449050601793499385225470404203377985639831015709022233910000220772509651332460444439191691460859682348212832462282927101269069741823484776754573489862542033926623518620867781366509696583146995271837448054012195366666049648269890827548115254721177330319675947383719393578106059230401890711349624673706841221794681074060891827669566711716683740590473936880953450489997047176390451343232377151032196515038246988883248709353994696082647818120566349467125784366645797409778483662049777748682765697087163192938512899314199518611673792654620563505951385713761697126872299805327673278710513763"
        self.__bits = len(self.__value)*3-1   # underestimate

    def _repr_(self, simplify=True):
        """
        EXAMPLES:
            sage: repr(khinchin)
            'khinchin'
            sage: khinchin._repr_()
            'khinchin'
        """
        return 'khinchin'

    def _mpfr_(self, R):
        """
        EXAMPLES:
            sage: RealField(100)(khinchin)
            2.6854520010653064453097148355
            sage: RealField(20000)(khinchin)
            Traceback (most recent call last):
            ...
            NotImplementedError: Khinchin's constant only available up to 3005 bits
        """
        if R.precision() <= self.__bits:
            return R(self.__value)
        raise NotImplementedError, "Khinchin's constant only available up to %s bits"%self.__bits

    def _real_double_(self, R):
        """
        EXAMPLES:
            sage: RDF(khinchin)
            2.68545200107
        """
	return R('2.685452001065306445309714835481795693820')

    def _real_rqdf_(self, R):
        """
        EXAMPLES:
            sage: RQDF(khinchin)
            2.685452001065306445309714835481795693820382293994462953051152345
        """
        return R(self.__value[:65])

    def __float__(self):
        """
        EXAMPLES:
            sage: float(khinchin)
            2.6854520010653062
        """
        return float(self.__value)

    def floor(self):
        """
        Returns the floor of self.

        EXAMPLES:
            sage: khinchin.floor()
            2
        """
        return Integer(2)

khinchin  = Khinchin()


class TwinPrime(Constant):
    r"""
    The Twin Primes constant is defined as $\prod 1 - 1/(p-1)^2$
    for primes $p > 2$.

    EXAMPLES:
	sage: float(twinprime)
        0.66016181584686962
        sage: R=RealField(200);R
        Real Field with 200 bits of precision
        sage: R(twinprime)
        0.66016181584686957392781211001455577843262336028473341331945
    """
    def __init__(self):
        """
        EXAMPLES:
            sage: loads(dumps(twinprime))
            twinprime
        """
        Constant.__init__(self,{'maxima':'twinprime'}) #Twin prime is not implemented in any other algebra systems.

        #digits come from http://www.gn-50uma.de/alula/essays/Moree/Moree-details.en.shtml

        self.__value = "0.660161815846869573927812110014555778432623360284733413319448423335405642304495277143760031413839867911779005226693304002965847755123366227747165713213986968741097620630214153735434853131596097803669932135255299767199302474590593101082978291553834469297505205916657133653611991532464281301172462306379341060056466676584434063501649322723528968010934966475600478812357962789459842433655749375581854814173628678098705969498703841243363386589311969079150040573717814371081810615401233104810577794415613125444598860988997585328984038108718035525261719887112136382808782349722374224097142697441764455225265548994829771790977784043757891956590649994567062907828608828395990394287082529070521554595671723599449769037800675978761690802426600295711092099633708272559284672129858001148697941855401824639887493941711828528382365997050328725708087980662201068630474305201992394282014311102297265141514194258422242375342296879836738796224286600285358098482833679152235700192585875285961205994728621007171131607980572"

        self.__bits = len(self.__value)*3-1   # underestimate

    def floor(self):
        """
        Returns the floor of self.

        EXAMPLES:
            sage: twinprime.floor()
            0
        """
        return Integer(0)

    def _repr_(self, simplify=True):
        """
        EXAMPLES:
            sage: repr(twinprime)
            'twinprime'
            sage: twinprime._repr_()
            'twinprime'
        """
        return 'twinprime'

    def _mpfr_(self, R):
        """
        EXAMPLES:
            sage: RealField(100)(twinprime)
            0.66016181584686957392781211001
            sage: RealField(20000)(twinprime)
            Traceback (most recent call last):
            ...
            NotImplementedError: Twin Prime constant only available up to 3011 bits
        """
        if R.precision() <= self.__bits:
            return R(self.__value)
        raise NotImplementedError, "Twin Prime constant only available up to %s bits"%self.__bits

    def _real_double_(self, R):
        """
        EXAMPLES:
            sage: RDF(twinprime)
            0.660161815847
        """
	return R('0.660161815846869573927812110014555778432')

    def _real_rqdf_(self, R):
        """
        EXAMPLES:
            sage: RQDF(twinprime)
            0.660161815846869573927812110014555778432623360284733413319448422
        """
        return R(self.__value[:65])

    def __float__(self):
        """
        EXAMPLES:
            sage: float(twinprime)
            0.66016181584686962
        """
	return 0.660161815846869573927812110014555778432

twinprime = TwinPrime()


class Merten(Constant):
    """
    The Merten constant is related to the Twin Primes constant
    and appears in Merten's second theorem.

    EXAMPLES:
        sage: float(merten)
        0.26149721284764277
        sage: R=RealField(200);R
        Real Field with 200 bits of precision
        sage: R(merten)
        0.26149721284764278375542683860869585905156664826119920619206
    """
    def __init__(self):
        """
        EXAMPLES:
            sage: loads(dumps(merten))
            merten
        """
        Constant.__init__(self,{'maxima':'merten'}) #Merten's constant is not implemented in any other algebra systems.

        # digits come from Sloane's tables at http://www.research.att.com/~njas/sequences/table?a=77761&fmt=0

        self.__value = "0.261497212847642783755426838608695859051566648261199206192064213924924510897368209714142631434246651051617"

        self.__bits = len(self.__value)*3-1   # underestimate

    def floor(self):
        """
        Returns the floor of self.

        EXAMPLES:
            sage: merten.floor()
            0
        """
        return Integer(0)

    def _repr_(self, simplify=True):
        """
        EXAMPLES:
            sage: merten._repr_()
            'merten'
            sage: repr(merten)
            'merten'
        """
        return 'merten'

    def _mpfr_(self, R):
        """
        EXAMPLES:
            sage: RealField(1000)(merten)
            Traceback (most recent call last):
            ...
            NotImplementedError: Merten's constant only available up to 320 bits
            sage: RealField(320)(merten)
            0.261497212847642783755426838608695859051566648261199206192064213924924510897368209714142631434247
        """
        if R.precision() <= self.__bits:
            return R(self.__value)
        raise NotImplementedError, "Merten's constant only available up to %s bits"%self.__bits

    def _real_double_(self, R):
        """
        EXAMPLES:
            sage: RDF(merten)
            0.261497212848
        """
        return R('0.261497212847642783755426838608695859051')

    def _real_rqdf_(self, R):
        """
        EXAMPLES:
            sage: RQDF(merten)
            0.261497212847642783755426838608695859051566648261199206192064212
        """
        return R(self.__value[:65])

    def __float__(self):
        """
        EXAMPLES:
            sage: float(merten)
            0.26149721284764277
        """
	return 0.261497212847642783755426838608695859051


merten = Merten()

class Brun(Constant):
    """
    Brun's constant is the sum of reciprocals of odd twin primes.

    It is not known to very high precision; calculating the number
    using twin primes up to $10^{16}$ (Sebah 2002) gives the number
    $1.9021605831040$.

    EXAMPLES:
	sage: float(brun)
        1.902160583104
        sage: R = RealField(41); R
        Real Field with 41 bits of precision
        sage: R(brun)
        1.90216058310
    """
    def __init__(self):
        """
        EXAMPLES:
            sage: loads(dumps(brun))
            brun
        """
        Constant.__init__(self,{'maxima':"brun"}) #Brun's constant is not implemented in any other algebra systems.

        # digits come from Sloane's tables at http://www.research.att.com/~njas/sequences/table?a=65421&fmt=0

        self.__value = "1.902160583104"

        self.__bits = len(self.__value)*3-1 # bits  -- todo: make more intelligent in a general function!!!

    def floor(self):
        """
        Return the floor of self.

        EXAMPLES:
            sage: brun.floor()
            1
        """
        return Integer(1)

    def _repr_(self, simplify=True):
        """
        EXAMPLES:
            sage: brun._repr_()
            'brun'
            sage: repr(brun)
            'brun'
        """
        return 'brun'

    def _mpfr_(self, R):
        """
        EXAMPLES:
            sage: RealField(53)(brun)
            Traceback (most recent call last):
            ...
            NotImplementedError: Brun's constant only available up to 41 bits
            sage: RealField(41)(brun)
            1.90216058310
        """
        if R.precision() <= self.__bits:
            return R(self.__value)
        raise NotImplementedError, "Brun's constant only available up to %s bits"%self.__bits

    def _real_double_(self, R):
        """
        EXAMPLES:
            sage: RDF(brun)
            1.9021605831
        """
        return R('1.9021605831040')

    def _real_rqdf_(self, R):
        """
        EXAMPLES:
            sage: RealField(53)(brun)
            Traceback (most recent call last):
            ...
            NotImplementedError: Brun's constant only available up to 41 bits
        """
        raise NotImplementedError, "Brun's constant only available up to %s bits"%self.__bits

    def __float__(self):
        """
        EXAMPLES:
            sage: float(brun)
            1.902160583104
        """
	return 1.9021605831040

brun=Brun()





