from sage.rings.all import CommutativeRing, RealField, is_Polynomial
import sage.rings.all
from sage.structure.all import RingElement, Element_cmp_
import operator
from sage.misc.latex import latex
import constants
from sage.interfaces.all import maxima

RR = RealField(53)

class FunctionRing_class(CommutativeRing):
    def __init__(self):
        self._default_precision = 53  # default bits of precision

    def _repr_(self):
        return "Ring of Mathematical Functions"

    def _latex_(self):
        return '\\text{RingOfFunctions}'

    def set_precision(self, prec):
        """
        Change the precision of the default real field used for coercing functions.

        EXAMPLES:
            sage: old_prec = FunctionRing.set_precision(200)
            sage: pi.str()
            '3.1415926535897932384626433832795028841971693993751058209749445'
            sage: gap(pi)
            "3.1415926535897932384626433832795028841971693993751058209749445"
            sage: _ = FunctionRing.set_precision(old_prec)
            sage: pi.str()
            '3.1415926535897931'

        Note that there seems to be no real numbers in GAP, which is why that
        coercion returns a GAP string.

        The precision is only used when there is no way to describe
        the function to higher precision (or exactly) to the
        interface:
             sage: maxima(pi)
             %pi

        If we coerce into a GP/PARI session, then the resulting number will
        have precision the precision of that session irregardless of the
        default precision of FunctionRing.  Note that GP/PARI precision is set in
        decimal digits, whereas in SAGE precisions are always in binary.

            sage: gp = Gp()  # new session so as not to affect overwrite global gp session
            sage: gp(pi)
            3.141592653589793238462643383             # 32-bit
            3.1415926535897932384626433832795028842   # 64-bit
            sage: gp.eval('\\p 100')
            '   realprecision = 105 significant digits (100 digits displayed)'      # 32-bit
            '   realprecision = 115 significant digits (100 digits displayed)'      # 64-bit
            sage: gp(pi)
            3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
            sage: _ = FunctionRing.set_precision(5)   # has no affect on gp(pi)
            sage: gp(pi)
            3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
        """
        old_prec = self._default_precision
        prec = int(prec)
        self._default_precision = prec
        global RR
        RR = sage.rings.all.RealField(prec)
        return old_prec

    def default_precision(self):
        """
        Get the precision of the default real field used for coercing functions.
        """
        return self._default_precision

    def __call__(self, x):
        try:
            return self._coerce_(x)
        except TypeError:
            return Function_gen(x)

    def _coerce_(self, x):
        if isinstance(x, Function):
            return x
        elif is_Polynomial(x):
            return Function_polynomial(x)
        #elif isinstance(x, sage.ext.element.Element):
        elif isinstance(x, (sage.ext.integer.Integer,
                            sage.ext.rational.Rational)):
            return constants.Constant_gen(x)
        raise TypeError

    def characteristic(self):
        return sage.rings.all.Integer(0)

FunctionRing = FunctionRing_class()

class Function(Element_cmp_, RingElement):
    def __init__(self, conversions={}):
        self.__conversions = conversions
        RingElement.__init__(self, FunctionRing)

    def __call__(self, x):
        if isinstance(x, Function) and not isinstance(x, constants.Constant):
            return Function_composition(self, x)

        elif is_Polynomial(x):
            return Function_composition(self, Function_polynomial(x))

        try:
            return self._call_(x)
        except AttributeError:
            pass
        return Function_at(self, x)

    def str(self, bits=None):
        """
        Return a string representation of self as a decimal number
        to at least the given bits of precision.
        """
        if not bits is None:
            R = sage.rings.all.RealField(bits)
        else:
            R = RR
        return str(self._mpfr_(R))

    def plot(self, n=50):
        n = int(n)
        import pylab
        pylab.plot([float(self(i/n)) for i in range(n)])

    def integral(self):
        raise NotImplementedError, "computation of integral of %s not implemented"%self

    def _interface_init_(self):
        """
        Default init string for interfaces.
        """
        return self.str()

    def _gap_init_(self):
        try:
            return self.__conversions['gap']
        except KeyError:
            return '"%s"'%self.str()  # no float numbers in GAP!?!

    def _kash_init_(self):
        try:
            return self.__conversions['kash']
        except KeyError:
            return self.str()

    def _maxima_init_(self):
        try:
            return self.__conversions['maxima']
        except KeyError:
            return self.str()

    def _mathematica_init_(self):
        try:
            return self.__conversions['mathematica']
        except KeyError:
            return self.str()

    def _maple_init_(self):
        try:
            return self.__conversions['maple']
        except KeyError:
            return self.str()

    def _octave_init_(self):
        try:
            return self.__conversions['octave']
        except KeyError:
            return self.str()

    def _pari_init_(self):
        try:
            return self.__conversions['pari']
        except KeyError:
            return self.str()

    def _singular_init_(self):
        """
        EXAMPLES:
            sage: singular(e)
            2.7182818284590451
            sage: P = e.parent()
            sage: old_prec = P.set_precision(200)
            sage: singular(e)
            2.7182818284590452353602874713526624977572470936999595749669679
            sage: _ = P.set_precision(old_prec)
        """
        try:
            return self.__conversions['singular']
        except KeyError:
            return '"%s"'%self.str()

    def _mpfr_(self, R):
        raise NotImplementedError, "computation of this function not implemented yet."

    def __float__(self):
        return float(self._mpfr_(RR))

    # The following adds formal arithmetic support for functions
    def _add_(self, right):
        return Function_arith(self, right, operator.add)

    def _sub_(self, right):
        return Function_arith(self, right, operator.sub)

    def _mul_(self, right):
        return Function_arith(self, right, operator.mul)

    def _div_(self, right):
        return Function_arith(self, right, operator.div)

    def __pow__(self, right):
        try:
            right = self.parent()._coerce_(right)
        except TypeError:
            raise TypeError, "computation of %s^%s not defined"%(self, right)
        return Function_arith(self, right, operator.pow)

    def _cmp_(self, right):
        """
        EXAMPLES:
            sage: s = e + pi
            sage: s == 0
            False
            sage: t = e^2  +pi
            sage: s == t
            False
            sage: s == s
            True
            sage: t == t
            True
            sage: s < t
            True
            sage: t < s
            False
        """
        if self is right:
            return 0
        R = RealField()
        c = cmp(R(self), R(right))
        if c: return c
        try:
	    return cmp(maxima(self),maxima(right))
	except TypeError:
	    pass
        raise NotImplementedError, "equality testing for general mathematical functions not yet implemented (can't prove equality)."


class Function_composition(Function):
    def __init__(self, f, g):
        self.__f = f
        self.__g = g
        Function.__init__(self)

    def _repr_(self):
        return '%s(%s)'%(self.__f, self.__g)

    def _latex_(self):
        return '%s(%s)'%(latex(self.__f), latex(self.__g))

    def _call_(self, x):
        return self.__f(self.__g(x))


#################################################################
#
# Support for arithmetic with constants.
#
#################################################################

symbols = {operator.add:' + ', operator.sub:' - ',
           operator.mul:'*', operator.div:'/',
           operator.pow:'^'}

class Function_arith(Function):
    """
    Generic arithmetic with functions is supported.

    EXAMPLES:
        sage: s = (pi + pi) * e + e
        sage: s
        (((pi + pi)*e) + e)
        sage: R(s)
        19.797750273806177
        sage: maxima(s)
        2*%e*%pi + %e

        sage: t = e^2 + pi + 2/3; t
        (((e^2) + pi) + 2/3)
        sage: R(t)
        11.197315419187108
        sage: maxima(t)
        %pi + %e^2 + 2/3
        sage: t^e
        ((((e^2) + pi) + 2/3)^e)
        sage: R(t^e)
        710.86524768885772
    """
    def __init__(self, x, y, op):
        if not isinstance(x, Function) or not isinstance(y, Function):
            raise TypeError
        Function.__init__(self, {})
        self.__x = x
        self.__y = y
        self.__op = op  # a binary operator, e.g., operator.sub

    def _repr_(self):
        """
        EXAMPLES:
            sage: log2 * e + pi^2/2
            ((log2*e) + ((pi^2)/2))
        """
        return '(%s%s%s)'%(self.__x, symbols[self.__op], self.__y)

    def _latex_(self):
        """
        EXAMPLES:
            sage: latex(log2 * e + pi^2/2)
            \log(2) \cdot e + \frac{\pi^{2}}{2}
            sage: latex(NaN^3 + 1/golden_ratio)
            \mbox{\rm NaN}^{3} + \frac{1}{\phi}
            sage: latex(log2 + euler_gamma + catalan + khinchin + twinprime + merten + brun)
            \log(2) + \gamma + K + \mbox{\rm khinchin} + \mbox{\rm twinprime} + \mbox{\rm merten} + \mbox{\rm brun}
        """
        if self.__op == operator.div:
            return '\\frac{%s}{%s}'%(latex(self.__x), latex(self.__y))
        elif self.__op == operator.mul:
            return '%s \\cdot %s'%(latex(self.__x), latex(self.__y))
        elif self.__op == operator.sub:
            return '%s - %s'%(latex(self.__x), latex(self.__y))
        elif self.__op == operator.add:
            return '%s + %s'%(latex(self.__x), latex(self.__y))
        elif self.__op == operator.pow:
            return '%s^{%s}'%(latex(self.__x), latex(self.__y))
        else:
            raise NotImplementedError, 'operator %s unknown'%self.__op

    def _call_(self, x):
        if isinstance(self, constants.Constant):
            return self
        return self.__op(self.__x(x), self.__y(x))

    def _gap_init_(self):
        """
        EXAMPLES:
            sage: gap(e + pi)
            "5.8598744820488378"
        """
        return '"%s"'%self.str()

    def _gp_(self, gp):
        """
        EXAMPLES:
            sage: gp(e + pi)
            5.859874482048838473822930855             # 32-bit
            5.8598744820488384738229308546321653819   # 64-bit
        """
        return self.__op(self.__x._gp_(gp), self.__y._gp_(gp))

    def _kash_(self, kash):
        """
        EXAMPLES:
            sage: kash(e + pi)                        # optional
            5.85987448204883847382293085463
        """
        return self.__op(self.__x._kash_(kash), self.__y._kash_(kash))

    def _maple_(self, maple):
        """
        EXAMPLES:
            sage: maple(e + pi)              # optional
            exp(1)+Pi
        """
        return self.__op(self.__x._maple_(maple), self.__y._maple_(maple))

    def _mathematica_(self, mathematica):
        """
        EXAMPLES:
            sage: mathematica(e + pi)        # optional
            E + Pi
        """
        return self.__op(self.__x._mathematica_(mathematica), self.__y._mathematica_(mathematica))

    def _maxima_(self, maxima):
        """
        EXAMPLES:
            sage: maxima(e + pi)
            %pi + %e
        """
        return self.__op(self.__x._maxima_(maxima), self.__y._maxima_(maxima))

    def _octave_(self, octave):
        """
        EXAMPLES:
            sage: octave(e + pi)                     # optional
            5.85987
        """
        return self.__op(self.__x._octave_(octave), self.__y._octave_(octave))

    def _pari_(self):
        """
        EXAMPLES:
            sage: pari(e + pi)
            5.859874482048838473822930855             # 32-bit
            5.8598744820488384738229308546321653819   # 64-bit
        """
        return self.__op(self.__x._pari_(), self.__y._pari_())

    def _singular_init_(self):
        """
        EXAMPLES:
            sage: singular(e + pi)
            5.8598744820488378
        """
        return '"%s"'%self.str()

    def _mpfr_(self, R):
        """
        EXAMPLES:
            sage: RealField(100)(e + pi)
            5.8598744820488384738229308546306
        """
        return self.__op(self.__x._mpfr_(R), self.__y._mpfr_(R))


class Function_gen(Function):
    """
    Function defined by a generic SAGE object.  This makes possible
    symbolic expressions like the following:

    EXAMPLES:
        sage: a = pi/2 + e
        sage: a
        ((pi/2) + e)
        sage: maxima(a)
        %pi/2 + %e
        sage: RR(a)
        4.2890781552539412
        sage: RealField(200)(a)
        4.2890781552539418545916091629924139398558317933875124854544427

        sage: b = e + 5/7
        sage: maxima(b)
        %e + 5/7
        sage: RR(b)
        3.4325675427447595
    """
    def __init__(self, x):
        Function.__init__(self)
        self.__x = x

    def obj(self):
        return self.__x

    def _repr_(self):
        return str(self.__x)

    def _latex_(self):
        return latex(self.__x)

    def _mpfr_(self, R):
        return R(self.__x)

    def str(self, bits=None):
        if bits is None:
            return str(self.__x)
        else:
            R = sage.rings.all.RealField(53)
            return str(R(self))

    def _gap_(self, gap):
        return gap(self.__x)

    def _gp_(self, gp):
        return gp(self.__x)

    def _kash_(self, kash):
        return kash(self.__x)

    def _maxima_(self, maxima):
        return maxima(self.__x)

    def _mathematica_(self, mathematica):
        return mathematica(self.__x)

    def _maple_(self, maple):
        return maple(self.__x)

    def _octave_(self, octave):
        return octave(self.__x)

    def _pari_(self):
        return pari(self.__x)

    def _singular_(self, singular):
        return singular(self.__x)


class Function_polynomial(Function):
    def __init__(self, f):
        self.__f = f
        Function.__init__(self)

    def _repr_(self):
        return str(self.__f)

    def _call_(self, x):
        return self.__f(x)

    def _latex_(self):
        return latex(self.__f)


class Function_at(Function):
    def __init__(self, f, x):
        self.__f = f
        self.__x = x
        Function.__init__(self)

    def _repr_(self):
        return '%s(%s)'%(self.__f, self.__x)

    def _mpfr_(self, R):
        x = R(self.__x)
        return self.__f(x)

    def _maxima_init_(self):
        try:
            return '%s(%s)'%(self.__f._maxima_init_(), self.__x._maxima_init_())
        except AttributeError:
            raise NotImplementedError, 'coercion of %s to maxima not implemented'%self



########################################################

class Function_sin(Function):
    def __init__(self):
        Function.__init__(self,
            {'maxima':'sin', 'mathematica':'Sin'})

    def _repr_(self):
        return "sin"

    def _latex_(self):
        return "\\sin"

    def _mpfr_(self, R):
        raise TypeError

    def _call_(self, x):
        return x.sin()

    def integral(self):
        return -cos

sin = Function_sin()

class Function_cos(Function):
    def __init__(self):
        Function.__init__(self,
            {'maxima':'cos', 'mathematica':'Cos'})

    def _repr_(self):
        return "cos"

    def _latex_(self):
        return "\\cos"

    def _mpfr_(self, R):
        raise TypeError

    def _call_(self, x):
        return x.cos()

    def integral(self):
        return -sin

cos = Function_cos()


class Function_maxima(Function):
    def __init__(self, x):
        Function.__init__(self, {'maxima':str(x)})
        self.__x = x
        self.__name = x.name()
        self.__p = x.parent()

    def _repr_(self):
        return str(self.__x)

    def _latex_(self):
        return x._latex_()

    def _mpfr_(self, R):
        raise TypeError

    def _call_(self, x):
        return float(self.__p.eval('%s(%s)'%(self.__name, x)))
        #return self.__x(x)

    def integral(self):
        return Function_maxima(self.__x.integrate())

#def maxima_function(s):
#    return Function_maxima(maxima(s))

