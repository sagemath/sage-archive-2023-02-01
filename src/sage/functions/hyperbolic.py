"""
Hyperbolic Functions
"""
from sage.symbolic.function import GinacFunction, BuiltinFunction
import math

class HyperbolicFunction(BuiltinFunction):
    r"""
    Abstract base class for the functions defined in this file.

    EXAMPLES::

        sage: from sage.functions.hyperbolic import HyperbolicFunction
        sage: f = HyperbolicFunction('foo', latex_name='\\foo', conversions={'mathematica':'Foo'},evalf_float=lambda x: 2*x)
        sage: f(x)
        foo(x)
        sage: f(0.5r)
        1.0
        sage: latex(f(x))
        \foo\left(x\right)
        sage: f(x)._mathematica_init_()
        'Foo[x]'
    """
    def __init__(self, name, latex_name=None, conversions=None,
            evalf_float=None):
        """
        Note that subclasses of HyperbolicFunction should be instantiated only
        once, since they inherit from BuiltinFunction which only uses the
        name and class to check if a function was already registered.

        EXAMPLES::

            sage: from sage.functions.hyperbolic import HyperbolicFunction
            sage: class Barh(HyperbolicFunction):
            ...  def __init__(self):
            ...     HyperbolicFunction.__init__(self, 'barh')
            sage: barh = Barh()
            sage: barh(x)
            barh(x)
        """
        self._evalf_float = evalf_float
        BuiltinFunction.__init__(self, name, latex_name=latex_name,
                conversions=conversions)

    def _evalf_(self, x, parent):
        """
        EXAMPLES::

            sage: from sage.functions.hyperbolic import HyperbolicFunction
            sage: class Fooh(HyperbolicFunction):
            ...  def __init__(self):
            ...     HyperbolicFunction.__init__(self, 'fooh',evalf_float=lambda x: 2*x)
            sage: fooh = Fooh()
            sage: fooh(float(5))
            10.0
            sage: fooh(0.5r)
            1.0
            sage: fooh(x).subs(x=.5r)
            1.0
            sage: fooh(x).n()
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.symbolic.expression.Expression' object has no attribute 'fooh'
        """
        if parent is float:
            return self._evalf_float(x)
        return getattr(x, self.name())()

class Function_sinh(GinacFunction):
    def __init__(self):
        r"""
        The hyperbolic sine function.

        EXAMPLES::

            sage: sinh(pi)
            sinh(pi)
            sage: sinh(3.1415)
            11.5476653707437
            sage: float(sinh(pi))
            11.54873935725774...
            sage: RR(sinh(pi))
            11.5487393572577

            sage: latex(sinh(x))
            \sinh\left(x\right)
        """
        GinacFunction.__init__(self, "sinh", latex_name=r"\sinh")

sinh = Function_sinh()

class Function_cosh(GinacFunction):
    def __init__(self):
        r"""
        The hyperbolic cosine function.

        EXAMPLES::

            sage: cosh(pi)
            cosh(pi)
            sage: cosh(3.1415)
            11.5908832931176
            sage: float(cosh(pi))
            11.591953275521519
            sage: RR(cosh(1/2))
            1.12762596520638

            sage: latex(cosh(x))
            \cosh\left(x\right)
        """
        GinacFunction.__init__(self, "cosh", latex_name=r"\cosh")

cosh = Function_cosh()


class Function_tanh(GinacFunction):
    def __init__(self):
        r"""
        The hyperbolic tangent function.

        EXAMPLES::

            sage: tanh(pi)
            tanh(pi)
            sage: tanh(3.1415)
            0.996271386633702
            sage: float(tanh(pi))
            0.996272076220749...
            sage: tan(3.1415/4)
            0.999953674278156
            sage: tanh(pi/4)
            tanh(1/4*pi)
            sage: RR(tanh(1/2))
            0.462117157260010

        ::

            sage: CC(tanh(pi + I*e))
            0.997524731976164 - 0.00279068768100315*I
            sage: ComplexField(100)(tanh(pi + I*e))
            0.99752473197616361034204366446 - 0.0027906876810031453884245163923*I
            sage: CDF(tanh(pi + I*e))
            0.997524731976 - 0.002790687681*I

        TESTS::

            sage: latex(tanh(x))
            \tanh\left(x\right)
        """
        GinacFunction.__init__(self, "tanh", latex_name=r"\tanh")

tanh = Function_tanh()

class Function_coth(HyperbolicFunction):
    def __init__(self):
        r"""
        The hyperbolic cotangent function.

        EXAMPLES::

            sage: coth(pi)
            coth(pi)
            sage: coth(3.1415)
            1.00374256795520
            sage: float(coth(pi))
            1.0037418731973213
            sage: RR(coth(pi))
            1.00374187319732

            sage: latex(coth(x))
            \coth\left(x\right)
        """
        HyperbolicFunction.__init__(self, "coth", latex_name=r"\coth",
                                   evalf_float=lambda x: 1/math.tanh(x))

    def _derivative_(self, *args, **kwds):
        """
        EXAMPLES::

            sage: bool(diff(coth(x), x) == diff(1/tanh(x), x))
            True
            sage: diff(coth(x), x)
            -csch(x)^2
        """
        x = args[0]
        return -csch(x)**2

coth = Function_coth()

class Function_sech(HyperbolicFunction):
    def __init__(self):
        r"""
        The hyperbolic secant function.

        EXAMPLES::

            sage: sech(pi)
            sech(pi)
            sage: sech(3.1415)
            0.0862747018248192
            sage: float(sech(pi))
            0.0862667383340544...
            sage: RR(sech(pi))
            0.0862667383340544

            sage: latex(sech(x))
            {\rm sech}\left(x\right)
        """
        HyperbolicFunction.__init__(self, "sech", latex_name=r"{\rm sech}",
                                   evalf_float=lambda x: 1/math.cosh(x))

    def _derivative_(self, *args, **kwds):
        """
        EXAMPLES::

            sage: bool(diff(sech(x), x) == diff(1/cosh(x), x))
            True
            sage: diff(sech(x), x)
            -tanh(x)*sech(x)
        """
        x = args[0]
        return -sech(x)*tanh(x)

sech = Function_sech()


class Function_csch(HyperbolicFunction):
    def __init__(self):
        """
        The hyperbolic cosecant function.

        EXAMPLES::

            sage: csch(pi)
            csch(pi)
            sage: csch(3.1415)
            0.0865975907592133
            sage: float(csch(pi))
            0.0865895375300469...
            sage: RR(csch(pi))
            0.0865895375300470

            sage: latex(csch(x))
            {\rm csch}\left(x\right)
        """
        HyperbolicFunction.__init__(self, "csch", latex_name=r"{\rm csch}",
                                   evalf_float=lambda x: 1/math.sinh(x))

    def _derivative_(self, *args, **kwds):
        """
        EXAMPLES::

            sage: bool(diff(csch(x), x) == diff(1/sinh(x), x))
            True
            sage: diff(csch(x), x)
            -coth(x)*csch(x)
        """
        x = args[0]
        return -csch(x)*coth(x)

csch = Function_csch()


##########################
# Inverse trig functions #
##########################

class Function_arcsinh(GinacFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic sine function.

        EXAMPLES::

            sage: arcsinh
            arcsinh
            sage: arcsinh(0.5)
            0.481211825059603
            sage: arcsinh(1/2)
            arcsinh(1/2)
            sage: arcsinh(1 + I*1.0)
            1.06127506190504 + 0.666239432492515*I

        TESTS::

            sage: arcsinh(x).operator()
            arcsinh
            sage: latex(arcsinh(x))
            {\rm arcsinh}\left(x\right)
        """
        GinacFunction.__init__(self, "arcsinh", latex_name=r"{\rm arcsinh}",
                conversions=dict(maxima='asinh'))

arcsinh = asinh = Function_arcsinh()

class Function_arccosh(GinacFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic cosine function.

        EXAMPLES::

            sage: arccosh(1/2)
            arccosh(1/2)
            sage: arccosh(1 + I*1.0)
            1.06127506190504 + 0.904556894302381*I
            sage: float(arccosh(2))
            1.3169578969248168
            sage: cosh(float(arccosh(2)))
            2.0

        Warning: If the input is real the output will be real or NaN::

            sage: arccosh(0.5)
            NaN

        But evaluation where the input is in the complex field yields a
        complex output::

            sage: arccosh(CC(0.5))
            1.04719755119660*I

        TESTS::

            sage: arccosh(x).operator()
            arccosh
            sage: latex(arccosh(x))
            {\rm arccosh}\left(x\right)
        """
        GinacFunction.__init__(self, "arccosh", latex_name=r"{\rm arccosh}",
                conversions=dict(maxima='acosh'))

arccosh = acosh = Function_arccosh()

class Function_arctanh(GinacFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic tangent function.

        EXAMPLES::

            sage: arctanh(0.5)
            0.549306144334055
            sage: arctanh(1/2)
            arctanh(1/2)
            sage: arctanh(1 + I*1.0)
            0.402359478108525 + 1.01722196789785*I

        TESTS::

            sage: arctanh(x).operator()
            arctanh
            sage: latex(arctanh(x))
            {\rm arctanh}\left(x\right)
        """
        GinacFunction.__init__(self, "arctanh", latex_name=r"{\rm arctanh}",
                conversions=dict(maxima='atanh'))

arctanh = atanh = Function_arctanh()

class Function_arccoth(HyperbolicFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic cotangent function.

        EXAMPLES::

            sage: arccoth(2.0)
            0.549306144334055
            sage: arccoth(2)
            arccoth(2)
            sage: arccoth(1 + I*1.0)
            0.402359478108525 - 0.553574358897045*I
            sage: arccoth(2).n(200)
            0.54930614433405484569762261846126285232374527891137472586735
            sage: float(arccoth(2))
            0.54930614433405489

            sage: latex(arccoth(x))
            {\rm arccoth}\left(x\right)
        """
        HyperbolicFunction.__init__(self, "arccoth",
                latex_name=r"{\rm arccoth}", conversions=dict(maxima='acoth'),
                evalf_float=lambda x: atanh(float(1/x)))

    def _derivative_(self, *args, **kwds):
        """
        EXAMPLES::
            sage: bool(diff(acoth(x), x) == diff(atanh(x), x))
            True
            sage: diff(acoth(x), x)
            -1/(x^2 - 1)
        """
        x = args[0]
        return -1/(x**2 - 1)

arccoth = acoth = Function_arccoth()

class Function_arcsech(HyperbolicFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic secant function.

        EXAMPLES::

            sage: arcsech(0.5)
            1.31695789692482
            sage: arcsech(1/2)
            arcsech(1/2)
            sage: arcsech(1 + I*1.0)
            -0.530637530952518 + 1.11851787964371*I
            sage: arcsech(1/2).n(200)
            1.3169578969248167086250463473079684440269819714675164797685
            sage: float(arcsech(1/2))
            1.3169578969248168

            sage: latex(arcsech(x))
            {\rm arcsech}\left(x\right)
        """
        HyperbolicFunction.__init__(self, "arcsech",
                latex_name=r"{\rm arcsech}",
                evalf_float=lambda x: acosh(float(1/x)),
                conversions=dict(maxima='asech'))

    def _derivative_(self, *args, **kwds):
        """
        EXAMPLES::

            sage: diff(asech(x), x)
            -1/((x + 1)*x*sqrt(-(x - 1)/(x + 1)))
        """
        x = args[0]
        return -1/(x * (x+1) * ( (1-x)/(1+x) ).sqrt())

arcsech = asech = Function_arcsech()

class Function_arccsch(HyperbolicFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic cosecant function.

        EXAMPLES::

            sage: arccsch(2.0)
            0.481211825059603
            sage: arccsch(2)
            arccsch(2)
            sage: arccsch(1 + I*1.0)
            0.530637530952518 - 0.452278447151191*I
            sage: arccsch(1).n(200)
            0.88137358701954302523260932497979230902816032826163541075330
            sage: float(arccsch(1))
            0.88137358701954305

            sage: latex(arccsch(x))
            {\rm arccsch}\left(x\right)
        """
        HyperbolicFunction.__init__(self, "arccsch",
                latex_name=r"{\rm arccsch}",
                evalf_float=lambda x: arcsinh(float(1/x)),
                conversions=dict(maxima='acsch'))

    def _derivative_(self, *args, **kwds):
        """
        EXAMPLES::

            sage: diff(acsch(x), x)
            -1/(sqrt(1/x^2 + 1)*x^2)
        """
        x = args[0]
        return -1/(x**2 * (1 + x**(-2)).sqrt())

arccsch = acsch = Function_arccsch()
