"""
Hyperbolic Functions
"""
from sage.symbolic.function import SFunction, PrimitiveFunction
from sage.libs.pari.gen import pari
import math

class HyperbolicFunction(PrimitiveFunction):
    def _evalf_(self, x, prec=0):
        """
        EXAMPLES::

            sage: arccosh._evalf_(2, 53)
            1.31695789692482
        """
        return getattr(x.n(prec), self.name())()

class Function_sinh(HyperbolicFunction):
    def __init__(self):
        """
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
        """
        PrimitiveFunction.__init__(self, "sinh", latex=r"\sinh",
                                   approx=math.sinh)

sinh = Function_sinh()

class Function_cosh(HyperbolicFunction):
    def __init__(self):
        """
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
        """
        PrimitiveFunction.__init__(self, "cosh", latex=r"\cosh",
                                   approx=math.cosh)

cosh = Function_cosh()


class Function_tanh(HyperbolicFunction):
    def __init__(self):
        """
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
        """
        PrimitiveFunction.__init__(self, "tanh", latex=r"\tanh",
                                   approx=math.tanh)

tanh = Function_tanh()

class Function_coth(HyperbolicFunction):
    def __init__(self):
        """
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
        """
        PrimitiveFunction.__init__(self, "coth", latex=r"\coth",
                                   approx=lambda x: 1/math.tanh(x))

coth = Function_coth()

class Function_sech(HyperbolicFunction):
    def __init__(self):
        """
        The hyperbolic secant function.

        EXAMPLES::

            sage: sech(pi)
            sech(pi)
            sage: sech(3.1415)
            0.0862747018248192
            sage: float(sech(pi))
            0.086266738334054419
            sage: RR(sech(pi))
            0.0862667383340544
        """
        PrimitiveFunction.__init__(self, "sech", latex=r"\sech",
                                   approx=lambda x: 1/math.cosh(x))

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
        """
        PrimitiveFunction.__init__(self, "csch", latex=r"\text{csch}",
                                   approx=lambda x: 1/math.sinh(x))

csch = Function_csch()


##########################
# Inverse trig functions #
##########################

class Function_arcsinh(HyperbolicFunction):
    def __init__(self):
        """
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
            sage: arcsinh._approx_(0.5)
            0.48121182505960347

        """
        PrimitiveFunction.__init__(self, "arcsinh", latex=r"\sinh^{-1}",
                                   approx=lambda x: float(pari(float(x)).asinh()),
                                   conversions=dict(maxima='asinh'))

arcsinh = asinh = Function_arcsinh()

class Function_arccosh(HyperbolicFunction):
    def __init__(self):
        """
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

        """
        PrimitiveFunction.__init__(self, "arccosh", latex=r"\cosh^{-1}",
                                   approx=lambda x: float(pari(float(x)).acosh()),
                                   conversions=dict(maxima='acosh'))

arccosh = acosh = Function_arccosh()

class Function_arctanh(HyperbolicFunction):
    def __init__(self):
        """
        The inverse of the hyperbolic tangent function.

        EXAMPLES::

            sage: arctanh(0.5)
            0.549306144334055
            sage: arctanh(1/2)
            arctanh(1/2)
            sage: arctanh(1 + I*1.0)
            0.402359478108525 + 1.01722196789785*I

        """
        PrimitiveFunction.__init__(self, "arctanh", latex=r"\tanh^{-1}",
                                   approx=lambda x: float(pari(float(x)).atanh()),
                                   conversions=dict(maxima='atanh'))

arctanh = atanh = Function_arctanh()

class Function_arccoth(HyperbolicFunction):
    def __init__(self):
        """
        The inverse of the hyperbolic cotangent function.

        EXAMPLES::

            sage: arccoth(2.0)
            0.549306144334055
            sage: arccoth(2)
            arccoth(2)
            sage: arccoth(1 + I*1.0)
            0.402359478108525 - 0.553574358897045*I
        """
        PrimitiveFunction.__init__(self, "arccoth", latex=r"\coth^{-1}",
                                   approx=lambda x: float(pari(float(1/x)).atanh()),
                                   conversions=dict(maxima='acoth'))

    def _evalf_(self, x, prec=0):
        """
        EXAMPLES::

            sage: arccoth(2).n(200)
            0.54930614433405484569762261846126285232374527891137472586735
        """
        return arctanh(1/x).n(prec)

arccoth = acoth = Function_arccoth()

class Function_arcsech(HyperbolicFunction):
    def __init__(self):
        """
        The inverse of the hyperbolic secant function.

        EXAMPLES::

            sage: arcsech(0.5)
            1.31695789692482
            sage: arcsech(1/2)
            arcsech(1/2)
            sage: arcsech(1 + I*1.0)
            -0.530637530952518 + 1.11851787964371*I
        """
        PrimitiveFunction.__init__(self, "arcsech", latex=r"\text{sech}^{-1}",
                                   approx=lambda x: float(pari(float(1/x)).acosh()),
                                   conversions=dict(maxima='asech'))

    def _evalf_(self, x, prec=0):
        """
        EXAMPLES::

            sage: arcsech(1/2).n(200)
            1.3169578969248167086250463473079684440269819714675164797685
        """
        return arccosh(1/x).n(prec)

arcsech = asech = Function_arcsech()

class Function_arccsch(HyperbolicFunction):
    def __init__(self):
        """
        The inverse of the hyperbolic cosecant function.

        EXAMPLES::

            sage: arccsch(2.0)
            0.481211825059603
            sage: arccsch(2)
            arccsch(2)
            sage: arccsch(1 + I*1.0)
            0.530637530952518 - 0.452278447151191*I
        """
        PrimitiveFunction.__init__(self, "arccsch", latex=r"\text{csch}^{-1}",
                                   approx=lambda x: float(pari(float(1/x)).arcsinh()),
                                   conversions=dict(maxima='acsch'))

    def _evalf_(self, x, prec=0):
        """
        EXAMPLES::

            sage: arccsch(1).n(200)
            0.88137358701954302523260932497979230902816032826163541075330
        """
        return arcsinh(1/x).n(prec)

arccsch = acsch = Function_arccsch()
