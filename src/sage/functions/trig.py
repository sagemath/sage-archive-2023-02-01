"""
Trigonometric Functions
"""
from sage.symbolic.function import SFunction, PrimitiveFunction
from sage.symbolic.expression import is_Expression
import math

class Function_sin(PrimitiveFunction):
    def __init__(self):
        """
        The sine function.

        EXAMPLES::

             sage: sin(0)
             0
             sage: sin(x).subs(x==0)
             0
             sage: sin(2).n(100)
             0.90929742682568169539601986591
             sage: loads(dumps(sin))
             sin
        """
        PrimitiveFunction.__init__(self, "sin", latex=r"\sin",
                                   conversions=dict(maxima='sin',mathematica='Sin'),
                                   approx=math.sin)
sin = Function_sin()

class Function_cos(PrimitiveFunction):
    def __init__(self):
        """
        The cosine function.

        EXAMPLES:

            sage: cos(pi)
            -1
            sage: cos(x).subs(x==pi)
            -1
            sage: cos(2).n(100)
            -0.41614683654714238699756822950
            sage: loads(dumps(cos))
            cos
        """
        PrimitiveFunction.__init__(self, "cos", latex=r"\cos",
                                   conversions=dict(maxima='cos',mathematica='Cos'),
                                   approx=math.cos,)

cos = Function_cos()

class Function_tan(PrimitiveFunction):
    def __init__(self):
        """
        The tangent function

        EXAMPLES::

            sage: tan(pi)
            0
            sage: tan(3.1415)
            -0.0000926535900581913
            sage: tan(3.1415/4)
            0.999953674278156
            sage: tan(pi/4)
            1
            sage: tan(1/2)
            tan(1/2)
            sage: RR(tan(1/2))
            0.546302489843790
        """
        PrimitiveFunction.__init__(self, "tan", latex=r"\tan",
                                   approx=math.tan)
tan = Function_tan()

class Function_sec(PrimitiveFunction):
    def __init__(self):
        """
        The secant function

        EXAMPLES::

            sage: sec(pi/4)
            sqrt(2)
            sage: RR(sec(pi/4))
            1.41421356237310
            sage: n(sec(pi/4),100)
            1.4142135623730950488016887242
            sage: sec(1/2)
            sec(1/2)
            sage: sec(0.5)
            1.13949392732455

            sage: latex(sec(x))
            \sec\left(x\right)
        """
        PrimitiveFunction.__init__(self, "sec", latex=r"\sec",
                                   approx=lambda x: 1/math.cos(x))

    def _evalf_(self, x, prec=0):
        """
        EXAMPLES::

            sage: n(sec(pi/4),100)
            1.4142135623730950488016887242
        """
        return (1 / x.cos()).n(prec)

    def _eval_(self, x):
        """
        EXAMPLES::

            sage: sec(pi/4)
            sqrt(2)
            sage: sec(x).subs(x==pi/4)
            sqrt(2)
            sage: sec(pi/7)
            sec(1/7*pi)
            sage: sec(x)
            sec(x)
        """
        cos_x = cos(x)
        if is_Expression(cos_x) and cos_x.operator() is cos:
            return None
        else:
            return 1/cos_x

sec = Function_sec()

class Function_csc(PrimitiveFunction):
    def __init__(self):
        """
        The cosecant function.

        EXAMPLES::

            sage: csc(pi/4)
            sqrt(2)
            sage: RR(csc(pi/4))
            1.41421356237310
            sage: n(csc(pi/4),100)
            1.4142135623730950488016887242
            sage: csc(1/2)
            csc(1/2)
            sage: csc(0.5)
            2.08582964293349

            sage: latex(csc(x))
            \csc\left(x\right)
        """
        PrimitiveFunction.__init__(self, "csc", latex=r"\csc",
                                   approx=lambda x: 1/math.sin(x))

    def _evalf_(self, x, prec=0):
        """
        EXAMPLES::

            sage: n(csc(pi/4),100)
            1.4142135623730950488016887242
        """
        return (1 / x.sin()).n(prec)

    def _eval_(self, x):
        """
        EXAMPLES::

            sage: csc(pi/4)
            sqrt(2)
            sage: csc(x).subs(x==pi/4)
            sqrt(2)
            sage: csc(pi/7)
            csc(1/7*pi)
            sage: csc(x)
            csc(x)
        """
        sin_x = sin(x)
        if is_Expression(sin_x) and sin_x.operator() is sin:
            return None
        else:
            return 1/sin_x


csc = Function_csc()

class Function_cot(PrimitiveFunction):
    def __init__(self):
        """
        The cotangent function.

        EXAMPLES::

            sage: cot(pi/4)
            1
            sage: RR(cot(pi/4))
            1.00000000000000
            sage: cot(1/2)
            cot(1/2)
            sage: cot(0.5)
            1.83048772171245

            sage: latex(cot(x))
            \cot\left(x\right)
        """
        PrimitiveFunction.__init__(self, "cot", latex=r"\cot",
                                   approx=lambda x: 1/math.tan(x))

    def _eval_(self, x):
        """
        EXAMPLES::

            sage: cot(pi/4)
            1
            sage: cot(x).subs(x==pi/4)
            1
            sage: cot(pi/7)
            cot(1/7*pi)
            sage: cot(x)
            cot(x)
        """
        tan_x = tan(x)
        if is_Expression(tan_x) and tan_x.operator() is tan:
            return None
        else:
            return 1/tan_x

    def _evalf_(self, x, prec=0):
        """
        EXAMPLES::

            sage: n(cot(pi/4),100)
            1.0000000000000000000000000000
        """
        return x.n(prec).cot()

cot = Function_cot()


###################################
# Inverse Trigonometric Functions #
###################################

class Function_arcsin(PrimitiveFunction):
    def __init__(self):
        """
        The arcsine function.

        EXAMPLES::

            sage: arcsin(0.5)
            0.523598775598299
            sage: arcsin(1/2)
            1/6*pi
            sage: arcsin(1 + 1.0*I)
            0.666239432492515 + 1.06127506190504*I

        TESTS::

            sage: arcsin(x).operator()
            arcsin
        """
        PrimitiveFunction.__init__(self, 'arcsin', latex=r"\sin^{-1}",
                approx=math.asin,
                conversions=dict(maxima='asin', ginac='asin'))
arcsin = asin = Function_arcsin()

class Function_arccos(PrimitiveFunction):
    def __init__(self):
        """
        The arccosine function

        EXAMPLES::

            sage: arccos(0.5)
            1.04719755119660
            sage: arccos(1/2)
            1/3*pi
            sage: arccos(1 + 1.0*I)
            0.904556894302381 - 1.06127506190504*I

        TESTS::

            sage: arccos(x).operator()
            arccos
        """
        PrimitiveFunction.__init__(self, 'arccos', latex=r"\cos^{-1}",
                approx=math.acos,
                conversions=dict(maxima='acos', ginac='acos'))

    def _evalf_(self, x, prec=0):
        """
        EXAMPLES::

            sage: arccos(3/4).n(100)
            0.72273424781341561117837735264
        """
        return x.n(prec).arccos()

arccos = acos = Function_arccos()

class Function_arctan(PrimitiveFunction):
    def __init__(self):
        """
        The arctangent function.

        EXAMPLES::

            sage: arctan(1/2)
            arctan(1/2)
            sage: RDF(arctan(1/2))
            0.463647609001
            sage: arctan(1 + I)
            arctan(I + 1)

        TESTS::

            sage: arctan(x).operator()
            arctan
        """
        PrimitiveFunction.__init__(self, "arctan", latex=r'\tan^{-1}',
                                   approx=math.atan,
                                   conversions=dict(maxima='atan',
                                                    ginac='atan'))

    def _evalf_(self, x, prec=0):
        """
        EXAMPLES::

            sage: arctan(1/2).n(100)
            0.46364760900080611621425623146
        """
        return x.n(prec).arctan()

arctan = atan = Function_arctan()

class Function_arccot(PrimitiveFunction):
    def __init__(self):
        """
        The arccotangent function.

        EXAMPLES::

            sage: arccot(1/2)
            arccot(1/2)
            sage: RDF(arccot(1/2))
            1.10714871779
            sage: arccot(1 + I)
            arccot(I + 1)
        """
        PrimitiveFunction.__init__(self, "arccot", latex=r'\cot^{-1}',
                                   approx=lambda x: math.pi/2 - math.atan(x),
                                   conversions=dict(maxima='acot'))

    def _evalf_(self, x, prec=0):
        """
        EXAMPLES::

            sage: arccot(1/2).n(100)
            1.1071487177940905030170654602
        """
        from sage.symbolic.constants import pi
        return (pi/2 - x.arctan()).n(prec)

arccot = acot = Function_arccot()

class Function_arccsc(PrimitiveFunction):
    def __init__(self):
        """
        The arccosecant function.

        EXAMPLES::

            sage: arccsc(2)
            arccsc(2)
            sage: RDF(arccsc(2))
            0.523598775598
            sage: arccsc(1 + I)
            arccsc(I + 1)
        """
        PrimitiveFunction.__init__(self, "arccsc", latex=r'\csc^{-1}',
                                   approx=lambda x: math.asin(1/x),
                                   conversions=dict(maxima='acsc'))

    def _evalf_(self, x, prec=0):
        """
        EXAMPLES::

            sage: arccsc(2).n(100)
            0.52359877559829887307710723055
        """
        return (1/x).arcsin().n(prec)

arccsc = acsc = Function_arccsc()

class Function_arcsec(PrimitiveFunction):
    def __init__(self):
        """
        The arcsecant function.

        EXAMPLES::

            sage: arcsec(2)
            arcsec(2)
            sage: RDF(arcsec(2))
            1.0471975512
            sage: arcsec(1 + I)
            arcsec(I + 1)
        """
        PrimitiveFunction.__init__(self, "arcsec", latex=r'\sec^{-1}',
                                   approx=lambda x: math.acos(1/x),
                                   conversions=dict(maxima='asec'))

    def _evalf_(self, x, prec=0):
        """
        EXAMPLES::

            sage: arcsec(2).n(100)
            1.0471975511965977461542144611
        """
        return (1/x).arccos().n(prec)

arcsec = asec = Function_arcsec()

class Function_arctan2(PrimitiveFunction):
    def __init__(self):
        """
        The modified arctangent function.

        Returns the arc tangent (measured in radians) of `y/x`, where
        unlike ``arctan(y/x)``, the signs of both ``x`` and ``y`` are
        considered.

        Note that the `y`-coordinate is by convention the first input.

        ``arctan2(y,x) = arctan(y/x)``

        EXAMPLES:

        Note the difference between the two functions::

            sage: arctan2(1,-1)
            3/4*pi
            sage: arctan(1/-1)
            -1/4*pi

        This is consistent with Python and Maxima::

            sage: maxima.atan2(1,-1)
            3*%pi/4
            sage: math.atan2(1,-1)
            2.3561944901923448

        More examples::

            sage: arctan2(1,0)
            1/2*pi
            sage: arctan2(2,3)
            arctan(2/3)
            sage: arctan2(-1,-1)
            -3/4*pi

        Of course we can approximate as well::

            sage: arctan2(-.5,1).n(100)
            -0.46364760900080611621425623146

        TESTS::

            sage: x,y = var('x,y')
            sage: arctan2(y,x).operator()
            arctan2
        """
        PrimitiveFunction.__init__(self, "arctan2", nargs=2, latex=r'\arctan',
                                   approx=math.atan2,
                                   conversions=dict(maxima='atan2',
                                                    ginac='atan2'))

    __call__ = SFunction.__call__

arctan2 = atan2 = Function_arctan2()
