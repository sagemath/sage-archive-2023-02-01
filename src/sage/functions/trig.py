"""
Trigonometric Functions
"""
from sage.symbolic.function import BuiltinFunction, GinacFunction
from sage.symbolic.expression import is_Expression
import math

class Function_sin(GinacFunction):
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
        GinacFunction.__init__(self, "sin", latex_name=r"\sin",
                conversions=dict(maxima='sin',mathematica='Sin'))

sin = Function_sin()

class Function_cos(GinacFunction):
    def __init__(self):
        """
        The cosine function.

        EXAMPLES::

            sage: cos(pi)
            -1
            sage: cos(x).subs(x==pi)
            -1
            sage: cos(2).n(100)
            -0.41614683654714238699756822950
            sage: loads(dumps(cos))
            cos
        """
        GinacFunction.__init__(self, "cos", latex_name=r"\cos",
                conversions=dict(maxima='cos',mathematica='Cos'))

cos = Function_cos()

class Function_tan(GinacFunction):
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
        GinacFunction.__init__(self, "tan", latex_name=r"\tan")

tan = Function_tan()

class Function_sec(BuiltinFunction):
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
        BuiltinFunction.__init__(self, "sec", latex_name=r"\sec")

    def _evalf_(self, x, parent=None):
        """
        EXAMPLES::

            sage: n(sec(pi/4),100)
            1.4142135623730950488016887242
            sage: float(sec(pi/4))
            1.4142135623730951
        """
        if parent is float:
            return 1/math.cos(x)
        return (1 / x.cos())

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

    def _derivative_(self, x, diff_param=None):
        """
        EXAMPLES::

            sage: bool(diff(sec(x), x) == diff(1/cos(x), x))
            True
            sage: diff(sec(x), x)
            tan(x)*sec(x)
        """
        return sec(x)*tan(x)

sec = Function_sec()

class Function_csc(BuiltinFunction):
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
        BuiltinFunction.__init__(self, "csc", latex_name=r"\csc")

    def _evalf_(self, x, parent=None):
        """
        EXAMPLES::

            sage: n(csc(pi/4),100)
            1.4142135623730950488016887242
            sage: float(csc(pi/4))
            1.4142135623730951
        """
        if parent is float:
            return 1/math.sin(x)
        return (1 / x.sin())

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

    def _derivative_(self, x, diff_param=None):
        """
        EXAMPLES::

            sage: bool(diff(csc(x), x) == diff(1/sin(x), x))
            True
            sage: diff(csc(x), x)
            -csc(x)*cot(x)
        """
        return -csc(x)*cot(x)

csc = Function_csc()

class Function_cot(BuiltinFunction):
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
        BuiltinFunction.__init__(self, "cot", latex_name=r"\cot")

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

    def _evalf_(self, x, parent=None):
        """
        EXAMPLES::

            sage: n(cot(pi/4),100)
            1.0000000000000000000000000000
            sage: float(cot(1))
            0.64209261593433...
        """
        if parent is float:
            return 1/math.tan(x)
        return x.cot()

    def _derivative_(self, x, diff_param=None):
        """
        EXAMPLES::

            sage: bool(diff(cot(x), x) == diff(1/tan(x), x))
            True
            sage: diff(cot(x), x)
            -csc(x)^2
        """
        return -csc(x)**2

cot = Function_cot()


###################################
# Inverse Trigonometric Functions #
###################################

class Function_arcsin(GinacFunction):
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
        GinacFunction.__init__(self, 'arcsin', latex_name=r"\arcsin",
                conversions=dict(maxima='asin'))

arcsin = asin = Function_arcsin()

class Function_arccos(GinacFunction):
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
            sage: arccos(3/4).n(100)
            0.72273424781341561117837735264

        TESTS::

            sage: arccos(x).operator()
            arccos
        """
        GinacFunction.__init__(self, 'arccos', latex_name=r"\arccos",
                conversions=dict(maxima='acos'))

arccos = acos = Function_arccos()

class Function_arctan(GinacFunction):
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
            sage: arctan(1/2).n(100)
            0.46364760900080611621425623146

        TESTS::

            sage: arctan(x).operator()
            arctan
        """
        GinacFunction.__init__(self, "arctan", latex_name=r'\arctan',
                conversions=dict(maxima='atan'))

arctan = atan = Function_arctan()

class Function_arccot(BuiltinFunction):
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
        BuiltinFunction.__init__(self, "arccot", latex_name=r'{\rm arccot}',
                                   conversions=dict(maxima='acot'))

    def _evalf_(self, x, parent=None):
        """
        EXAMPLES::

            sage: arccot(1/2).n(100)
            1.1071487177940905030170654602
            sage: float(arccot(1/2))
            1.1071487177940904
        """
        if parent is float:
            return math.pi/2 - math.atan(x)
        from sage.symbolic.constants import pi
        return parent(pi/2 - x.arctan())

    def _derivative_(self, *args, **kwds):
        """
        EXAMPLES::

            sage: bool(diff(acot(x), x) == -diff(atan(x), x))
            True
            sage: diff(acot(x), x)
            -1/(x^2 + 1)
        """
        x = args[0]
        return -1/(x**2 + 1)

arccot = acot = Function_arccot()

class Function_arccsc(BuiltinFunction):
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
        BuiltinFunction.__init__(self, "arccsc", latex_name=r'{\rm arccsc}',
                                   conversions=dict(maxima='acsc'))

    def _evalf_(self, x, parent=None):
        """
        EXAMPLES::

            sage: arccsc(2).n(100)
            0.52359877559829887307710723055
            sage: float(arccsc(2))
            0.52359877559829...
        """
        if parent is float:
            return math.asin(1/x)
        return (1/x).arcsin()

    def _derivative_(self, x, diff_param=None):
        """
        EXAMPLES::

            sage: diff(acsc(x), x)
            -1/(sqrt(-1/x^2 + 1)*x^2)
        """
        return -1/(x**2 * (1 - x**(-2)).sqrt())

arccsc = acsc = Function_arccsc()

class Function_arcsec(BuiltinFunction):
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
        BuiltinFunction.__init__(self, "arcsec", latex_name=r'{\rm arcsec}',
                                   conversions=dict(maxima='asec'))

    def _evalf_(self, x, parent=None):
        """
        EXAMPLES::

            sage: arcsec(2).n(100)
            1.0471975511965977461542144611
        """
        if parent is float:
            return math.acos(1/x)
        return (1/x).arccos()

    def _derivative_(self, x, diff_param=None):
        """
        EXAMPLES::

            sage: diff(asec(x), x)
            1/(sqrt(-1/x^2 + 1)*x^2)
        """
        return 1/(x**2 * (1 - x**(-2)).sqrt())

arcsec = asec = Function_arcsec()

class Function_arctan2(GinacFunction):
    def __init__(self):
        """
        The modified arctangent function.

        Returns the arc tangent (measured in radians) of `y/x`, where
        unlike ``arctan(y/x)``, the signs of both ``x`` and ``y`` are
        considered.  In particular, this function measures the angle
        of a ray through the origin and `(x,y)`, with the positive
        `x`-axis the zero mark, and with output angle `\theta`
        being between `-\pi<\theta<=\pi`.

        Hence, ``arctan2(y,x) = arctan(y/x)`` only for `x>0`.  One
        may consider the usual arctan to measure angles of lines
        through the origin, while the modified function measures
        rays through the origin.

        Note that the `y`-coordinate is by convention the first input.


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

            sage: arctan2(-1/2,1).n(100)
            -0.46364760900080611621425623146
            sage: arctan2(2,3).n(100)
            0.58800260354756755124561108063

        TESTS::

            sage: x,y = var('x,y')
            sage: arctan2(y,x).operator()
            arctan2

        Check if #8565 is fixed::

            sage: atan2(-pi,0)
            -1/2*pi
        """
        GinacFunction.__init__(self, "arctan2", nargs=2, latex_name=r'\arctan',
                conversions=dict(maxima='atan2'))

arctan2 = atan2 = Function_arctan2()
