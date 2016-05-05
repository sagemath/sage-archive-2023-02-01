r"""
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

        We can prevent evaluation using the ``hold`` parameter::

            sage: sin(0,hold=True)
            sin(0)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = sin(0,hold=True); a.simplify()
            0

        TESTS::

            sage: conjugate(sin(x))
            sin(conjugate(x))
            sage: sin(complex(1,1))     # rel tol 1e-15
            (1.2984575814159773+0.6349639147847361j)

            sage: sin(pi/5)
            1/4*sqrt(-2*sqrt(5) + 10)
            sage: sin(pi/8)
            1/2*sqrt(-sqrt(2) + 2)
            sage: sin(pi/24)
            1/4*sqrt(-2*sqrt(6) - 2*sqrt(2) + 8)
            sage: sin(pi/30)
            -1/8*sqrt(5) + 1/4*sqrt(-3/2*sqrt(5) + 15/2) - 1/8
            sage: cos(pi/8)
            1/2*sqrt(sqrt(2) + 2)
            sage: cos(pi/10)
            1/2*sqrt(1/2*sqrt(5) + 5/2)
            sage: cos(pi/12)
            1/12*sqrt(6)*(sqrt(3) + 3)
            sage: cos(pi/15)
            1/8*sqrt(5) + 1/4*sqrt(3/2*sqrt(5) + 15/2) - 1/8
            sage: cos(pi/24)
            1/4*sqrt(2*sqrt(6) + 2*sqrt(2) + 8)
            sage: tan(pi/5)
            sqrt(-2*sqrt(5) + 5)
            sage: tan(pi/8)
            sqrt(2) - 1
            sage: tan(pi/10)
            sqrt(-2/5*sqrt(5) + 1)
            sage: tan(pi/16)
            -sqrt(2) + sqrt(2*sqrt(2) + 4) - 1
            sage: tan(pi/20)
            sqrt(5) - 1/2*sqrt(8*sqrt(5) + 20) + 1
            sage: tan(pi/24)
            sqrt(6) - sqrt(3) + sqrt(2) - 2

            sage: all(sin(rat*pi).n(200)-sin(rat*pi,hold=True).n(200) < 1e-30 for rat in [1/5,2/5,1/30,7/30,11/30,13/30,1/8,3/8,1/24,5/24,7/24,11/24])
            True
            sage: all(cos(rat*pi).n(200)-cos(rat*pi,hold=True).n(200) < 1e-30 for rat in [1/10,3/10,1/12,5/12,1/15,2/15,4/15,7/15,1/8,3/8,1/24,5/24,11/24])
            True
            sage: all(tan(rat*pi).n(200)-tan(rat*pi,hold=True).n(200) < 1e-30 for rat in [1/5,2/5,1/10,3/10,1/20,3/20,7/20,9/20,1/8,3/8,1/16,3/16,5/16,7/16,1/24,5/24,7/24,11/24])
            True
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

        We can prevent evaluation using the ``hold`` parameter::

            sage: cos(0,hold=True)
            cos(0)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = cos(0,hold=True); a.simplify()
            1

        TESTS::

            sage: conjugate(cos(x))
            cos(conjugate(x))
            sage: cos(complex(1,1))     # rel tol 1e-15
            (0.8337300251311491-0.9888977057628651j)

        """
        GinacFunction.__init__(self, "cos", latex_name=r"\cos",
                conversions=dict(maxima='cos',mathematica='Cos'))

cos = Function_cos()

class Function_tan(GinacFunction):
    def __init__(self):
        """
        The tangent function.

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

        We can prevent evaluation using the ``hold`` parameter::

            sage: tan(pi/4,hold=True)
            tan(1/4*pi)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = tan(pi/4,hold=True); a.simplify()
            1

        TESTS::

            sage: conjugate(tan(x))
            tan(conjugate(x))
            sage: tan(complex(1,1))     # rel tol 1e-15
            (0.2717525853195118+1.0839233273386946j)

        Check that :trac:`19791` is fixed::

            sage: tan(2+I).imag().n()
            1.16673625724092
        """
        GinacFunction.__init__(self, "tan", latex_name=r"\tan")

tan = Function_tan()

class Function_cot(GinacFunction):
    def __init__(self):
        r"""
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

        We can prevent evaluation using the ``hold`` parameter::

            sage: cot(pi/4,hold=True)
            cot(1/4*pi)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = cot(pi/4,hold=True); a.simplify()
            1

        EXAMPLES::

            sage: cot(pi/4)
            1
            sage: cot(x).subs(x==pi/4)
            1
            sage: cot(pi/7)
            cot(1/7*pi)
            sage: cot(x)
            cot(x)

            sage: n(cot(pi/4),100)
            1.0000000000000000000000000000
            sage: float(cot(1))
            0.64209261593433...
            sage: bool(diff(cot(x), x) == diff(1/tan(x), x))
            True
            sage: diff(cot(x), x)
            -cot(x)^2 - 1

        TESTS:

        Test complex input::

            sage: cot(complex(1,1))     # rel tol 1e-15
            (0.21762156185440273-0.8680141428959249j)
            sage: cot(1.+I)
            0.217621561854403 - 0.868014142895925*I
        """
        GinacFunction.__init__(self, "cot", latex_name=r"\cot")

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

             sage: import numpy
             sage: a = numpy.arange(2, 5)
             sage: cot(a)
             array([-0.45765755, -7.01525255,  0.86369115])
        """
        return 1 / tan(x)

cot = Function_cot()


class Function_sec(GinacFunction):
    def __init__(self):
        r"""
        The secant function.

        EXAMPLES::

            sage: sec(pi/4)
            sqrt(2)
            sage: sec(x).subs(x==pi/4)
            sqrt(2)
            sage: sec(pi/7)
            sec(1/7*pi)
            sage: sec(x)
            sec(x)
            sage: RR(sec(pi/4))
            1.41421356237310
            sage: n(sec(pi/4),100)
            1.4142135623730950488016887242
            sage: float(sec(pi/4))
            1.4142135623730951
            sage: sec(1/2)
            sec(1/2)
            sage: sec(0.5)
            1.13949392732455

            sage: bool(diff(sec(x), x) == diff(1/cos(x), x))
            True
            sage: diff(sec(x), x)
            sec(x)*tan(x)
            sage: latex(sec(x))
            \sec\left(x\right)

        We can prevent evaluation using the ``hold`` parameter::

            sage: sec(pi/4,hold=True)
            sec(1/4*pi)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = sec(pi/4,hold=True); a.simplify()
            sqrt(2)

        TESTS:

        Test complex input::

            sage: sec(complex(1,1))     # rel tol 1e-15
            (0.49833703055518686+0.5910838417210451j)
        """
        GinacFunction.__init__(self, "sec", latex_name=r"\sec")

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy
            sage: a = numpy.arange(2, 5)
            sage: sec(a)
            array([-2.40299796, -1.01010867, -1.52988566])
        """
        return 1 / cos(x)

sec = Function_sec()

class Function_csc(GinacFunction):
    def __init__(self):
        r"""
        The cosecant function.

        EXAMPLES::

            sage: csc(pi/4)
            sqrt(2)
            sage: csc(x).subs(x==pi/4)
            sqrt(2)
            sage: csc(pi/7)
            csc(1/7*pi)
            sage: csc(x)
            csc(x)
            sage: RR(csc(pi/4))
            1.41421356237310
            sage: n(csc(pi/4),100)
            1.4142135623730950488016887242
            sage: float(csc(pi/4))
            1.4142135623730951
            sage: csc(1/2)
            csc(1/2)
            sage: csc(0.5)
            2.08582964293349

            sage: bool(diff(csc(x), x) == diff(1/sin(x), x))
            True
            sage: diff(csc(x), x)
            -cot(x)*csc(x)
            sage: latex(csc(x))
            \csc\left(x\right)

        We can prevent evaluation using the ``hold`` parameter::

            sage: csc(pi/4,hold=True)
            csc(1/4*pi)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = csc(pi/4,hold=True); a.simplify()
            sqrt(2)

        TESTS:

        Test complex input::

            sage: csc(complex(1,1))     # rel tol 1e-15
            (0.6215180171704284-0.30393100162842646j)
        """
        GinacFunction.__init__(self, "csc", latex_name=r"\csc")

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy
            sage: a = numpy.arange(2, 5)
            sage: csc(a)
            array([ 1.09975017,  7.0861674 , -1.32134871])
        """
        return 1 / sin(x)

csc = Function_csc()

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

        We can delay evaluation using the ``hold`` parameter::

            sage: arcsin(0,hold=True)
            arcsin(0)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = arcsin(0,hold=True); a.simplify()
            0

        ``conjugate(arcsin(x))==arcsin(conjugate(x))``, unless on the branch
        cuts which run along the real axis outside the interval [-1, +1].::

            sage: conjugate(arcsin(x))
            conjugate(arcsin(x))
            sage: var('y', domain='positive')
            y
            sage: conjugate(arcsin(y))
            conjugate(arcsin(y))
            sage: conjugate(arcsin(y+I))
            conjugate(arcsin(y + I))
            sage: conjugate(arcsin(1/16))
            arcsin(1/16)
            sage: conjugate(arcsin(2))
            conjugate(arcsin(2))
            sage: conjugate(arcsin(-2))
            -conjugate(arcsin(2))

        TESTS::

            sage: arcsin(x).operator()
            arcsin
            sage: asin(complex(1,1))
            (0.6662394324925152+1.0612750619050357j)
        """
        GinacFunction.__init__(self, 'arcsin', latex_name=r"\arcsin",
                conversions=dict(maxima='asin', sympy='asin'))

arcsin = asin = Function_arcsin()

class Function_arccos(GinacFunction):
    def __init__(self):
        """
        The arccosine function.

        EXAMPLES::

            sage: arccos(0.5)
            1.04719755119660
            sage: arccos(1/2)
            1/3*pi
            sage: arccos(1 + 1.0*I)
            0.904556894302381 - 1.06127506190504*I
            sage: arccos(3/4).n(100)
            0.72273424781341561117837735264

        We can delay evaluation using the ``hold`` parameter::

            sage: arccos(0,hold=True)
            arccos(0)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = arccos(0,hold=True); a.simplify()
            1/2*pi

        ``conjugate(arccos(x))==arccos(conjugate(x))``, unless on the branch
        cuts, which run along the real axis outside the interval [-1, +1].::

            sage: conjugate(arccos(x))
            conjugate(arccos(x))
            sage: var('y', domain='positive')
            y
            sage: conjugate(arccos(y))
            conjugate(arccos(y))
            sage: conjugate(arccos(y+I))
            conjugate(arccos(y + I))
            sage: conjugate(arccos(1/16))
            arccos(1/16)
            sage: conjugate(arccos(2))
            conjugate(arccos(2))
            sage: conjugate(arccos(-2))
            pi - conjugate(arccos(2))

        TESTS::

            sage: arccos(x).operator()
            arccos
            sage: acos(complex(1,1))
            (0.9045568943023814-1.0612750619050357j)
        """
        GinacFunction.__init__(self, 'arccos', latex_name=r"\arccos",
                conversions=dict(maxima='acos', sympy='acos'))

arccos = acos = Function_arccos()

class Function_arctan(GinacFunction):
    def __init__(self):
        """
        The arctangent function.

        EXAMPLES::

            sage: arctan(1/2)
            arctan(1/2)
            sage: RDF(arctan(1/2))  # rel tol 1e-15
            0.46364760900080615
            sage: arctan(1 + I)
            arctan(I + 1)
            sage: arctan(1/2).n(100)
            0.46364760900080611621425623146

        We can delay evaluation using the ``hold`` parameter::

            sage: arctan(0,hold=True)
            arctan(0)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = arctan(0,hold=True); a.simplify()
            0

        ``conjugate(arctan(x))==arctan(conjugate(x))``, unless on the branch
        cuts which run along the imaginary axis outside the interval [-I, +I].::

            sage: conjugate(arctan(x))
            conjugate(arctan(x))
            sage: var('y', domain='positive')
            y
            sage: conjugate(arctan(y))
            arctan(y)
            sage: conjugate(arctan(y+I))
            conjugate(arctan(y + I))
            sage: conjugate(arctan(1/16))
            arctan(1/16)
            sage: conjugate(arctan(-2*I))
            conjugate(arctan(-2*I))
            sage: conjugate(arctan(2*I))
            conjugate(arctan(2*I))
            sage: conjugate(arctan(I/2))
            arctan(-1/2*I)

        TESTS::

            sage: arctan(x).operator()
            arctan
            sage: atan(complex(1,1))
            (1.0172219678978514+0.4023594781085251j)

        Check that :trac:`19918` is fixed::

            sage: arctan(-x).subs(x=oo)
            -1/2*pi
            sage: arctan(-x).subs(x=-oo)
            1/2*pi
        """
        GinacFunction.__init__(self, "arctan", latex_name=r'\arctan',
                conversions=dict(maxima='atan', sympy='atan'))

arctan = atan = Function_arctan()

class Function_arccot(GinacFunction):
    def __init__(self):
        """
        The arccotangent function.

        EXAMPLES::

            sage: arccot(1/2)
            arccot(1/2)
            sage: RDF(arccot(1/2))  # abs tol 2e-16
            1.1071487177940906
            sage: arccot(1 + I)
            arccot(I + 1)
            sage: arccot(1/2).n(100)
            1.1071487177940905030170654602
            sage: float(arccot(1/2))  # abs tol 2e-16
            1.1071487177940906
            sage: bool(diff(acot(x), x) == -diff(atan(x), x))
            True
            sage: diff(acot(x), x)
            -1/(x^2 + 1)

        We can delay evaluation using the ``hold`` parameter::

            sage: arccot(1,hold=True)
            arccot(1)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = arccot(1,hold=True); a.simplify()
            1/4*pi

        TESTS:

        Test complex input::

            sage: arccot(complex(1,1))  # rel tol 1e-15
            (0.5535743588970452-0.4023594781085251j)
            sage: arccot(1.+I)
            0.553574358897045 - 0.402359478108525*I

        """
        GinacFunction.__init__(self, "arccot", latex_name=r'{\rm arccot}',
                conversions=dict(maxima='acot', sympy='acot'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy
            sage: a = numpy.arange(2, 5)
            sage: arccot(a)
            array([ 0.46364761,  0.32175055,  0.24497866])
        """
        return math.pi/2 - arctan(x)

arccot = acot = Function_arccot()

class Function_arccsc(GinacFunction):
    def __init__(self):
        """
        The arccosecant function.

        EXAMPLES::

            sage: arccsc(2)
            arccsc(2)
            sage: RDF(arccsc(2))  # rel tol 1e-15
            0.5235987755982988
            sage: arccsc(2).n(100)
            0.52359877559829887307710723055
            sage: float(arccsc(2))
            0.52359877559829...
            sage: arccsc(1 + I)
            arccsc(I + 1)
            sage: diff(acsc(x), x)
            -1/(sqrt(x^2 - 1)*x)

        We can delay evaluation using the ``hold`` parameter::

            sage: arccsc(1,hold=True)
            arccsc(1)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = arccsc(1,hold=True); a.simplify()
            1/2*pi

        TESTS:

        Test complex input::

            sage: arccsc(complex(1,1))  # rel tol 1e-15
            (0.45227844715119064-0.5306375309525178j)
        """
        GinacFunction.__init__(self, "arccsc", latex_name=r'{\rm arccsc}',
                                   conversions=dict(maxima='acsc'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy
            sage: a = numpy.arange(2, 5)
            sage: arccsc(a)
            array([ 0.52359878,  0.33983691,  0.25268026])
        """
        return arcsin(1.0/x)

arccsc = acsc = Function_arccsc()

class Function_arcsec(GinacFunction):
    def __init__(self):
        """
        The arcsecant function.

        EXAMPLES::

            sage: arcsec(2)
            arcsec(2)
            sage: arcsec(2.0)
            1.04719755119660
            sage: arcsec(2).n(100)
            1.0471975511965977461542144611
            sage: arcsec(1/2).n(100)
            NaN
            sage: RDF(arcsec(2))  # abs tol 1e-15
            1.0471975511965976
            sage: arcsec(1 + I)
            arcsec(I + 1)
            sage: diff(asec(x), x)
            1/(sqrt(x^2 - 1)*x)

        We can delay evaluation using the ``hold`` parameter::

            sage: arcsec(1,hold=True)
            arcsec(1)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = arcsec(1,hold=True); a.simplify()
            0

        TESTS:

        Test complex input::

            sage: arcsec(complex(1,1))  # rel tol 1e-15
            (1.118517879643706+0.5306375309525178j)
        """
        GinacFunction.__init__(self, "arcsec", latex_name=r'{\rm arcsec}',
                                   conversions=dict(maxima='asec'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy
            sage: a = numpy.arange(2, 5)
            sage: arcsec(a)
            array([ 1.04719755,  1.23095942,  1.31811607])
        """
        return arccos(1.0/x)

arcsec = asec = Function_arcsec()

class Function_arctan2(GinacFunction):
    def __init__(self):
        r"""
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
            2.356194490192345

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

        We can delay evaluation using the ``hold`` parameter::

            sage: arctan2(-1/2,1,hold=True)
            arctan2(-1/2, 1)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: arctan2(-1/2,1,hold=True).simplify()
            -arctan(1/2)

        The function also works with numpy arrays as input::

            sage: import numpy
            sage: a = numpy.linspace(1, 3, 3)
            sage: b = numpy.linspace(3, 6, 3)
            sage: atan2(a, b)
            array([ 0.32175055,  0.41822433,  0.46364761])

            sage: atan2(1,a)
            array([ 0.78539816,  0.46364761,  0.32175055])

            sage: atan2(a, 1)
            array([ 0.78539816,  1.10714872,  1.24904577])

        TESTS::

            sage: x,y = var('x,y')
            sage: arctan2(y,x).operator()
            arctan2

        Check if :trac:`8565` is fixed::

            sage: atan2(-pi,0)
            -1/2*pi

        Check if :trac:`8564` is fixed::

            sage: arctan2(x,x)._sympy_()
            atan2(x, x)

        Check if numerical evaluation works :trac:`9913`::

            sage: arctan2(0, -log(2)).n()
            3.14159265358979

        Check if atan2(0,0) throws error of :trac:`11423`::

            sage: atan2(0,0)
            Traceback (most recent call last):
            ...
            RuntimeError: arctan2_eval(): arctan2(0,0) encountered

            sage: atan2(0,0,hold=True)
            arctan2(0, 0)

            sage: atan2(0,0,hold=True).n()
            Traceback (most recent call last):
            ...
            ValueError: arctan2(0,0) undefined

        Check if :trac:`10062` is fixed, this was caused by
        ``(I*I).is_positive()`` returning ``True``::

            sage: arctan2(0, I*I)
            pi
        """
        GinacFunction.__init__(self, "arctan2", nargs=2, latex_name=r'\arctan',
                conversions=dict(maxima='atan2', sympy='atan2'))

arctan2 = atan2 = Function_arctan2()
