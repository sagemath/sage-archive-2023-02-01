"""
Logarithmic functions
"""
from sage.symbolic.function import GinacFunction, BuiltinFunction
from sage.libs.pari.gen import pari, PariError
from sage.rings.all import CDF, Integer, Rational, RealField, ComplexField
from sage.symbolic.all import SR
import math

class Function_exp(GinacFunction):
    def __init__(self):
        r"""
        The exponential function, `\exp(x) = e^x`.

        EXAMPLES::

            sage: exp(-1)
            e^(-1)
            sage: exp(2)
            e^2
            sage: exp(2).n(100)
            7.3890560989306502272304274606
            sage: exp(x^2 + log(x))
            e^(x^2 + log(x))
            sage: exp(x^2 + log(x)).simplify()
            x*e^(x^2)
            sage: exp(2.5)
            12.1824939607035
            sage: exp(float(2.5))
            12.182493960703473
            sage: exp(RDF('2.5'))
            12.1824939607

        ::

            sage: exp(pi*I/2)
            I
            sage: exp(pi*I)
            -1
            sage: exp(8*pi*I)
            1
            sage: exp(7*pi*I/2)
            -I

        TEST::

            sage: latex(exp(x))
            e^{x}
            sage: latex(exp(sqrt(x)))
            e^{\sqrt{x}}
            sage: latex(exp)
            \exp
            sage: latex(exp(sqrt(x))^x)
            \left(e^{\sqrt{x}}\right)^{x}
            sage: latex(exp(sqrt(x)^x))
            e^{\left(\sqrt{x}^{x}\right)}

        Test simplifications when taking powers of exp, #7264::

            sage: var('a,b,c,I')
            (a, b, c, I)
            sage: model_exp = exp(I)**a*(b)
            sage: sol1_l={b: 5.0, a: 1.1}
            sage: model_exp.subs(sol1_l)
            5.00000000000000*(e^I)^1.10000000000000

        ::

            sage: exp(3)^I*exp(x)
            (e^3)^I*e^x
            sage: exp(x)*exp(x)
            e^(2*x)
            sage: exp(x)*exp(a)
            e^(a + x)
            sage: exp(x)*exp(a)^2
            e^(2*a + x)

        Another instance of the same problem, #7394::

            sage: 2*sqrt(e)
            2*sqrt(e)
        """
        GinacFunction.__init__(self, "exp", latex_name=r"\exp",
                                   conversions=dict(maxima='exp'))

    def __call__(self, x, coerce=True, hold=False, prec=None,
            dont_call_method_on_arg=False):
        """
        Note that the ``prec`` argument is deprecated. The precision for
        the result is deduced from the precision of the input. Convert
        the input to a higher precision explicitly if a result with higher
        precision is desired.::

            sage: t = exp(RealField(100)(2)); t
            7.3890560989306502272304274606
            sage: t.prec()
            100

        TESTS::

            sage: exp(2,prec=100)
            doctest:...: DeprecationWarning: The prec keyword argument is deprecated. Explicitly set the precision of the input, for example exp(RealField(300)(1)), or use the prec argument to .n() for exact inputs, e.g., exp(1).n(300), instead.
            7.3890560989306502272304274606
        """
        if prec is not None:
            from sage.misc.misc import deprecation
            deprecation("The prec keyword argument is deprecated. Explicitly set the precision of the input, for example exp(RealField(300)(1)), or use the prec argument to .n() for exact inputs, e.g., exp(1).n(300), instead.")
            x = GinacFunction.__call__(self, x, coerce=coerce, hold=hold,
                    dont_call_method_on_arg=dont_call_method_on_arg)
            return x.n(prec)
        return GinacFunction.__call__(self, x, coerce=coerce, hold=hold,
                dont_call_method_on_arg=dont_call_method_on_arg)

exp = Function_exp()

class Function_log(GinacFunction):
    def __init__(self):
        r"""
        The natural logarithm of x.  See `log?` for
        more information about its behavior.

        EXAMPLES::

            sage: ln(e^2)
            2
            sage: ln(2)
            log(2)
            sage: ln(10)
            log(10)

        ::

            sage: ln(RDF(10))
            2.30258509299
            sage: ln(2.718)
            0.999896315728952
            sage: ln(2.0)
            0.693147180559945
            sage: ln(float(-1))
            3.1415926535897931j
            sage: ln(complex(-1))
            3.1415926535897931j

        TESTS::

            sage: latex(x.log())
            \log\left(x\right)
            sage: latex(log(1/4))
            \log\left(\frac{1}{4}\right)
            sage: loads(dumps(ln(x)+1))
            log(x) + 1

        Check if float arguments are handled properly.::

            sage: from sage.functions.log import function_log as log
            sage: log(float(5))
            1.6094379124341003
            sage: log(float(0))
            -inf
            sage: log(float(-1))
            3.1415926535897931j
            sage: log(x).subs(x=float(-1))
            3.1415926535897931j
        """
        GinacFunction.__init__(self, 'log', latex_name=r'\log',
                                   conversions=dict(maxima='log'))

ln = function_log = Function_log()

def log(x, base=None):
    """
    Return the logarithm of x to the given base.

    Calls the ``log`` method of the object x when computing
    the logarithm, thus allowing use of logarithm on any object
    containing a ``log`` method. In other words, log works
    on more than just real numbers.

    EXAMPLES::

        sage: log(e^2)
        2
        sage: log(1024, 2); RDF(log(1024, 2))
        10
        10.0
        sage: log(10, 4); RDF(log(10, 4))
        log(10)/log(4)
        1.66096404744

    ::

        sage: log(10, 2)
        log(10)/log(2)
        sage: n(log(10, 2))
        3.32192809488736
        sage: log(10, e)
        log(10)
        sage: n(log(10, e))
        2.30258509299405

    The log function works for negative numbers, complex
    numbers, and symbolic numbers too, picking the branch
    with angle between `-pi` and `pi`::

        sage: log(-1+0*I)
        I*pi
        sage: log(CC(-1))
        3.14159265358979*I
        sage: log(-1.0)
        3.14159265358979*I

    For input zero, the following behavior occurs::

        sage: log(0)
        -Infinity
        sage: log(CC(0))
        -infinity
        sage: log(0.0)
        -infinity

    The log function also works in finite fields as long as the base is
    generator of the multiplicative group::

        sage: F = GF(13); g = F.multiplicative_generator(); g
        2
        sage: a = F(8)
        sage: log(a,g); g^log(a,g)
        3
        8
        sage: log(a,3)
        Traceback (most recent call last):
        ...
        ValueError: base (=3) for discrete log must generate multiplicative group

    The log function also works for p-adics (see documentation for
    p-adics for more information)::

        sage: R = Zp(5); R
        5-adic Ring with capped relative precision 20
        sage: a = R(16); a
        1 + 3*5 + O(5^20)
        sage: log(a)
        3*5 + 3*5^2 + 3*5^4 + 3*5^5 + 3*5^6 + 4*5^7 + 2*5^8 + 5^9 + 5^11 + 2*5^12 + 5^13 + 3*5^15 + 2*5^16 + 4*5^17 + 3*5^18 + 3*5^19 + O(5^20)
    """
    if base is None:
        try:
            return x.log()
        except AttributeError:
            return ln(x)
    else:
        try:
            return x.log(base)
        except (AttributeError, TypeError):
            return log(x) / log(base)


class Function_polylog(GinacFunction):
    def __init__(self):
        r"""
        The polylog function
        `\text{Li}_n(z) = \sum_{k=1}^{\infty} z^k / k^n`.

        INPUT:

        -  ``n`` - object
        -  ``z`` - object

        EXAMPLES::

            sage: polylog(1, x)
            -log(-x + 1)
            sage: polylog(2,1)
            1/6*pi^2
            sage: polylog(2,x^2+1)
            polylog(2, x^2 + 1)
            sage: polylog(4,0.5)
            polylog(4, 0.500000000000000)

            sage: f = polylog(4, 1); f
            1/90*pi^4
            sage: f.n()
            1.08232323371114

            sage: polylog(4, 2).n()
            2.42786280675470 - 0.174371300025453*I
            sage: complex(polylog(4,2))
            (2.4278628067547032-0.17437130002545306j)
            sage: float(polylog(4,0.5))
            0.51747906167389934

            sage: z = var('z')
            sage: polylog(2,z).series(z==0, 5)
            1*z + 1/4*z^2 + 1/9*z^3 + 1/16*z^4 + Order(z^5)

            sage: loads(dumps(polylog))
            polylog

            sage: latex(polylog(5, x))
            {\rm Li}_{5}(x)

        TESTS:

        Check if #8459 is fixed::

            sage: t = maxima(polylog(5,x)).sage(); t
            polylog(5, x)
            sage: t.operator() == polylog
            True
            sage: t.subs(x=.5).n()
            0.508400579242269
        """
        GinacFunction.__init__(self, "polylog", nargs=2)

    def _maxima_init_evaled_(self, *args):
        """
        EXAMPLES:

        These are indirect doctests for this function.::

            sage: polylog(2, x)._maxima_()
            li[2](x)
            sage: polylog(4, x)._maxima_()
            polylog(4,x)
        """
        n, x = args
        if int(n) in [1,2,3]:
            return 'li[%s](%s)'%(n, x)
        else:
            return 'polylog(%s, %s)'%(n, x)


polylog = Function_polylog()

class Function_dilog(GinacFunction):
    def __init__(self):
        r"""
        The dilogarithm function
        `\text{Li}_2(z) = \sum_{k=1}^{\infty} z^k / k^2`.

        This is simply an alias for polylog(2, z).

        EXAMPLES::

            sage: dilog(1)
            1/6*pi^2
            sage: dilog(1/2)
            1/12*pi^2 - 1/2*log(2)^2
            sage: dilog(x^2+1)
            dilog(x^2 + 1)
            sage: float(dilog(1))
            1.6449340668482262
            sage: var('z')
            z
            sage: dilog(z).diff(z, 2)
            log(-z + 1)/z^2 - 1/((z - 1)*z)
            sage: dilog(z).series(z==1/2, 3)
            (1/12*pi^2 - 1/2*log(2)^2) + (-2*log(1/2))*(z - 1/2) + (2*log(1/2) + 2)*(z - 1/2)^2 + Order(1/8*(2*z - 1)^3)

            sage: latex(dilog(z))
            {\rm Li}_2\left(z\right)
        """
        GinacFunction.__init__(self, 'dilog',
                conversions=dict(maxima='li[2]'))

dilog = Function_dilog()

