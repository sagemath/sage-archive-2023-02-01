"""
Logarithmic functions
"""
from sage.symbolic.function import SFunction, PrimitiveFunction
from sage.libs.pari.gen import pari, PariError
from sage.rings.all import CDF, Integer, Rational, RealField, ComplexField
from sage.symbolic.all import SR
import math

class Function_exp(PrimitiveFunction):
    def __init__(self):
        r"""
        The exponential function, `\exp(x) = e^x`.

        EXAMPLES::

            sage: exp(-1)
            e^(-1)
            sage: exp(2)
            e^2
            sage: exp(2,prec=100)
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
        """
        PrimitiveFunction.__init__(self, "exp", latex=r"\exp",
                                   conversions=dict(maxima='exp'),
                                   approx=math.exp)

    def __call__(self, x, prec=None):
        """
        INPUT:

        -  ``x`` - a number

        -  ``prec`` - integer (default: None): if None, returns
           an exact exponential; otherwise returns a numerical exponential if
           necessary, to the given bits of precision.

        EXAMPLES::

            sage: exp(pi*I/2)
            I
            sage: exp(pi*I)
            -1
            sage: exp(8*pi*I)
            1
            sage: exp(7*pi*I/2)
            -I
        """
        if isinstance(x, float):
            return self._approx_(x)
        if not isinstance(x, (Integer, Rational)):
            try:
                return x.exp()
            except AttributeError:
                pass
        if prec:
            return RealField(prec)(x).exp()
        return SFunction.__call__(self, SR(x))

exp = Function_exp()

class Function_log(PrimitiveFunction):
    def __init__(self):
        """
        The log function.

        EXAMPLES::

            sage: log(e^2)
            2
            sage: log(1024, 2)
            10
            sage: log(10, 4)
            log(10)/log(4)

        ::

            sage: RDF(log(10,2))
            3.32192809489
            sage: RDF(log(8, 2))
            3.0
            sage: log(RDF(10))
            2.30258509299
            sage: log(2.718)
            0.999896315728952
        """
        PrimitiveFunction.__init__(self, 'log', latex=r'\log',
                                   conversions=dict(maxima='log'),
                                   approx=math.log)

function_log = Function_log()

def ln(x):
    """
    The natural logarithm of x.

    INPUT:

    -  ``x`` - positive real number

    OUTPUT:

    -  ``ln(x)`` - real number

    EXAMPLES::

        sage: ln(e^2)
        2
        sage: ln(2)
        log(2)
        sage: ln(2.0)
        0.693147180559945
    """
    return function_log(x)

def log(x, base=None):
    """
    Return the logarithm of x to the given base.

    Calls the ``log`` method of the object x when computing
    the logarithm, thus allowing use of logarithm on any object
    containing a ``log`` method. In other words, log works
    on more than just real numbers.

    TODO: Add p-adic log example.

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
    """
    if base is None:
        try:
            return x.log()
        except AttributeError:
            return ln(x)
    else:
        try:
            return x.log(base)
        except AttributeError:
            return ln(x) / ln(base)


class Function_polylog(PrimitiveFunction):
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
            polylog(2,x^2 + 1)
            sage: polylog(4,0.5)
            polylog(4,0.500000000000000)

            sage: f = polylog(4, 1); f
            1/90*pi^4
            sage: f.n()
            1.08232323371114

            sage: polylog(4, 2).n()
            2.42786280675470 - 0.174371300025453*I
            sage: complex(polylog(4,2))
            (2.4278628067547001-0.17437130002545301j)
            sage: float(polylog(4,0.5))
            0.51747906167389901

            sage: z = var('z')
            sage: polylog(2,z).series(z==0, 5)
            1*z + 1/4*z^2 + 1/9*z^3 + 1/16*z^4 + Order(z^5)

            sage: loads(dumps(polylog))
            polylog
        """
        PrimitiveFunction.__init__(self, "polylog", nargs=2)

    __call__ = SFunction.__call__

    def _latex_composition(self, n, x):
        """
        Return Latex representation of this polylogarithm.

        EXAMPLES::

            sage: latex(polylog(5, x))
            \text{Li}_{5}(x)
        """
        return "\\text{Li}_{%s}(%s)"%(n, x)

    def _maxima_init_evaled_(self, *args):
        """
        EXAMPLES:

        These are indirect doctests for this function.

            sage: polylog(2, x)._maxima_init_()
            'li[2](x)'
            sage: polylog(4, x)._maxima_init_()
            'polylog(4, x)'
        """
        n, x = args
        if int(n) in [1,2,3]:
            return 'li[%s](%s)'%(n, x)
        else:
            return 'polylog(%s, %s)'%(n, x)


polylog = Function_polylog()

class Function_dilog(PrimitiveFunction):
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
        """
        PrimitiveFunction.__init__(self, 'dilog',
                                   conversions=dict(ginac='dilog', maxima='li[2]'))

dilog = Function_dilog()

