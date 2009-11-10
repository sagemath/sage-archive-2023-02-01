r"""
Generalized Functions

Sage implements several generalized functions (also known as
distributions) such as Dirac delta, Heaviside step functions. These
generalized functions can be manipulated within Sage like any other
symbolic functions.


AUTHORS:

- Golam Mortuza Hossain (2009-06-26): initial version


EXAMPLES:

Dirac delta function::

    sage: dirac_delta(x)
    dirac_delta(x)

Heaviside step function::

    sage: heaviside(x)
    heaviside(x)

Unit step function::

    sage: unit_step(x)
    unit_step(x)

Signum (sgn) function::

    sage: sgn(x)
    sgn(x)

Kronecker delta function::

    sage: m,n=var('m,n')
    sage: kronecker_delta(m,n)
    kronecker_delta(m, n)

"""

##############################################################################
#
#       Copyright (C) 2009 Golam Mortuza Hossain <gmhossain@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL v2+)
#                  http://www.gnu.org/licenses/
#
##############################################################################

from sage.symbolic.function import SFunction, PrimitiveFunction
from sage.rings.complex_interval_field import ComplexIntervalField

class FunctionDiracDelta(PrimitiveFunction):
    r"""
    The Dirac delta (generalized) function, `\delta(x)` (``dirac_delta(x)``).

    INPUT:

    -  ``x`` - a real number or a symbolic expression

    DEFINITION:

    Dirac delta function `\delta(x)`, is defined in Sage as:

        `\delta(x) = 0` for real `x \ne 0` and
        `\int_{-\infty}^{\infty} \delta(x) dx = 1`

    Its alternate definition with respect to an arbitrary test
    function `f(x)` is

        `\int_{-\infty}^{\infty} f(x) \delta(x-a) dx = f(a)`

    EXAMPLES::

        sage: dirac_delta(1)
        0
        sage: dirac_delta(0)
        dirac_delta(0)
        sage: dirac_delta(x)
        dirac_delta(x)

    REFERENCES:

    -  http://en.wikipedia.org/wiki/Dirac_delta_function

    """
    def __init__(self):
        r"""
        The Dirac delta (generalized) function, ``dirac_delta(x)``.

        INPUT:

        -  ``x`` - a real number or a symbolic expression

        EXAMPLES::

            sage: dirac_delta(1)
            0
            sage: dirac_delta(0)
            dirac_delta(0)
            sage: dirac_delta(x)
            dirac_delta(x)
        """
        PrimitiveFunction.__init__(self, "dirac_delta", latex=r"\delta",
                                   conversions=dict(maxima='delta',
                                    mathematica='DiracDelta'))

    __call__ = SFunction.__call__

    def _eval_(self, x):
        """
        INPUT:

        -  ``x`` - a real number or a symbolic expression

        EXAMPLES::

            sage: dirac_delta(1)
            0
            sage: dirac_delta(0)
            dirac_delta(0)
            sage: dirac_delta(x)
            dirac_delta(x)
            sage: dirac_delta(exp(-10000000000000000000))
            0

        Evaluation test::

            sage: dirac_delta(x).subs(x=1)
            0
        """
        try:
            approx_x = ComplexIntervalField()(x)
            if bool(approx_x.imag() == 0):      # x is real
                if bool(approx_x.real() == 0):  # x is zero
                    return None
                else:
                    return 0
        except:                     # x is symbolic
            pass
        return None

dirac_delta = FunctionDiracDelta()

class FunctionHeaviside(PrimitiveFunction):
    r"""
    The Heaviside step function, `H(x)` (``heaviside(x)``).

    INPUT:

    -  ``x`` - a real number or a symbolic expression

    DEFINITION:

    The Heaviside step function, `H(x)` is defined in Sage as:

        `H(x) = 0` for `x < 0` and `H(x) = 1` for `x > 0`

    EXAMPLES::

        sage: heaviside(-1)
        0
        sage: heaviside(1)
        1
        sage: heaviside(0)
        heaviside(0)
        sage: heaviside(x)
        heaviside(x)

    REFERENCES:

    -  http://en.wikipedia.org/wiki/Heaviside_function

    """
    def __init__(self):
        r"""
        The Heaviside step function, ``heaviside(x)``.

        INPUT:

        -  ``x`` - a real number or a symbolic expression

        EXAMPLES::

            sage: heaviside(-1)
            0
            sage: heaviside(1)
            1
            sage: heaviside(0)
            heaviside(0)
            sage: heaviside(x)
            heaviside(x)
        """
        PrimitiveFunction.__init__(self, "heaviside", latex="H",
                                   conversions=dict(mathematica='HeavisideTheta'))

    __call__ = SFunction.__call__

    def _eval_(self, x):
        """
        INPUT:

        -  ``x`` - a real number or a symbolic expression

        EXAMPLES::

            sage: heaviside(-1/2)
            0
            sage: heaviside(1)
            1
            sage: heaviside(0)
            heaviside(0)
            sage: heaviside(x)
            heaviside(x)
            sage: heaviside(exp(-1000000000000000000000))
            1

        Evaluation test::

            sage: heaviside(x).subs(x=1)
            1
            sage: heaviside(x).subs(x=-1)
            0
        """
        try:
            approx_x = ComplexIntervalField()(x)
            if bool(approx_x.imag() == 0):      # x is real
                if bool(approx_x.real() == 0):  # x is zero
                    return None
                # Now we have a non-zero real
                if bool((approx_x**(0.5)).imag() == 0): # Check: x > 0
                    return 1
                else:
                    return 0
        except:                     # x is symbolic
            pass
        return None

    def _derivative_(self, *args, **kwds):
        """
        Derivative of Heaviside step function

        EXAMPLES::

            sage: heaviside(x).diff(x)
            dirac_delta(x)
        """
        diff_param = kwds['diff_param']
        assert diff_param == 0
        x = args[diff_param]
        return dirac_delta(x)

heaviside = FunctionHeaviside()

class FunctionUnitStep(PrimitiveFunction):
    r"""
    The unit step function, `\mathrm{u}(x)` (``unit_step(x)``).

    INPUT:

    -  ``x`` - a real number or a symbolic expression

    DEFINITION:

    The unit step function, `\mathrm{u}(x)` is defined in Sage as:

        `\mathrm{u}(x) = 0` for `x < 0` and `\mathrm{u}(x) = 1` for `x \geq 0`

    EXAMPLES::

        sage: unit_step(-1)
        0
        sage: unit_step(1)
        1
        sage: unit_step(0)
        1
        sage: unit_step(x)
        unit_step(x)
    """
    def __init__(self):
        r"""
        The unit step function, ``unit_step(x)``.

        INPUT:

        -  ``x`` - a real number or a symbolic expression

        EXAMPLES:

            sage: unit_step(-1)
            0
            sage: unit_step(1)
            1
            sage: unit_step(0)
            1
            sage: unit_step(x)
            unit_step(x)
        """
        PrimitiveFunction.__init__(self, "unit_step", latex=r"\mathrm{u}",
                                   conversions=dict(mathematica='UnitStep'))

    __call__ = SFunction.__call__

    def _eval_(self, x):
        """
        INPUT:

        -  ``x`` - a real number or a symbolic expression

        EXAMPLES::

            sage: unit_step(-1)
            0
            sage: unit_step(1)
            1
            sage: unit_step(0)
            1
            sage: unit_step(x)
            unit_step(x)
            sage: unit_step(-exp(-10000000000000000000))
            0

        Evaluation test::

            sage: unit_step(x).subs(x=1)
            1
            sage: unit_step(x).subs(x=0)
            1
        """
        try:
            approx_x = ComplexIntervalField()(x)
            if bool(approx_x.imag() == 0):      # x is real
                if bool(approx_x.real() == 0):  # x is zero
                    return 1
                # Now we have a non-zero real
                if bool((approx_x**(0.5)).imag() == 0): # Check: x > 0
                    return 1
                else:
                    return 0
        except:                     # x is symbolic
            pass
        return None

    def _derivative_(self, *args, **kwds):
        """
        Derivative of unit step function

        EXAMPLES::

            sage: unit_step(x).diff(x)
            dirac_delta(x)
        """
        diff_param = kwds['diff_param']
        assert diff_param == 0
        x = args[diff_param]
        return dirac_delta(x)

unit_step = FunctionUnitStep()

class FunctionSignum(PrimitiveFunction):
    r"""
    The signum or sgn function `\mathrm{sgn}(x)` (``sgn(x)``).

    INPUT:

    -  ``x`` - a real number or a symbolic expression

    DEFINITION:

    The sgn function, `\mathrm{sgn}(x)` is defined as:

        `\mathrm{sgn}(x) =  1` for `x > 0`,
        `\mathrm{sgn}(x) =  0` for `x = 0` and
        `\mathrm{sgn}(x) = -1` for `x < 0`

    EXAMPLES::

        sage: sgn(-1)
        -1
        sage: sgn(1)
        1
        sage: sgn(0)
        0
        sage: sgn(x)
        sgn(x)

    REFERENCES:

    -  http://en.wikipedia.org/wiki/Sign_function

    """
    def __init__(self):
        r"""
        The sgn function, ``sgn(x)``.

        EXAMPLES:

            sage: sgn(-1)
            -1
            sage: sgn(1)
            1
            sage: sgn(0)
            0
            sage: sgn(x)
            sgn(x)
        """
        PrimitiveFunction.__init__(self, "sgn", latex=r"\mathrm{sgn}",
                                   conversions=dict(maxima='signum',mathematica='Sign'))

    __call__  = SFunction.__call__

    def _eval_(self, x, **kwds):
        """

        EXAMPLES::

            sage: sgn(-1)
            -1
            sage: sgn(1)
            1
            sage: sgn(0)
            0
            sage: sgn(x)
            sgn(x)
            sage: sgn(-exp(-10000000000000000000))
            -1

        Evaluation test::

            sage: sgn(x).subs(x=1)
            1
            sage: sgn(x).subs(x=0)
            0
            sage: sgn(x).subs(x=-1)
            -1
        """
        try:
            approx_x = ComplexIntervalField()(x)
            if bool(approx_x.imag() == 0):      # x is real
                if bool(approx_x.real() == 0):  # x is zero
                    return 0
                # Now we have a non-zero real
                if bool((approx_x**(0.5)).imag() == 0): # Check: x > 0
                    return 1
                else:
                    return -1
        except:                     # x is symbolic
            pass
        return None

    def _derivative_(self, *args, **kwds):
        """
        Derivative of sgn function

        EXAMPLES::

            sage: sgn(x).diff(x)
            2*dirac_delta(x)
        """
        diff_param = kwds['diff_param']
        assert diff_param == 0
        x = args[diff_param]
        return 2*dirac_delta(x)

sgn = FunctionSignum()

class FunctionKroneckerDelta(PrimitiveFunction):
    r"""
    The Kronecker delta function `\delta_{m,n}` (``kronecker_delta(m, n)``).

    INPUT:

    -  ``m`` - a number or a symbolic expression
    -  ``n`` - a number or a symbolic expression

    DEFINITION:

    Kronecker delta function `\delta_{m,n}` is defined as:

        `\delta_{m,n} = 0` for `m \ne n` and
        `\delta_{m,n} = 1` for `m = n`

    EXAMPLES::

        sage: kronecker_delta(1,2)
        0
        sage: kronecker_delta(1,1)
        1
        sage: m,n=var('m,n')
        sage: kronecker_delta(m,n)
        kronecker_delta(m, n)

    REFERENCES:

    - http://en.wikipedia.org/wiki/Kronecker_delta

    """
    def __init__(self):
        r"""
        The Kronecker delta function.

        EXAMPLES::

            sage: kronecker_delta(1,2)
            0
            sage: kronecker_delta(1,1)
            1
        """
        PrimitiveFunction.__init__(self, "kronecker_delta", nargs=2,
                                        conversions=dict(maxima='kron_delta',
                                        mathematica='KroneckerDelta'))

    def __call__(self, m, n, **kwds):
        """
        The Kronecker delta function.

        EXAMPLES::

            sage: x,y=var('x,y')
            sage: kronecker_delta(x, y)
            kronecker_delta(x, y)
            sage: kronecker_delta(y, x)
            kronecker_delta(x, y)
        """
        # It is a symmetric function. So lets keep arguments sorted. This
        # ensures symbolic (k_d(m, n) - k_d(n, m) cancels automatically.
        if bool(repr(m) > repr(n)):
            return SFunction.__call__(self, n, m)
        return SFunction.__call__(self, m, n)

    def _eval_(self, m, n, **kwds):
        """
        The Kronecker delta function.

        EXAMPLES::

            sage: kronecker_delta(1,2)
            0
            sage: kronecker_delta(1,1)
            1
            sage: kronecker_delta(x,2*x)
            kronecker_delta(2*x, x)

        Evaluation test::

            sage: kronecker_delta(1,x).subs(x=1)
            1
        """
        x = m - n
        try:
            approx_x = ComplexIntervalField()(x)
            if bool(approx_x.imag() == 0):      # x is real
                if bool(approx_x.real() == 0):  # x is zero
                    return 1
                else:
                    return 0
            else:
                return 0            # x is complex
        except:                     # x is symbolic
            pass
        return None

    def _derivative_(self, *args, **kwds):
        """
        Derivative of Kronecker delta

        EXAMPLES::

            sage: kronecker_delta(x,1).diff(x)
            0
        """
        # Kronecker delta is non-zero (but finite) only in the set of
        # zero-measure unlike Dirac delta. Consequently, it is null
        # for the purpose of integration/differentiation. For *discrete sum*
        # Kronecker delta is however non-trivial.
        return 0

    def _print_latex_(self, m, n, **kwds):
        """
        Return latex expression

        EXAMPLES::

            sage: from sage.misc.latex import latex
            sage: m,n=var('m,n')
            sage: latex(kronecker_delta(m,n))
            \delta_{m,n}

        """
        from sage.misc.latex import latex
        return "\\delta_{%s,%s}"%(latex(m), latex(n))

kronecker_delta = FunctionKroneckerDelta()

