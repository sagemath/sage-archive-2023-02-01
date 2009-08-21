"""
Symbolic Integration
"""

##############################################################################
#
#       Copyright (C) 2009 Golam Mortuza Hossain <gmhossain@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL v2+)
#                  http://www.gnu.org/licenses/
#
##############################################################################

from sage.symbolic.ring import SR, is_SymbolicVariable
from sage.symbolic.function import SFunction, sfunctions_funcs

# Following thin wrappers provide uniform access to different algorithm
# from _eval_ method. Probably in future these wrapper themselves can
# be made stand-alone algorithm interface by re-factoring "integral" codes
# of calculus.calculus

def _maxima_integrator(*args, **kwds):
    """
    Thin wrapper around maxima symbolic integration. It raises
    NotImplementedError if it fails.

    EXAMPLES::

        sage: from sage.symbolic.integration.integration import _maxima_integrator
        sage: _maxima_integrator(sin(x), x)
        -cos(x)
        sage: f(x) = function('f', x)
        sage: _maxima_integrator(f(x), x)
        Traceback (most recent call last):
        ...
        NotImplementedError: Maxima failed to integrate

    """
    from sage.calculus.calculus import integral as maxima_integral
    kwds['raise_error'] = True
    kwds['algorithm'] = 'maxima'
    return maxima_integral(*args, **kwds)

def _sympy_integrator(*args, **kwds):
    """
    Thin wrapper around sympy symbolic integration. It raises
    NotImplementedError if it fails.

    EXAMPLES::

        sage: from sage.symbolic.integration.integration import _sympy_integrator
        sage: _sympy_integrator(sin(x), x)
        -cos(x)
        sage: f(x) = function('f', x)
        sage: _sympy_integrator(f(x), x)
        Traceback (most recent call last):
        ...
        NotImplementedError: Sympy failed to integrate

    """
    from sage.calculus.calculus import integral as sympy_integral
    kwds['algorithm'] = 'sympy'
    # Check whether sympy could evaluate the integral
    try:
        return sympy_integral(*args, **kwds)
    except:
        raise NotImplementedError, "Sympy failed to integrate"

def _mathematica_free_integrator(*args, **kwds):
    """
    Thin wrapper around mathematica_free symbolic integration. It raises
    NotImplementedError if it fails.

    EXAMPLES::

        sage: from sage.symbolic.integration.integration import _mathematica_free_integrator
        sage: _mathematica_free_integrator(sin(x), x)
        -cos(x)
        sage: f(x) = function('f', x)
        sage: _mathematica_free_integrator(f(x), x)
        Traceback (most recent call last):
        ...
        NotImplementedError: mathematica_free failed to integrate

    """
    from sage.calculus.calculus import integral as mathematica_free_integral
    kwds['algorithm'] = 'mathematica_free'
    # Check whether mathematica_free could evaluate the integral
    try:
        return mathematica_free_integral(*args, **kwds)
    except:
        raise NotImplementedError, "mathematica_free failed to integrate"

##################################################################
#  Table of available integration algorithm:
##################################################################

# To add new integration algorithmas in Sage, list them here.
algorithm_table = {}

algorithm_table['maxima'] = _maxima_integrator
algorithm_table['sympy'] = _sympy_integrator
algorithm_table['mathematica_free'] = _mathematica_free_integrator

# List here the integrators that should be tried in a given
# order by default
default_integrators = \
    [
    _maxima_integrator,
    ]

######################################################
#
# Class implementing symbolic integration
#
######################################################

class SymbolicIntegration(SFunction):
    r"""
    Symbolic integration class for Sage. It uses either user-specified
    algorithm or default set of integrators for integration.

    EXAMPLES::

        sage: (x^2).integral(x)
        1/3*x^3
        sage: (4*x*log(x)).integral(x)
        2*x^2*log(x) - x^2
        sage: cos(x).integral(x, algorithm='sympy')
        sin(x)
        sage: exp(x).integral(2*x)
        2*e^x

    """
    def __init__(self, *args, **kwds):
        """
        Init method for symbolic integration

        EXAMPLES::

            sage: log(x).integral(x)
            x*log(x) - x

        """
        kwds['built_in_function'] = True
        # Ask pynac not to apply chain rule
        kwds['apply_chain_rule'] = False
        for name in sfunctions_funcs:
            if hasattr(self, "_%s_"%name):
                kwds['%s_func'%name] = getattr(self, "_%s_"%name)

        SFunction.__init__(self, "integrate", *args, **kwds)

    def __call__(self, f, x, a=None, b=None, **kwds):
        """
        Call method of symbolic integration

        EXAMPLES::

            sage: log(x).integral(x)
            x*log(x) - x
            sage: cos(x).integral(x, algorithm='maxima')
            sin(x)
            sage: exp(x).integral(x,0,1)
            e - 1

        """
        # Check algorithm
        if kwds.has_key('algorithm'):
            if kwds['algorithm'] not in algorithm_table.keys():
                raise ValueError, "Unknown algorithm: %s" % kwds['algorithm']
            self._algorithm = kwds['algorithm']
        else:
            self._algorithm = None
        # Check arguments
        if isinstance(x, tuple):
            if len(x) == 3 and a is None and b is None:
                a = x[1]
                b = x[2]
                x = x[0]
            else:
                raise TypeError, "Limits of integration should be " \
                    "specified in the form: (x, a, b)"
        if a is not None:
            if b is not None:
                args = (f,x,a,b)
            else:
                raise TypeError, 'Only one endpoint is given'
        else:
            args = (f,x)

        return SFunction.__call__(self, *args)

    def _eval_(self, f, x, a=None, b=None, **kwds):
        """
        Returns the results of symbolic evaluation of the integral

        EXAMPLES::

            sage: exp(x).integral(x)
            e^x
            sage: exp(x).integral(x^2)
            2*(x - 1)*e^x
            sage: exp(x).integral(x,0,1)
            e - 1

        """
        # Check for x
        if not is_SymbolicVariable(x):
            if len(x.variables()) == 1:
                nx = x.variables()[0]
                f = f*x.diff(nx)
                x = nx
            else:
                return None
        if a is None:
            args = (f,x)
        else:
            args = (f,x,a,b)
        if self._algorithm is not None:
            try:
                return algorithm_table[self._algorithm](*args, **kwds)
            except NotImplementedError:
                return None
        # Default behavior is to try all listed integration algorithm
        for integrator in default_integrators:
            try:
                return integrator(*args, **kwds)
            except NotImplementedError:
                pass
        return None

    def _evalf_(self, f, x, a=None, b=None, **kwds):
        """
        Returns numerical approximation of the integral

        EXAMPLES::

            sage: h = (sin(x)/x^2).integral(x, 1, 2); h
            integrate(sin(x)/x^2, x, 1, 2)
            sage: h.n()
            0.472399177268953
        """
        if a is None:
            raise TypeError, "For numerical evaluation you need " \
                "to specify limits of integration"
        from sage.gsl.integration import numerical_integral
        # Currently, tuple doesn't coerce to SR. So we return
        # the result alone and not with the errors
        return numerical_integral(f, a, b)[0]

    def _derivative_(self, f, x, a=None, b=None, **kwds):
        """
        Returns derivative of symbolic integration

        EXAMPLES::

            sage: f(x) = function('f', x); a,b=var('a,b')
            sage: h = f(x).integral(x)
            sage: h.diff(x)
            f(x)
            sage: h.diff(a)
            0
            sage: h = f(x).integral(x,a,b)
            sage: h.diff(x)
            0
            sage: h.diff(a)
            -f(a)
            sage: h.diff(b)
            f(b)
        """
        dvar = kwds['diff_param'] # variable of differentiation
        x_diff = x.diff(dvar)
        if a is None:
            if not x_diff.is_zero():
                return f*x_diff
            else:
                return f.diff(dvar).integral(x)
        else:
            if x_diff.is_zero(): # dummy variable *x* is not same as *dvar*
                ans = f.diff(dvar).integral(x,a,b)
            else:
                ans = SR(0)
            return ans + f.subs(x==b)*b.diff(dvar) \
                        - f.subs(x==a)*a.diff(dvar)

    def _print_latex_(self, f, x, a=None, b=None, **kwds):
        r"""
        Returns LaTeX expression for integration of a symbolic function.

        EXAMPLES::

            sage: from sage.symbolic.integration.integration import integral
            sage: _ilatex = integral._print_latex_
            sage: var('x,a,b')
            (x, a, b)
            sage: f(x) = function('f',x)
            sage: _ilatex(f(x),x)
            '\\int f\\left(x\\right)\\,{d x}'
            sage: _ilatex(f(x),x,a,b)
            '\\int_{a}^{b} f\\left(x\\right)\\,{d x}'
            sage: latex(integrate(1/(1+sqrt(x)),x,0,1))
            \int_{0}^{1} \frac{1}{\sqrt{x} + 1}\,{d x}
        """
        from sage.misc.latex import latex
        if not is_SymbolicVariable(x):
            dx_str = "{d \\left(%s\\right)}"%(latex(x))
        else:
            dx_str = "{d %s}"%(latex(x))
        # Check whether its a definite integral
        if a is not None:
            return "\\int_{%s}^{%s} %s\\,%s"%(latex(a), latex(b), latex(f), dx_str)
        # Typeset as indefinite integral
        return "\\int %s\\,%s"%(latex(f), dx_str)

integral = SymbolicIntegration()
