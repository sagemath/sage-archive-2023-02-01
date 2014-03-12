r"""
Bessel Functions

This module provides symbolic Bessel Functions. These functions use the
`mpmath library`_ for numerical evaluation and Maxima, GiNaC, Pynac for
symbolics.

The main objects which are exported from this module are:

 * ``bessel_J`` -- The Bessel J function
 * ``bessel_Y`` -- The Bessel Y function
 * ``bessel_I`` -- The Bessel I function
 * ``bessel_K`` -- The Bessel K function
 * ``Bessel``   -- A factory function for producing Bessel functions of
   various kinds and orders

-  Bessel functions, first defined by the Swiss mathematician
   Daniel Bernoulli and named after Friedrich Bessel, are canonical
   solutions y(x) of Bessel's differential equation:

   .. math::

         x^2 \frac{d^2 y}{dx^2} + x \frac{dy}{dx} + \left(x^2 - \nu^2\right)y =
         0,

   for an arbitrary complex number `\nu` (the order).

-  In this module, `J_\nu` denotes the unique solution of Bessel's equation
   which is non-singular at `x = 0`. This function is known as the Bessel
   Function of the First Kind. This function also arises as a special case
   of the hypergeometric function `{}_0F_1`:

   .. math::

        J_\nu(x) = \frac{x^n}{2^\nu \Gamma(\nu + 1)} {}_0F_1(\nu +
        1, -\frac{x^2}{4}).

-  The second linearly independent solution to Bessel's equation (which is
   singular at `x=0`) is denoted by `Y_\nu` and is called the Bessel
   Function of the Second Kind:

   .. math::

        Y_\nu(x) = \frac{ J_\nu(x) \cos(\pi \nu) -
        J_{-\nu}(x)}{\sin(\pi \nu)}.

-  There are also two commonly used combinations of the Bessel J and Y
   Functions. The Bessel I Function, or the Modified Bessel Function of the
   First Kind, is defined by:

   .. math::

       I_\nu(x) = i^{-\nu} J_\nu(ix).

   The Bessel K Function, or the Modified Bessel Function of the Second Kind,
   is defined by:

   .. math::

       K_\nu(x) = \frac{\pi}{2} \cdot \frac{I_{-\nu}(x) -
       I_n(x)}{\sin(\pi \nu)}.

   We should note here that the above formulas for Bessel Y and K functions
   should be understood as limits when `\nu` is an integer.

-  It follows from Bessel's differential equation that the derivative of
   `J_n(x)` with respect to `x` is:

   .. math::

       \frac{d}{dx} J_n(x) = \frac{1}{x^n} \left(x^n J_{n-1}(x) - n x^{n-1}
       J_n(z) \right)

-  Another important formulation of the two linearly independent
   solutions to Bessel's equation are the Hankel functions
   `H_\nu^{(1)}(x)` and `H_\nu^{(2)}(x)`,
   defined by:

   .. math::

         H_\nu^{(1)}(x) = J_\nu(x) + i Y_\nu(x)

   .. math::

         H_\nu^{(2)}(x) = J_\nu(x) - i Y_\nu(x)

   where `i` is the imaginary unit (and `J_*` and
   `Y_*` are the usual J- and Y-Bessel functions). These
   linear combinations are also known as Bessel functions of the third
   kind; they are also two linearly independent solutions of Bessel's
   differential equation. They are named for Hermann Hankel.

EXAMPLES:

    Evaluate the Bessel J function symbolically and numerically::

        sage: bessel_J(0, x)
        bessel_J(0, x)
        sage: bessel_J(0, 0)
        bessel_J(0, 0)
        sage: bessel_J(0, x).diff(x)
        -1/2*bessel_J(1, x) + 1/2*bessel_J(-1, x)

        sage: N(bessel_J(0, 0), digits = 20)
        1.0000000000000000000
        sage: find_root(bessel_J(0,x), 0, 5)
        2.404825557695773

    Plot the Bessel J function::

        sage: f(x) = Bessel(0)(x); f
        x |--> bessel_J(0, x)
        sage: plot(f, (x, 1, 10))

    Visualize the Bessel Y function on the complex plane::

        sage: complex_plot(bessel_Y(0, x), (-5, 5), (-5, 5))

    Evaluate a combination of Bessel functions::

        sage: f(x) = bessel_J(1, x) - bessel_Y(0, x)
        sage: f(pi)
        bessel_J(1, pi) - bessel_Y(0, pi)
        sage: f(pi).n()
        -0.0437509653365599
        sage: f(pi).n(digits=50)
        -0.043750965336559909054985168023342675387737118378169

    Symbolically solve a second order differential equation with initial
    conditions `y(1) = a` and `y'(1) = b` in terms of Bessel functions::

        sage: y = function('y', x)
        sage: a, b = var('a, b')
        sage: diffeq = x^2*diff(y,x,x) + x*diff(y,x) + x^2*y == 0
        sage: f = desolve(diffeq, y, [1, a, b]); f
        (a*bessel_Y(1, 1) + b*bessel_Y(0, 1))*bessel_J(0, x)/(bessel_J(0,
        1)*bessel_Y(1, 1) - bessel_J(1, 1)*bessel_Y(0, 1)) -
        (a*bessel_J(1, 1) + b*bessel_J(0, 1))*bessel_Y(0, x)/(bessel_J(0,
        1)*bessel_Y(1, 1) - bessel_J(1, 1)*bessel_Y(0, 1))


    For more examples, see the docstring for :meth:`Bessel`.

AUTHORS:

    - Benjamin Jones (2012-12-27): initial version

    - Some of the documentation here has been adapted from David Joyner's
      original documentation of Sage's special functions module (2006).

REFERENCES:

    - Abramowitz and Stegun: Handbook of Mathematical Functions,
      http://www.math.sfu.ca/~cbm/aands/

    - http://en.wikipedia.org/wiki/Bessel_function

    - mpmath Library `Bessel Functions`_

.. _`mpmath Library`: http://code.google.com/p/mpmath/
.. _`Bessel Functions`: http://mpmath.googlecode.com/svn/trunk/doc/build/functions/bessel.html

"""

#*****************************************************************************
#       Copyright (C) 2013 Benjamin Jones <benjaminfjones@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.functions.other import sqrt
from sage.functions.log import exp
from sage.functions.hyperbolic import sinh, cosh
from sage.libs.mpmath import utils as mpmath_utils
from sage.misc.latex import latex
from sage.rings.all import RR, Integer
from sage.structure.coerce import parent
from sage.structure.element import get_coercion_model
from sage.symbolic.constants import pi
from sage.symbolic.function import BuiltinFunction, is_inexact
from sage.symbolic.expression import Expression

# remove after deprecation period
from sage.calculus.calculus import maxima
from sage.functions.other import real, imag
from sage.misc.sage_eval import sage_eval
from sage.rings.real_mpfr import RealField
from sage.plot.plot import plot
from sage.rings.all import ZZ


class Function_Bessel_J(BuiltinFunction):
    r"""
    The Bessel J Function, denoted by bessel_J(`\nu`, x) or `J_\nu(x)`.
    As a Taylor series about `x=0` it is equal to:

    .. math::

        J_\nu(x) = \sum_{k=0}^\infty \frac{(-1)^k}{k! \Gamma(k+\nu+1)}
        \left(\frac{x}{2}\right)^{2k+\nu}

    The parameter `\nu` is called the order and may be any real or
    complex number; however, integer and half-integer values are most
    common. It is defined for all complex numbers `x` when `\nu`
    is an integer or greater than zero and it diverges as `x \to 0`
    for negative non-integer values of `\nu`.

    For integer orders `\nu = n` there is an integral representation:

    .. math::

        J_n(x) = \frac{1}{\pi} \int_0^\pi \cos(n t - x \sin(t)) \; dt

    This function also arises as a special case of the hypergeometric
    function `{}_0F_1`:

    .. math::

        J_\nu(x) = \frac{x^n}{2^\nu \Gamma(\nu + 1)} {}_0F_1\left(\nu +
        1, -\frac{x^2}{4}\right).

    EXAMPLES::

        sage: bessel_J(1.0, 1.0)
        0.440050585744933
        sage: bessel_J(2, I).n(digits=30)
        -0.135747669767038281182852569995

        sage: bessel_J(1, x)
        bessel_J(1, x)
        sage: n = var('n')
        sage: bessel_J(n, x)
        bessel_J(n, x)

    Examples of symbolic manipulation::

        sage: a = bessel_J(pi, bessel_J(1, I)); a
        bessel_J(pi, bessel_J(1, I))
        sage: N(a, digits=20)
        0.00059023706363796717363 - 0.0026098820470081958110*I

        sage: f = bessel_J(2, x)
        sage: f.diff(x)
        -1/2*bessel_J(3, x) + 1/2*bessel_J(1, x)

    Comparison to a well-known integral representation of `J_1(1)`::

        sage: A = numerical_integral(1/pi*cos(x - sin(x)), 0, pi)
        sage: A[0]  # abs tol 1e-14
        0.44005058574493355
        sage: bessel_J(1.0, 1.0) - A[0] < 1e-15
        True

    Currently, integration is not supported (directly) since we cannot
    yet convert hypergeometric functions to and from Maxima::

        sage: f = bessel_J(2, x)
        sage: f.integrate(x)
        Traceback (most recent call last):
        ...
        TypeError: cannot coerce arguments: no canonical coercion from <type 'list'> to Symbolic Ring

        sage: m = maxima(bessel_J(2, x))
        sage: m.integrate(x)
        hypergeometric([3/2],[5/2,3],-x^2/4)*x^3/24

    Visualization::

        sage: plot(bessel_J(1,x), (x,0,5), color='blue')
        sage: complex_plot(bessel_J(1, x), (-5, 5), (-5, 5)) # long time

    ALGORITHM:

        Numerical evaluation is handled by the mpmath library. Symbolics are
        handled by a combination of Maxima and Sage (Ginac/Pynac).
    """
    def __init__(self):
        """
        See the docstring for :meth:`Function_Bessel_J`.

        EXAMPLES::

            sage: sage.functions.bessel.Function_Bessel_J()
            bessel_J
        """
        BuiltinFunction.__init__(self, "bessel_J", nargs=2,
                                 conversions=dict(mathematica='BesselJ',
                                                  maxima='bessel_j',
                                                  sympy='besselj'))

    # remove after deprecation period
    def __call__(self, *args, **kwds):
        """
        Custom ``__call__`` method which uses the old Bessel function code if
        the ``algorithm`` or ``prec`` arguments are used. This should be
        removed after the deprecation period.

        EXAMPLES::

            sage: bessel_J(0, 1.0, "maxima", 53)
            doctest:1: DeprecationWarning: precision argument is deprecated; algorithm argument is currently deprecated, but will be available as a named keyword in the future
            See http://trac.sagemath.org/4102 for details.
            .7651976865579666
        """
        if len(args) > 2 or len(kwds) > 0:
            from sage.misc.superseded import deprecation
            deprecation(4102, 'precision argument is deprecated; algorithm '
                              'argument is currently deprecated, but will be '
                              'available as a named keyword in the future')
            return _bessel_J(*args, **kwds)
        else:
            return super(BuiltinFunction, self).__call__(*args, **kwds)

    def _eval_(self, n, x):
        """
        EXAMPLES::

            sage: a, b = var('a, b')
            sage: bessel_J(a, b)
            bessel_J(a, b)
            sage: bessel_J(1.0, 1.0)
            0.440050585744933
        """
        if (not isinstance(n, Expression) and
                not isinstance(x, Expression) and
                (is_inexact(n) or is_inexact(x))):
            coercion_model = get_coercion_model()
            n, x = coercion_model.canonical_coercion(n, x)
            return self._evalf_(n, x, parent(n))

        return None

    def _evalf_(self, n, x, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: bessel_J(0.0, 1.0)
            0.765197686557966
            sage: bessel_J(0, 1).n(digits=20)
            0.76519768655796655145
        """
        import mpmath
        return mpmath_utils.call(mpmath.besselj, n, x, parent=parent)

    def _derivative_(self, n, x, diff_param):
        """
        Return the derivative of the Bessel J function.

        EXAMPLES::

            sage: f(z) = bessel_J(10, z)
            sage: derivative(f, z)
            z |--> -1/2*bessel_J(11, z) + 1/2*bessel_J(9, z)
            sage: nu = var('nu')
            sage: bessel_J(nu, z).diff(nu)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative with respect to order

        """
        if diff_param == 1:
            return (bessel_J(n - 1, x) - bessel_J(n + 1, x)) / Integer(2)
        else:
            raise NotImplementedError('derivative with respect to order')

    def _print_latex_(self, n, z):
        """
        Custom _print_latex_ method.

        EXAMPLES::

            sage: latex(bessel_J(1, x))
            \operatorname{J_{1}}(x)
        """
        return r"\operatorname{J_{%s}}(%s)" % (latex(n), latex(z))

bessel_J = Function_Bessel_J()


class Function_Bessel_Y(BuiltinFunction):
    r"""
    The Bessel Y functions, also known as the Bessel functions of the second
    kind, Weber functions, or Neumann functions.

    `Y_\nu(z)` is a holomorphic function of `z` on the complex plane,
    cut along the negative real axis. It is singular at `z = 0`. When `z`
    is fixed, `Y_\nu(z)` is an entire function of the order `\nu`.

    DEFINITION:

    .. math::

        Y_n(z) = \frac{J_\nu(z) \cos(\nu z) -
        J_{-\nu}(z)}{\sin(\nu z)}

    Its derivative with respect to `z` is:

    .. math::

        \frac{d}{dz} Y_n(z) = \frac{1}{z^n} \left(z^n Y_{n-1}(z) - n z^{n-1}
        Y_n(z) \right)

    EXAMPLES::

        sage: bessel_Y(1, x)
        bessel_Y(1, x)
        sage: bessel_Y(1.0, 1.0)
        -0.781212821300289
        sage: n = var('n')
        sage: bessel_Y(n, x)
        bessel_Y(n, x)
        sage: bessel_Y(2, I).n()
        1.03440456978312 - 0.135747669767038*I
        sage: bessel_Y(0, 0).n()
        -infinity

    Examples of symbolic manipulation::

        sage: a = bessel_Y(pi, bessel_Y(1, I)); a
        bessel_Y(pi, bessel_Y(1, I))
        sage: N(a, digits=20)
        4.2059146571791095708 + 21.307914215321993526*I

        sage: f = bessel_Y(2, x)
        sage: f.diff(x)
        -1/2*bessel_Y(3, x) + 1/2*bessel_Y(1, x)

    High precision and complex valued inputs (see :trac:`4230`)::

        sage: bessel_Y(0, 1).n(128)
        0.088256964215676957982926766023515162828
        sage: bessel_Y(0, RealField(200)(1))
        0.088256964215676957982926766023515162827817523090675546711044
        sage: bessel_Y(0, ComplexField(200)(0.5+I))
        0.077763160184438051408593468823822434235010300228009867784073 + 1.0142336049916069152644677682828326441579314239591288411739*I

    Visualization::

        sage: plot(bessel_Y(1,x), (x,0,5), color='blue')
        sage: complex_plot(bessel_Y(1, x), (-5, 5), (-5, 5)) # long time

    ALGORITHM:

        Numerical evaluation is handled by the mpmath library. Symbolics are
        handled by a combination of Maxima and Sage (Ginac/Pynac).
    """
    def __init__(self):
        """
        See the docstring for :meth:`Function_Bessel_Y`.

        EXAMPLES::

            sage: sage.functions.bessel.Function_Bessel_Y()(0, x)
            bessel_Y(0, x)
        """
        BuiltinFunction.__init__(self, "bessel_Y", nargs=2,
                                 conversions=dict(mathematica='BesselY',
                                                  maxima='bessel_y',
                                                  sympy='bessely'))

    # remove after deprecation period
    def __call__(self, *args, **kwds):
        """
        Custom ``__call__`` method which uses the old Bessel function code if
        the ``algorithm`` or ``prec`` arguments are used. This should be
        removed after the deprecation period.

        EXAMPLES::

            sage: bessel_Y(0, 1, "maxima", 53)
            doctest:1: DeprecationWarning: precision argument is deprecated; algorithm argument is currently deprecated, but will be available as a named keyword in the future
            See http://trac.sagemath.org/4102 for details.
            0.0882569642156769
        """
        if len(args) > 2 or len(kwds) > 0:
            from sage.misc.superseded import deprecation
            deprecation(4102, 'precision argument is deprecated; algorithm '
                              'argument is currently deprecated, but will be '
                              'available as a named keyword in the future')
            return _bessel_Y(*args, **kwds)
        else:
            return super(BuiltinFunction, self).__call__(*args, **kwds)

    def _eval_(self, n, x):
        """
        EXAMPLES::

            sage: a,b = var('a, b')
            sage: bessel_Y(a, b)
            bessel_Y(a, b)
            sage: bessel_Y(0, 1).n(128)
            0.088256964215676957982926766023515162828
        """
        if (not isinstance(n, Expression) and not isinstance(x, Expression) and
                (is_inexact(n) or is_inexact(x))):
            coercion_model = get_coercion_model()
            n, x = coercion_model.canonical_coercion(n, x)
            return self._evalf_(n, x, parent(n))

        return None  # leaves the expression unevaluated

    def _evalf_(self, n, x, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: bessel_Y(1.0+2*I, 3.0+4*I)
            0.699410324467538 + 0.228917940896421*I
            sage: bessel_Y(0, 1).n(256)
            0.08825696421567695798292676602351516282781752309067554671104384761199978932351
        """
        import mpmath
        return mpmath_utils.call(mpmath.bessely, n, x, parent=parent)

    def _derivative_(self, n, x, diff_param):
        """
        Return the derivative of the Bessel Y function.

        EXAMPLES::

            sage: f(x) = bessel_Y(10, x)
            sage: derivative(f, x)
            x |--> -1/2*bessel_Y(11, x) + 1/2*bessel_Y(9, x)
            sage: nu = var('nu')
            sage: bessel_Y(nu, x).diff(nu)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative with respect to order
        """
        if diff_param == 1:
            return (bessel_Y(n - 1, x) - bessel_Y(n + 1, x)) / Integer(2)
        else:
            raise NotImplementedError('derivative with respect to order')

    def _print_latex_(self, n, z):
        """
        Custom _print_latex_ method.

        EXAMPLES::

            sage: latex(bessel_Y(1, x))
            \operatorname{Y_{1}}(x)
        """
        return r"\operatorname{Y_{%s}}(%s)" % (latex(n), latex(z))

bessel_Y = Function_Bessel_Y()


class Function_Bessel_I(BuiltinFunction):
    r"""
    The Bessel I function, or the Modified Bessel Function of the First Kind.

    DEFINITION:

    .. math::

        I_\nu(x) = i^{-\nu} J_\nu(ix)

    EXAMPLES::

        sage: bessel_I(1, x)
        bessel_I(1, x)
        sage: bessel_I(1.0, 1.0)
        0.565159103992485
        sage: n = var('n')
        sage: bessel_I(n, x)
        bessel_I(n, x)
        sage: bessel_I(2, I).n()
        -0.114903484931900

    Examples of symbolic manipulation::

        sage: a = bessel_I(pi, bessel_I(1, I))
        sage: N(a, digits=20)
        0.00026073272117205890528 - 0.0011528954889080572266*I

        sage: f = bessel_I(2, x)
        sage: f.diff(x)
        1/2*bessel_I(3, x) + 1/2*bessel_I(1, x)

    Special identities that bessel_I satisfies::

        sage: bessel_I(1/2, x)
        sqrt(2)*sqrt(1/(pi*x))*sinh(x)
        sage: eq = bessel_I(1/2, x) == bessel_I(0.5, x)
        sage: eq.test_relation()
        True
        sage: bessel_I(-1/2, x)
        sqrt(2)*sqrt(1/(pi*x))*cosh(x)
        sage: eq = bessel_I(-1/2, x) == bessel_I(-0.5, x)
        sage: eq.test_relation()
        True

    Examples of asymptotic behavior::

        sage: limit(bessel_I(0, x), x=oo)
        +Infinity
        sage: limit(bessel_I(0, x), x=0)
        1

    High precision and complex valued inputs::

        sage: bessel_I(0, 1).n(128)
        1.2660658777520083355982446252147175376
        sage: bessel_I(0, RealField(200)(1))
        1.2660658777520083355982446252147175376076703113549622068081
        sage: bessel_I(0, ComplexField(200)(0.5+I))
        0.80644357583493619472428518415019222845373366024179916785502 + 0.22686958987911161141397453401487525043310874687430711021434*I

    Visualization::

        sage: plot(bessel_I(1,x), (x,0,5), color='blue')
        sage: complex_plot(bessel_I(1, x), (-5, 5), (-5, 5)) # long time

    ALGORITHM:

        Numerical evaluation is handled by the mpmath library. Symbolics are
        handled by a combination of Maxima and Sage (Ginac/Pynac).
    """
    def __init__(self):
        """
        See the docstring for :meth:`Function_Bessel_I`.

        EXAMPLES::

            sage: bessel_I(1,x)
            bessel_I(1, x)
        """
        BuiltinFunction.__init__(self, "bessel_I", nargs=2,
                                 conversions=dict(mathematica='BesselI',
                                                  maxima='bessel_i',
                                                  sympy='besseli'))

    # remove after deprecation period
    def __call__(self, *args, **kwds):
        """
        Custom ``__call__`` method which uses the old Bessel function code if
        the ``algorithm`` or ``prec`` arguments are used. This should be
        removed after the deprecation period.

        EXAMPLES::

            sage: bessel_I(0, 1, "maxima", 53)
            doctest:1: DeprecationWarning: precision argument is deprecated; algorithm argument is currently deprecated, but will be available as a named keyword in the future
            See http://trac.sagemath.org/4102 for details.
            1.266065877752009
        """
        if len(args) > 2 or len(kwds) > 0:
            from sage.misc.superseded import deprecation
            deprecation(4102, 'precision argument is deprecated; algorithm '
                              'argument is currently deprecated, but will be '
                              'available as a named keyword in the future')
            return _bessel_I(*args, **kwds)
        else:
            return super(BuiltinFunction, self).__call__(*args, **kwds)

    def _eval_(self, n, x):
        """
        EXAMPLES::

            sage: y=var('y')
            sage: bessel_I(y,x)
            bessel_I(y, x)
            sage: bessel_I(0.0, 1.0)
            1.26606587775201
            sage: bessel_I(1/2, 1)
            sqrt(2)*sinh(1)/sqrt(pi)
            sage: bessel_I(-1/2, pi)
            sqrt(2)*cosh(pi)/pi
        """
        if (not isinstance(n, Expression) and not isinstance(x, Expression) and
                (is_inexact(n) or is_inexact(x))):
            coercion_model = get_coercion_model()
            n, x = coercion_model.canonical_coercion(n, x)
            return self._evalf_(n, x, parent(n))

        # special identities
        if n == Integer(1) / Integer(2):
            return sqrt(2 / (pi * x)) * sinh(x)
        elif n == -Integer(1) / Integer(2):
            return sqrt(2 / (pi * x)) * cosh(x)

        return None  # leaves the expression unevaluated

    def _evalf_(self, n, x, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: bessel_I(1,3).n(digits=20)
            3.9533702174026093965
        """
        import mpmath
        return mpmath_utils.call(mpmath.besseli, n, x, parent=parent)

    def _derivative_(self, n, x, diff_param):
        """
        Return the derivative of the Bessel I function `I_n(x)` with respect
        to `x`.

        EXAMPLES::

            sage: f(z) = bessel_I(10, x)
            sage: derivative(f, x)
            z |--> 1/2*bessel_I(11, x) + 1/2*bessel_I(9, x)
            sage: nu = var('nu')
            sage: bessel_I(nu, x).diff(nu)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative with respect to order
        """
        if diff_param == 1:
            return (bessel_I(n - 1, x) + bessel_I(n + 1, x)) / Integer(2)
        else:
            raise NotImplementedError('derivative with respect to order')

    def _print_latex_(self, n, z):
        """
        Custom _print_latex_ method.

        EXAMPLES::

            sage: latex(bessel_I(1, x))
            \operatorname{I_{1}}(x)
        """
        return r"\operatorname{I_{%s}}(%s)" % (latex(n), latex(z))

bessel_I = Function_Bessel_I()


class Function_Bessel_K(BuiltinFunction):
    r"""
    The Bessel K function, or the modified Bessel function of the second kind.

    DEFINITION:

    .. math::

        K_\nu(x) = \frac{\pi}{2} \frac{I_{-\nu}(x)-I_\nu(x)}{\sin(\nu \pi)}

    EXAMPLES::

        sage: bessel_K(1, x)
        bessel_K(1, x)
        sage: bessel_K(1.0, 1.0)
        0.601907230197235
        sage: n = var('n')
        sage: bessel_K(n, x)
        bessel_K(n, x)
        sage: bessel_K(2, I).n()
        -2.59288617549120 + 0.180489972066962*I

    Examples of symbolic manipulation::

        sage: a = bessel_K(pi, bessel_K(1, I)); a
        bessel_K(pi, bessel_K(1, I))
        sage: N(a, digits=20)
        3.8507583115005220157 + 0.068528298579883425792*I

        sage: f = bessel_K(2, x)
        sage: f.diff(x)
        1/2*bessel_K(3, x) + 1/2*bessel_K(1, x)

        sage: bessel_K(1/2, x)
        bessel_K(1/2, x)
        sage: bessel_K(1/2, -1)
        bessel_K(1/2, -1)
        sage: bessel_K(1/2, 1)
        sqrt(1/2)*sqrt(pi)*e^(-1)

    Examples of asymptotic behavior::

        sage: bessel_K(0, 0.0)
        +infinity
        sage: limit(bessel_K(0, x), x=0)
        +Infinity
        sage: limit(bessel_K(0, x), x=oo)
        0

    High precision and complex valued inputs::

        sage: bessel_K(0, 1).n(128)
        0.42102443824070833333562737921260903614
        sage: bessel_K(0, RealField(200)(1))
        0.42102443824070833333562737921260903613621974822666047229897
        sage: bessel_K(0, ComplexField(200)(0.5+I))
        0.058365979093103864080375311643360048144715516692187818271179 - 0.67645499731334483535184142196073004335768129348518210260256*I

    Visualization::

        sage: plot(bessel_K(1,x), (x,0,5), color='blue')
        sage: complex_plot(bessel_K(1, x), (-5, 5), (-5, 5)) # long time

    ALGORITHM:

        Numerical evaluation is handled by the mpmath library. Symbolics are
        handled by a combination of Maxima and Sage (Ginac/Pynac).

    TESTS:

    Verify that :trac:`3426` is fixed:

    The Bessel K function can be evaluated numerically at complex orders::

        sage: bessel_K(10 * I, 10).n()
        9.82415743819925e-8

    For a fixed imaginary order and increasing, real, second component the
    value of Bessel K is exponentially decaying::

        sage: for x in [10, 20, 50, 100, 200]: print bessel_K(5*I, x).n()
        5.27812176514912e-6
        3.11005908421801e-10
        2.66182488515423e-23 - 8.59622057747552e-58*I
        4.11189776828337e-45 - 1.01494840019482e-80*I
        1.15159692553603e-88 - 6.75787862113718e-125*I
    """
    def __init__(self):
        """
        See the docstring for :meth:`Function_Bessel_K`.

        EXAMPLES::

            sage: sage.functions.bessel.Function_Bessel_K()
            bessel_K
        """
        BuiltinFunction.__init__(self, "bessel_K", nargs=2,
                                 conversions=dict(mathematica='BesselK',
                                                  maxima='bessel_k',
                                                  sympy='besselk'))

    # remove after deprecation period
    def __call__(self, *args, **kwds):
        """
        Custom ``__call__`` method which uses the old Bessel function code if
        the ``algorithm`` or ``prec`` arguments are used. This should be
        removed after the deprecation period.

        EXAMPLES::

            sage: bessel_K(0, 1, "maxima", 53)
            doctest:1: DeprecationWarning: precision argument is deprecated; algorithm argument is currently deprecated, but will be available as a named keyword in the future
            See http://trac.sagemath.org/4102 for details.
            0.0882569642156769
        """
        if len(args) > 2 or len(kwds) > 0:
            from sage.misc.superseded import deprecation
            deprecation(4102, 'precision argument is deprecated; algorithm '
                              'argument is currently deprecated, but will be '
                              'available as a named keyword in the future')
            return _bessel_Y(*args, **kwds)
        else:
            return super(BuiltinFunction, self).__call__(*args, **kwds)

    def _eval_(self, n, x):
        """
        EXAMPLES::

            sage: bessel_K(1,0)
            bessel_K(1, 0)
            sage: bessel_K(1.0, 0.0)
            +infinity
            sage: bessel_K(-1, 1).n(128)
            0.60190723019723457473754000153561733926
        """
        if (not isinstance(n, Expression) and not isinstance(x, Expression) and
                (is_inexact(n) or is_inexact(x))):
            coercion_model = get_coercion_model()
            n, x = coercion_model.canonical_coercion(n, x)
            return self._evalf_(n, x, parent(n))

        # special identity
        if n == Integer(1) / Integer(2) and x > 0:
            return sqrt(pi / 2) * exp(-x) * x ** (-Integer(1) / Integer(2))

        return None  # leaves the expression unevaluated

    def _evalf_(self, n, x, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: bessel_K(0.0, 1.0)
            0.421024438240708
            sage: bessel_K(0, RealField(128)(1))
            0.42102443824070833333562737921260903614
        """
        import mpmath
        return mpmath_utils.call(mpmath.besselk, n, x, parent=parent)

    def _derivative_(self, n, x, diff_param):
        """
        Return the derivative of the Bessel K function.

        EXAMPLES::

            sage: f(x) = bessel_K(10, x)
            sage: derivative(f, x)
            x |--> 1/2*bessel_K(11, x) + 1/2*bessel_K(9, x)
            sage: nu = var('nu')
            sage: bessel_K(nu, x).diff(nu)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative with respect to order
        """
        if diff_param == 1:
            return (bessel_K(n - 1, x) + bessel_K(n + 1, x)) / Integer(2)
        else:
            raise NotImplementedError('derivative with respect to order')

    def _print_latex_(self, n, z):
        """
        Custom _print_latex_ method.

        EXAMPLES::

            sage: latex(bessel_K(1, x))
            \operatorname{K_{1}}(x)
        """
        return r"\operatorname{K_{%s}}(%s)" % (latex(n), latex(z))

bessel_K = Function_Bessel_K()


# dictionary used in Bessel
bessel_type_dict = {'I': bessel_I, 'J': bessel_J, 'K': bessel_K, 'Y': bessel_Y}


def Bessel(*args, **kwds):
    """
    A function factory that produces symbolic I, J, K, and Y Bessel functions.
    There are several ways to call this function:

        - ``Bessel(order, type)``
        - ``Bessel(order)`` -- type defaults to 'J'
        - ``Bessel(order, typ=T)``
        - ``Bessel(typ=T)`` -- order is unspecified, this is a 2-parameter
          function
        - ``Bessel()`` -- order is unspecified, type is 'J'

    where ``order`` can be any integer and T must be one of the strings 'I',
    'J', 'K', or 'Y'.

    See the EXAMPLES below.

    EXAMPLES:

    Construction of Bessel functions with various orders and types::

        sage: Bessel()
        bessel_J
        sage: Bessel(1)(x)
        bessel_J(1, x)
        sage: Bessel(1, 'Y')(x)
        bessel_Y(1, x)
        sage: Bessel(-2, 'Y')(x)
        bessel_Y(-2, x)
        sage: Bessel(typ='K')
        bessel_K
        sage: Bessel(0, typ='I')(x)
        bessel_I(0, x)

    Evaluation::

        sage: f = Bessel(1)
        sage: f(3.0)
        0.339058958525936
        sage: f(3)
        bessel_J(1, 3)
        sage: f(3).n(digits=50)
        0.33905895852593645892551459720647889697308041819801

        sage: g = Bessel(typ='J')
        sage: g(1,3)
        bessel_J(1, 3)
        sage: g(2, 3+I).n()
        0.634160370148554 + 0.0253384000032695*I
        sage: abs(numerical_integral(1/pi*cos(3*sin(x)), 0.0, pi)[0] - Bessel(0, 'J')(3.0)) < 1e-15
        True

    Symbolic calculus::

        sage: f(x) = Bessel(0, 'J')(x)
        sage: derivative(f, x)
        x |--> -1/2*bessel_J(1, x) + 1/2*bessel_J(-1, x)
        sage: derivative(f, x, x)
        x |--> 1/4*bessel_J(2, x) - 1/2*bessel_J(0, x) + 1/4*bessel_J(-2, x)

    Verify that `J_0` satisfies Bessel's differential equation numerically
    using the ``test_relation()`` method::

        sage: y = bessel_J(0, x)
        sage: diffeq = x^2*derivative(y,x,x) + x*derivative(y,x) + x^2*y == 0
        sage: diffeq.test_relation(proof=False)
        True

    Conversion to other systems::

        sage: x,y = var('x,y')
        sage: f = maxima(Bessel(typ='K')(x,y))
        sage: f.derivative('x')
        %pi*csc(%pi*x)*('diff(bessel_i(-x,y),x,1)-'diff(bessel_i(x,y),x,1))/2-%pi*bessel_k(x,y)*cot(%pi*x)
        sage: f.derivative('y')
        -(bessel_k(x+1,y)+bessel_k(x-1,y))/2

    Compute the particular solution to Bessel's Differential Equation that
    satisfies `y(1) = 1` and `y'(1) = 1`, then verify the initial conditions
    and plot it::

        sage: y = function('y', x)
        sage: diffeq = x^2*diff(y,x,x) + x*diff(y,x) + x^2*y == 0
        sage: f = desolve(diffeq, y, [1, 1, 1]); f
        (bessel_Y(1, 1) + bessel_Y(0, 1))*bessel_J(0, x)/(bessel_J(0,
        1)*bessel_Y(1, 1) - bessel_J(1, 1)*bessel_Y(0, 1)) - (bessel_J(1,
        1) + bessel_J(0, 1))*bessel_Y(0, x)/(bessel_J(0, 1)*bessel_Y(1, 1)
        - bessel_J(1, 1)*bessel_Y(0, 1))
        sage: f.subs(x=1).n() # numerical verification
        1.00000000000000
        sage: fp = f.diff(x)
        sage: fp.subs(x=1).n()
        1.00000000000000

        sage: f.subs(x=1).simplify_full() # symbolic verification
        1
        sage: fp = f.diff(x)
        sage: fp.subs(x=1).simplify_full()
        1

        sage: plot(f, (x,0,5))

    Plotting::

        sage: f(x) = Bessel(0)(x); f
        x |--> bessel_J(0, x)
        sage: plot(f, (x, 1, 10))

        sage: plot([ Bessel(i, 'J') for i in range(5) ], 2, 10)

        sage: G = Graphics()
        sage: G += sum([ plot(Bessel(i), 0, 4*pi, rgbcolor=hue(sin(pi*i/10))) for i in range(5) ])
        sage: show(G)

    A recreation of Abramowitz and Stegun Figure 9.1::

        sage: G  = plot(Bessel(0, 'J'), 0, 15, color='black')
        sage: G += plot(Bessel(0, 'Y'), 0, 15, color='black')
        sage: G += plot(Bessel(1, 'J'), 0, 15, color='black', linestyle='dotted')
        sage: G += plot(Bessel(1, 'Y'), 0, 15, color='black', linestyle='dotted')
        sage: show(G, ymin=-1, ymax=1)

    """
    # Determine the order and type of function from the arguments and keywords.
    # These are recored in local variables: _type, _order, _system, _nargs.
    _type = None
    if len(args) == 0:    # no order specified
        _order = None
        _nargs = 2
    elif len(args) == 1:  # order is specified
        _order = args[0]
        _nargs = 1
    elif len(args) == 2:  # both order and type are positional arguments
        _order = args[0]
        _type = args[1]
        _nargs = 1
    else:
        from sage.misc.superseded import deprecation
        deprecation(4102, 'precision argument is deprecated; algorithm '
                          'argument is currently deprecated, but will be '
                          'available as a named keyword in the future')
        return _Bessel(*args, **kwds)

    # check for type inconsistency
    if _type is not None and 'typ' in kwds and _type != kwds['typ']:
        raise ValueError("inconsistent types given")
    # record the function type
    if _type is None:
        if 'typ' in kwds:
            _type = kwds['typ']
        else:
            _type = 'J'
    if not (_type in ['I', 'J', 'K', 'Y']):
        raise ValueError("type must be one of I, J, K, Y")
    # record the numerical evaluation system
    if 'algorithm' in kwds:
        _system = kwds['algorithm']
    else:
        _system = 'mpmath'

    # return the function
    _f = bessel_type_dict[_type]
    if _nargs == 1:
        return lambda x: _f(_order, x)
    else:
        return _f

####################################################
###  to be removed after the deprecation period  ###
####################################################


def _bessel_I(nu,z,algorithm = "pari",prec=53):
    r"""
    Implements the "I-Bessel function", or "modified Bessel function,
    1st kind", with index (or "order") nu and argument z.

    INPUT:


    -  ``nu`` - a real (or complex, for pari) number

    -  ``z`` - a real (positive) algorithm - "pari" or
       "maxima" or "scipy" prec - real precision (for PARI only)


    DEFINITION::

            Maxima:
                             inf
                            ====   - nu - 2 k  nu + 2 k
                            \     2          z
                             >    -------------------
                            /     k! Gamma(nu + k + 1)
                            ====
                            k = 0

            PARI:

                             inf
                            ====   - 2 k  2 k
                            \     2      z    Gamma(nu + 1)
                             >    -----------------------
                            /       k! Gamma(nu + k + 1)
                            ====
                            k = 0



    Sometimes ``bessel_I(nu,z)`` is denoted
    ``I_nu(z)`` in the literature.

    .. warning::

       In Maxima (the manual says) i0 is deprecated but
       ``bessel_i(0,*)`` is broken. (Was fixed in recent CVS patch
       though.)

    EXAMPLES::

        sage: from sage.functions.bessel import _bessel_I
        sage: _bessel_I(1,1,"pari",500)
        0.565159103992485027207696027609863307328899621621092009480294489479255640964371134092664997766814410064677886055526302676857637684917179812041131208121
        sage: _bessel_I(1,1)
        0.565159103992485
        sage: _bessel_I(2,1.1,"maxima")
        0.16708949925104...
        sage: _bessel_I(0,1.1,"maxima")
        1.32616018371265...
        sage: _bessel_I(0,1,"maxima")
        1.2660658777520...
        sage: _bessel_I(1,1,"scipy")
        0.565159103992...

    Check whether the return value is real whenever the argument is real (#10251)::

        sage: _bessel_I(5, 1.5, algorithm='scipy') in RR
        True

    """
    if algorithm=="pari":
        from sage.libs.pari.all import pari
        try:
            R = RealField(prec)
            nu = R(nu)
            z = R(z)
        except TypeError:
            C = ComplexField(prec)
            nu = C(nu)
            z = C(z)
            K = C
        K = z.parent()
        return K(pari(nu).besseli(z, precision=prec))
    elif algorithm=="scipy":
        if prec != 53:
            raise ValueError, "for the scipy algorithm the precision must be 53"
        import scipy.special
        ans = str(scipy.special.iv(float(nu),complex(real(z),imag(z))))
        ans = ans.replace("(","")
        ans = ans.replace(")","")
        ans = ans.replace("j","*I")
        ans = sage_eval(ans)
        return real(ans) if z in RR else ans # Return real value when arg is real
    elif algorithm == "maxima":
        if prec != 53:
            raise ValueError, "for the maxima algorithm the precision must be 53"
        return sage_eval(maxima.eval("bessel_i(%s,%s)"%(float(nu),float(z))))
    else:
        raise ValueError, "unknown algorithm '%s'"%algorithm

def _bessel_J(nu,z,algorithm="pari",prec=53):
    r"""
    Return value of the "J-Bessel function", or "Bessel function, 1st
    kind", with index (or "order") nu and argument z.

    ::

            Defn:
            Maxima:
                             inf
                            ====          - nu - 2 k  nu + 2 k
                            \     (-1)^k 2           z
                             >    -------------------------
                            /        k! Gamma(nu + k + 1)
                            ====
                            k = 0

            PARI:

                             inf
                            ====          - 2k    2k
                            \     (-1)^k 2      z    Gamma(nu + 1)
                             >    ----------------------------
                            /         k! Gamma(nu + k + 1)
                            ====
                            k = 0


    Sometimes bessel_J(nu,z) is denoted J_nu(z) in the literature.

    .. warning::

       Inaccurate for small values of z.

    EXAMPLES::

        sage: from sage.functions.bessel import _bessel_J
        sage: _bessel_J(2,1.1)
        0.136564153956658
        sage: _bessel_J(0,1.1)
        0.719622018527511
        sage: _bessel_J(0,1)
        0.765197686557967
        sage: _bessel_J(0,0)
        1.00000000000000
        sage: _bessel_J(0.1,0.1)
        0.777264368097005

    We check consistency of PARI and Maxima::

        sage: n(_bessel_J(3,10,"maxima"))
        0.0583793793051...
        sage: n(_bessel_J(3,10,"pari"))
        0.0583793793051868
        sage: _bessel_J(3,10,"scipy")
        0.0583793793052...

    Check whether the return value is real whenever the argument is real (#10251)::
        sage: _bessel_J(5, 1.5, algorithm='scipy') in RR
        True
    """

    if algorithm=="pari":
        from sage.libs.pari.all import pari
        try:
            R = RealField(prec)
            nu = R(nu)
            z = R(z)
        except TypeError:
            C = ComplexField(prec)
            nu = C(nu)
            z = C(z)
            K = C
        if nu == 0:
            nu = ZZ(0)
        K = z.parent()
        return K(pari(nu).besselj(z, precision=prec))
    elif algorithm=="scipy":
        if prec != 53:
            raise ValueError, "for the scipy algorithm the precision must be 53"
        import scipy.special
        ans = str(scipy.special.jv(float(nu),complex(real(z),imag(z))))
        ans = ans.replace("(","")
        ans = ans.replace(")","")
        ans = ans.replace("j","*I")
        ans = sage_eval(ans)
        return real(ans) if z in RR else ans
    elif algorithm == "maxima":
        if prec != 53:
            raise ValueError, "for the maxima algorithm the precision must be 53"
        f = maxima.function('n,z', 'bessel_j(n, z)')
        return f(nu, z)
    else:
        raise ValueError, "unknown algorithm '%s'"%algorithm

def _bessel_K(nu,z,algorithm="pari",prec=53):
    r"""
    Implements the "K-Bessel function", or "modified Bessel function,
    2nd kind", with index (or "order") nu and argument z. Defn::

                    pi*(bessel_I(-nu, z) - bessel_I(nu, z))
                   ----------------------------------------
                                2*sin(pi*nu)


    if nu is not an integer and by taking a limit otherwise.

    Sometimes bessel_K(nu,z) is denoted K_nu(z) in the literature. In
    PARI, nu can be complex and z must be real and positive.

    EXAMPLES::

        sage: from sage.functions.bessel import _bessel_K
        sage: _bessel_K(3,2,"scipy")
        0.64738539094...
        sage: _bessel_K(3,2)
        0.64738539094...
        sage: _bessel_K(1,1)
        0.60190723019...
        sage: _bessel_K(1,1,"pari",10)
        0.60
        sage: _bessel_K(1,1,"pari",100)
        0.60190723019723457473754000154

    TESTS::

        sage: _bessel_K(2,1.1, algorithm="maxima")
        Traceback (most recent call last):
        ...
        NotImplementedError: The K-Bessel function is only implemented for the pari and scipy algorithms

        Check whether the return value is real whenever the argument is real (#10251)::

        sage: _bessel_K(5, 1.5, algorithm='scipy') in RR
        True

    """
    if algorithm=="scipy":
        if prec != 53:
            raise ValueError, "for the scipy algorithm the precision must be 53"
        import scipy.special
        ans = str(scipy.special.kv(float(nu),float(z)))
        ans = ans.replace("(","")
        ans = ans.replace(")","")
        ans = ans.replace("j","*I")
        ans = sage_eval(ans)
        return real(ans) if z in RR else ans
    elif algorithm == 'pari':
        from sage.libs.pari.all import pari
        try:
            R = RealField(prec)
            nu = R(nu)
            z = R(z)
        except TypeError:
            C = ComplexField(prec)
            nu = C(nu)
            z = C(z)
            K = C
        K = z.parent()
        return K(pari(nu).besselk(z, precision=prec))
    elif algorithm == 'maxima':
        raise NotImplementedError, "The K-Bessel function is only implemented for the pari and scipy algorithms"
    else:
        raise ValueError, "unknown algorithm '%s'"%algorithm


def _bessel_Y(nu,z,algorithm="maxima", prec=53):
    r"""
    Implements the "Y-Bessel function", or "Bessel function of the 2nd
    kind", with index (or "order") nu and argument z.

    .. note::

       Currently only prec=53 is supported.

    Defn::

                    cos(pi n)*bessel_J(nu, z) - bessel_J(-nu, z)
                   -------------------------------------------------
                                     sin(nu*pi)

    if nu is not an integer and by taking a limit otherwise.

    Sometimes bessel_Y(n,z) is denoted Y_n(z) in the literature.

    This is computed using Maxima by default.

    EXAMPLES::

        sage: from sage.functions.bessel import _bessel_Y
        sage: _bessel_Y(2,1.1,"scipy")
        -1.4314714939...
        sage: _bessel_Y(2,1.1)
        -1.4314714939590...
        sage: _bessel_Y(3.001,2.1)
        -1.0299574976424...

    TESTS::

        sage: _bessel_Y(2,1.1, algorithm="pari")
        Traceback (most recent call last):
        ...
        NotImplementedError: The Y-Bessel function is only implemented for the maxima and scipy algorithms
    """
    if algorithm=="scipy":
        if prec != 53:
            raise ValueError, "for the scipy algorithm the precision must be 53"
        import scipy.special
        ans = str(scipy.special.yv(float(nu),complex(real(z),imag(z))))
        ans = ans.replace("(","")
        ans = ans.replace(")","")
        ans = ans.replace("j","*I")
        ans = sage_eval(ans)
        return real(ans) if z in RR else ans
    elif algorithm == "maxima":
        if prec != 53:
            raise ValueError, "for the maxima algorithm the precision must be 53"
        return RR(maxima.eval("bessel_y(%s,%s)"%(float(nu),float(z))))
    elif algorithm == "pari":
        raise NotImplementedError, "The Y-Bessel function is only implemented for the maxima and scipy algorithms"
    else:
        raise ValueError, "unknown algorithm '%s'"%algorithm

class _Bessel():
    """
    A class implementing the I, J, K, and Y Bessel functions.

    EXAMPLES::

        sage: from sage.functions.bessel import _Bessel
        sage: g = _Bessel(2); g
        J_{2}
        sage: print g
        J-Bessel function of order 2
        sage: g.plot(0,10)

    ::

        sage: _Bessel(2, typ='I')(pi)
        2.61849485263445
        sage: _Bessel(2, typ='J')(pi)
        0.485433932631509
        sage: _Bessel(2, typ='K')(pi)
        0.0510986902537926
        sage: _Bessel(2, typ='Y')(pi)
        -0.0999007139289404
    """
    def __init__(self, nu, typ = "J", algorithm = None, prec = 53):
        """
        Initializes new instance of the Bessel class.

        INPUT:

         - ``typ`` -- (default: J) the type of Bessel function: 'I', 'J', 'K'
           or 'Y'.

         - ``algorithm`` -- (default: maxima for type Y, pari for other types)
           algorithm to use to compute the Bessel function: 'pari', 'maxima' or
           'scipy'.  Note that type K is not implemented in Maxima and type Y
           is not implemented in PARI.

         - ``prec`` -- (default: 53) precision in bits of the Bessel function.
           Only supported for the PARI algorithm.

        EXAMPLES::

            sage: from sage.functions.bessel import _Bessel
            sage: g = _Bessel(2); g
            J_{2}
            sage: _Bessel(1,'I')
            I_{1}
            sage: _Bessel(6, prec=120)(pi)
            0.014545966982505560573660369604001804
            sage: _Bessel(6, algorithm="pari")(pi)
            0.0145459669825056

        For the Bessel J-function, Maxima returns a symbolic result.  For
        types I and Y, we always get a numeric result::

            sage: b = _Bessel(6, algorithm="maxima")(pi); b
            bessel_j(6,pi)
            sage: b.n(53)
            0.0145459669825056
            sage: _Bessel(6, typ='I', algorithm="maxima")(pi)
            0.0294619840059568
            sage: _Bessel(6, typ='Y', algorithm="maxima")(pi)
            -4.33932818939038

        SciPy usually gives less precise results::

            sage: _Bessel(6, algorithm="scipy")(pi)
            0.0145459669825000...

        TESTS::

            sage: _Bessel(1,'Z')
            Traceback (most recent call last):
            ...
            ValueError: typ must be one of I, J, K, Y
        """
        if not (typ in ['I', 'J', 'K', 'Y']):
            raise ValueError, "typ must be one of I, J, K, Y"

        # Did the user ask for the default algorithm?
        if algorithm is None:
            if typ == 'Y':
                algorithm = 'maxima'
            else:
                algorithm = 'pari'

        self._system = algorithm
        self._order = nu
        self._type = typ
        prec = int(prec)
        if prec < 0:
            raise ValueError, "prec must be a positive integer"
        self._prec = int(prec)

    def __str__(self):
        """
        Returns a string representation of this Bessel object.

        TEST::

            sage: from sage.functions.bessel import _Bessel
            sage: a = _Bessel(1,'I')
            sage: str(a)
            'I-Bessel function of order 1'
        """
        return self.type()+"-Bessel function of order "+str(self.order())

    def __repr__(self):
        """
        Returns a string representation of this Bessel object.

        TESTS::

            sage: from sage.functions.bessel import _Bessel
            sage: _Bessel(1,'I')
            I_{1}
        """
        return self.type()+"_{"+str(self.order())+"}"

    def type(self):
        """
        Returns the type of this Bessel object.

        TEST::

            sage: from sage.functions.bessel import _Bessel
            sage: a = _Bessel(3,'K')
            sage: a.type()
            'K'
        """
        return self._type

    def prec(self):
        """
        Returns the precision (in number of bits) used to represent this
        Bessel function.

        TESTS::

            sage: from sage.functions.bessel import _Bessel
            sage: a = _Bessel(3,'K')
            sage: a.prec()
            53
            sage: B = _Bessel(20,prec=100); B
            J_{20}
            sage: B.prec()
            100
        """
        return self._prec

    def order(self):
        """
        Returns the order of this Bessel function.

        TEST::

            sage: from sage.functions.bessel import _Bessel
            sage: a = _Bessel(3,'K')
            sage: a.order()
            3
        """
        return self._order

    def system(self):
        """
        Returns the package used, e.g. Maxima, PARI, or SciPy, to compute with
        this Bessel function.

        TESTS::

            sage: from sage.functions.bessel import _Bessel
            sage: _Bessel(20,algorithm='maxima').system()
            'maxima'
            sage: _Bessel(20,prec=100).system()
            'pari'
        """
        return self._system

    def __call__(self,z):
        """
        Implements evaluation of all the Bessel functions directly
        from the Bessel class. This essentially allows a Bessel object to
        behave like a function that can be invoked.

        TESTS::

            sage: from sage.functions.bessel import _Bessel
            sage: _Bessel(3,'K')(5.0)
            0.00829176841523093
            sage: _Bessel(20,algorithm='maxima')(5.0)
            27.703300521289436e-12
            sage: _Bessel(20,prec=100)(5.0101010101010101)
            2.8809188227195382093062257967e-11
            sage: B = _Bessel(2,'Y',algorithm='scipy',prec=50)
            sage: B(2.0)
            Traceback (most recent call last):
            ...
            ValueError: for the scipy algorithm the precision must be 53
        """
        nu = self.order()
        t = self.type()
        s = self.system()
        p = self.prec()
        if t == "I":
            return _bessel_I(nu,z,algorithm=s,prec=p)
        if t == "J":
            return _bessel_J(nu,z,algorithm=s,prec=p)
        if t == "K":
            return _bessel_K(nu,z,algorithm=s,prec=p)
        if t == "Y":
            return _bessel_Y(nu,z,algorithm=s,prec=p)

    def plot(self,a,b):
        """
        Enables easy plotting of all the Bessel functions directly
        from the Bessel class.

        TESTS::

            sage: from sage.functions.bessel import _Bessel
            sage: plot(_Bessel(2),3,4)
            sage: _Bessel(2).plot(3,4)
            sage: P = _Bessel(2,'I').plot(1,5)
            sage: P += _Bessel(2,'J').plot(1,5)
            sage: P += _Bessel(2,'K').plot(1,5)
            sage: P += _Bessel(2,'Y').plot(1,5)
            sage: show(P)
        """
        nu = self.order()
        s = self.system()
        t = self.type()
        if t == "I":
            f = lambda z: _bessel_I(nu,z,s)
            P = plot(f,a,b)
        if t == "J":
            f = lambda z: _bessel_J(nu,z,s)
            P = plot(f,a,b)
        if t == "K":
            f = lambda z: _bessel_K(nu,z,s)
            P = plot(f,a,b)
        if t == "Y":
            f = lambda z: _bessel_Y(nu,z,s)
            P = plot(f,a,b)
        return P
