r"""
Exponential Integrals

AUTHORS:

- Benjamin Jones (2011-06-12)

This module provides easy access to many exponential integral
special functions. It utilizes Maxima's `special functions package`_ and
the `mpmath library`_.

REFERENCES:

- [AS]_ Abramowitz and Stegun: *Handbook of Mathematical Functions*
- Wikipedia Entry: http://en.wikipedia.org/wiki/Exponential_integral
- Online Encyclopedia of Special Function: http://algo.inria.fr/esf/index.html
- NIST Digital Library of Mathematical Functions: http://dlmf.nist.gov/
- Maxima `special functions package`_
- `mpmath library`_

.. [AS] 'Handbook of Mathematical Functions', Milton Abramowitz and Irene
   A. Stegun, National Bureau of Standards Applied Mathematics Series, 55.
   See also http://www.math.sfu.ca/~cbm/aands/.
.. _`special functions package`: http://maxima.sourceforge.net/docs/manual/en/maxima_15.html
.. _`mpmath library`: http://code.google.com/p/mpmath/

AUTHORS:

- Benjamin Jones

    Implementations of the classes ``Function_exp_integral_*``.

- David Joyner and William Stein

    Authors of the code which was moved from special.py and trans.py.
    Implementation of :meth:`exp_int` (from sage/functions/special.py).
    Implementation of :meth:`exponential_integral_1` (from
    sage/functions/transcendental.py).

"""

#*****************************************************************************
#       Copyright (C) 2011 Benjamin Jones <benjaminfjones@gmail.com>
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

import sage.interfaces.all
from sage.misc.sage_eval import sage_eval
from sage.symbolic.function import BuiltinFunction, is_inexact
from sage.calculus.calculus import maxima
from sage.symbolic.expression import Expression
from sage.structure.parent import Parent
from sage.structure.coerce import parent
from sage.libs.mpmath import utils as mpmath_utils
mpmath_utils_call = mpmath_utils.call # eliminate some overhead in _evalf_

from sage.rings.rational_field import RationalField
from sage.rings.real_mpfr import RealField
from sage.rings.complex_field import ComplexField
from sage.rings.all import ZZ, QQ, RR, RDF
from sage.functions.log import exp, log
from sage.functions.trig import sin, cos
from sage.functions.hyperbolic import sinh, cosh


class Function_exp_integral_e(BuiltinFunction):
    r"""
    The generalized complex exponential integral `E_n(z)` defined by

    .. math::

        \operatorname{E_n}(z) = \int_1^{\infty} \frac{e^{-z t}}{t^n} \; dt

    for complex numbers `n` and `z`, see [AS]_ 5.1.4.

    The special case where `n = 1` is denoted in Sage by
    ``exp_integral_e1``.

    EXAMPLES:

    Numerical evaluation is handled using mpmath::

        sage: N(exp_integral_e(1,1))
        0.219383934395520
        sage: exp_integral_e(1, RealField(100)(1))
        0.21938393439552027367716377546

    We can compare this to PARI's evaluation of
    :meth:`exponential_integral_1`::

        sage: N(exponential_integral_1(1))
        0.219383934395520

    We can verify one case of [AS]_ 5.1.45, i.e.
    `E_n(z) = z^{n-1}\Gamma(1-n,z)`::

        sage: N(exp_integral_e(2, 3+I))
        0.00354575823814662 - 0.00973200528288687*I
        sage: N((3+I)*gamma(-1, 3+I))
        0.00354575823814662 - 0.00973200528288687*I

    Maxima returns the following improper integral as a multiple of
    ``exp_integral_e(1,1)``::

        sage: uu = integral(e^(-x)*log(x+1),x,0,oo)
        sage: uu
        e*exp_integral_e(1, 1)
        sage: uu.n(digits=30)
        0.596347362323194074341078499369

    Symbolic derivatives and integrals are handled by Sage and Maxima::

        sage: x = var('x')
        sage: f = exp_integral_e(2,x)
        sage: f.diff(x)
        -exp_integral_e(1, x)

        sage: f.integrate(x)
        -exp_integral_e(3, x)

        sage: f = exp_integral_e(-1,x)
        sage: f.integrate(x)
        Ei(-x) - gamma(-1, x)

    Some special values of ``exp_integral_e`` can be simplified.
    [AS]_ 5.1.23::

        sage: exp_integral_e(0,x)
        e^(-x)/x

    [AS]_ 5.1.24::

        sage: exp_integral_e(6,0)
        1/5
        sage: nn = var('nn')
        sage: assume(nn > 1)
        sage: f = exp_integral_e(nn,0)
        sage: f.simplify()
        1/(nn - 1)


    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    """
    def __init__(self):
        """
        See the docstring for :meth:`Function_exp_integral_e`.

        EXAMPLES::

            sage: exp_integral_e(1,0)
            exp_integral_e(1, 0)

        """
        BuiltinFunction.__init__(self, "exp_integral_e", nargs=2,
                                 latex_name=r'exp_integral_e',
                                 conversions=dict(maxima='expintegral_e'))

    def _eval_(self, n, z):
        """
        EXAMPLES::

            sage: exp_integral_e(1.0, x)
            exp_integral_e(1.00000000000000, x)
            sage: exp_integral_e(x, 1.0)
            exp_integral_e(x, 1.00000000000000)
            sage: exp_integral_e(1.0, 1.0)
            0.219383934395520

        """
        if not isinstance(n, Expression) and not isinstance(z, Expression) and \
               (is_inexact(n) or is_inexact(z)):
            coercion_model = sage.structure.element.get_coercion_model()
            n, z = coercion_model.canonical_coercion(n, z)
            return self._evalf_(n, z, parent(n))

        z_zero = False
        # special case: z == 0 and n > 1
        if isinstance(z, Expression):
            if z.is_trivial_zero():
                z_zero = True # for later
                if n > 1:
                    return 1/(n-1)
        else:
            if not z:
                z_zero = True
                if n > 1:
                    return 1/(n-1)

        # special case: n == 0
        if isinstance(n, Expression):
            if n.is_trivial_zero():
                if z_zero:
                    return None
                else:
                    return exp(-z)/z
        else:
            if not n:
                if z_zero:
                    return None
                else:
                    return exp(-z)/z

        return None # leaves the expression unevaluated

    def _evalf_(self, n, z, parent=None):
        """
        EXAMPLES::

            sage: N(exp_integral_e(1, 1+I))
            0.000281624451981418 - 0.179324535039359*I
            sage: exp_integral_e(1, RealField(100)(1))
            0.21938393439552027367716377546

        """
        import mpmath
        return mpmath_utils.call(mpmath.expint, n, z, parent=parent)

    def _derivative_(self, n, z, diff_param=None):
        """
        If `n` is an integer strictly larger than 0, then the derivative of
        `E_n(z)` with respect to `z` is
        `-E_{n-1}(z)`. See [AS]_ 5.1.26.

        EXAMPLES::

            sage: x = var('x')
            sage: f = exp_integral_e(2,x)
            sage: f.diff(x)
            -exp_integral_e(1, x)

            sage: f = exp_integral_e(2,sqrt(x))
            sage: f.diff(x)
            -1/2*exp_integral_e(1, sqrt(x))/sqrt(x)

        """
        if n in ZZ and n > 0:
            return -1*exp_integral_e(n-1,z)
        else:
            raise NotImplementedError("The derivative of this function is only implemented for n = 1, 2, 3, ...")

exp_integral_e = Function_exp_integral_e()


class Function_exp_integral_e1(BuiltinFunction):
    r"""
    The generalized complex exponential integral `E_1(z)` defined by

    .. math::

        \operatorname{E_1}(z) = \int_z^\infty \frac{e^{-t}}{t}\; dt

    see [AS]_ 5.1.4.

    EXAMPLES:

    Numerical evaluation is handled using mpmath::

        sage: N(exp_integral_e1(1))
        0.219383934395520
        sage: exp_integral_e1(RealField(100)(1))
        0.21938393439552027367716377546

    We can compare this to PARI's evaluation of
    :meth:`exponential_integral_1`::

        sage: N(exp_integral_e1(2.0))
        0.0489005107080611
        sage: N(exponential_integral_1(2.0))
        0.0489005107080611

    Symbolic derivatives and integrals are handled by Sage and Maxima::

        sage: x = var('x')
        sage: f = exp_integral_e1(x)
        sage: f.diff(x)
        -e^(-x)/x

        sage: f.integrate(x)
        -exp_integral_e(2, x)

    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    """
    def __init__(self):
        """
        See the docstring for :class:`Function_exp_integral_e1`.

        EXAMPLES::

            sage: exp_integral_e1(1)
            exp_integral_e1(1)

        """
        BuiltinFunction.__init__(self, "exp_integral_e1", nargs=1,
                                 latex_name=r'exp_integral_e1',
                                 conversions=dict(maxima='expintegral_e1'))

    def _eval_(self, z):
        """
        EXAMPLES::

            sage: exp_integral_e1(x)
            exp_integral_e1(x)
            sage: exp_integral_e1(1.0)
            0.219383934395520

        """
        if not isinstance(z, Expression) and is_inexact(z):
            return self._evalf_(z, parent(z))

        return None # leaves the expression unevaluated

    def _evalf_(self, z, parent=None):
        """
        EXAMPLES::

            sage: N(exp_integral_e1(1+I))
            0.000281624451981418 - 0.179324535039359*I
            sage: exp_integral_e1(RealField(200)(0.5))
            0.55977359477616081174679593931508523522684689031635351524829

        """
        import mpmath
        return mpmath_utils_call(mpmath.e1, z, parent=parent)

    def _derivative_(self, z, diff_param=None):
        """
        The derivative of `E_1(z)` is `-e^{-z}/z`. See [AS], 5.1.26.

        EXAMPLES::

            sage: x = var('x')
            sage: f = exp_integral_e1(x)
            sage: f.diff(x)
            -e^(-x)/x

            sage: f = exp_integral_e1(x^2)
            sage: f.diff(x)
            -2*e^(-x^2)/x

        """
        return -exp(-z)/z

exp_integral_e1 = Function_exp_integral_e1()


class Function_log_integral(BuiltinFunction):
    r"""
    The logarithmic integral `\operatorname{li}(z)` defined by

    .. math::

        \operatorname{li}(x) = \int_0^z \frac{dt}{\ln(t)} = \operatorname{Ei}(\ln(x))

    for x > 1 and by analytic continuation for complex arguments z (see [AS]_ 5.1.3).

    EXAMPLES:

    Numerical evaluation for real and complex arguments is handled using mpmath::

        sage: N(log_integral(3))
        2.16358859466719
        sage: N(log_integral(3), digits=30)
        2.16358859466719197287692236735
        sage: log_integral(ComplexField(100)(3+I))
        2.2879892769816826157078450911 + 0.87232935488528370139883806779*I
        sage: log_integral(0)
        0

    Symbolic derivatives and integrals are handled by Sage and Maxima::

        sage: x = var('x')
        sage: f = log_integral(x)
        sage: f.diff(x)
        1/log(x)

        sage: f.integrate(x)
        x*log_integral(x) - Ei(2*log(x))

    Here is a test from the mpmath documentation. There are
    1,925,320,391,606,803,968,923 many prime numbers less than 1e23. The
    value of ``log_integral(1e23)`` is very close to this::

        sage: log_integral(1e23)
        1.92532039161405e21

    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    REFERENCES:

    - http://en.wikipedia.org/wiki/Logarithmic_integral_function
    - mpmath documentation: `logarithmic-integral`_

    .. _`logarithmic-integral`: http://mpmath.googlecode.com/svn/trunk/doc/build/functions/expintegrals.html#logarithmic-integral


    """
    def __init__(self):
        """
        See the docstring for ``Function_log_integral``.

        EXAMPLES::

            sage: log_integral(3)
            log_integral(3)

        """
        BuiltinFunction.__init__(self, "log_integral", nargs=1,
                                 latex_name=r'log_integral',
                                 conversions=dict(maxima='expintegral_li'))

    def _eval_(self, z):
        """
        EXAMPLES::

            sage: z = var('z')
            sage: log_integral(z)
            log_integral(z)
            sage: log_integral(3.0)
            2.16358859466719
            sage: log_integral(0)
            0

        """
        if isinstance(z, Expression):
            if z.is_trivial_zero():         # special case: z = 0
                return z
        else:
            if is_inexact(z):
                return self._evalf_(z, parent(z))
            elif not z:
                return z
        return None # leaves the expression unevaluated

    def _evalf_(self, z, parent=None):
        """
        EXAMPLES::

            sage: N(log_integral(1e6))
            78627.5491594622
            sage: log_integral(RealField(200)(1e6))
            78627.549159462181919862910747947261161321874382421767074759

        """
        import mpmath
        return mpmath_utils_call(mpmath.li, z, parent=parent)

    def _derivative_(self, z, diff_param=None):
        r"""
        The derivative of `\operatorname{li}(z) is `1/log(z)`.

        EXAMPLES::

            sage: x = var('x')
            sage: f = log_integral(x)
            sage: f.diff(x)
            1/log(x)

            sage: f = log_integral(x^2)
            sage: f.diff(x)
            2*x/log(x^2)

        """
        return 1/log(z)

li = log_integral = Function_log_integral()


class Function_sin_integral(BuiltinFunction):
    r"""
    The trigonometric integral `\operatorname{Si}(z)` defined by

    .. math::

        \operatorname{Si}(z) = \int_0^z \frac{\sin(t)}{t}\; dt,

    see [AS]_ 5.2.1.

    EXAMPLES:

    Numerical evaluation for real and complex arguments is handled using mpmath::

        sage: sin_integral(0)
        0
        sage: sin_integral(0.0)
        0.000000000000000
        sage: sin_integral(3.0)
        1.84865252799947
        sage: N(sin_integral(3), digits=30)
        1.84865252799946825639773025111
        sage: sin_integral(ComplexField(100)(3+I))
        2.0277151656451253616038525998 + 0.015210926166954211913653130271*I

    The alias `Si` can be used instead of `sin_integral`::

        sage: Si(3.0)
        1.84865252799947

    The limit of `\operatorname{Si}(z)` as `z \to \infty` is `\pi/2`::

        sage: N(sin_integral(1e23))
        1.57079632679490
        sage: N(pi/2)
        1.57079632679490

    At 200 bits of precision `\operatorname{Si}(10^{23})` agrees with `\pi/2` up to
    `10^{-24}`::

        sage: sin_integral(RealField(200)(1e23))
        1.5707963267948966192313288218697837425815368604836679189519
        sage: N(pi/2, prec=200)
        1.5707963267948966192313216916397514420985846996875529104875

    The exponential sine integral is analytic everywhere::

        sage: sin_integral(-1.0)
        -0.946083070367183
        sage: sin_integral(-2.0)
        -1.60541297680269
        sage: sin_integral(-1e23)
        -1.57079632679490

    Symbolic derivatives and integrals are handled by Sage and Maxima::

        sage: x = var('x')
        sage: f = sin_integral(x)
        sage: f.diff(x)
        sin(x)/x

        sage: f.integrate(x)
        x*sin_integral(x) + cos(x)

        sage: integrate(sin(x)/x, x)
        1/2*I*Ei(-I*x) - 1/2*I*Ei(I*x)

    Compare values of the functions `\operatorname{Si}(x)` and
    `f(x) = (1/2)i \cdot \operatorname{Ei}(-ix) - (1/2)i \cdot
    \operatorname{Ei}(ix) - \pi/2`, which are both anti-derivatives of
    `\sin(x)/x`, at some random positive real numbers::

        sage: f(x) = 1/2*I*Ei(-I*x) - 1/2*I*Ei(I*x) - pi/2
        sage: g(x) = sin_integral(x)
        sage: R = [ abs(RDF.random_element()) for i in range(100) ]
        sage: all(abs(f(x) - g(x)) < 1e-10 for x in R)
        True

    The Nielsen spiral is the parametric plot of (Si(t), Ci(t))::

        sage: x=var('x')
        sage: f(x) = sin_integral(x)
        sage: g(x) = cos_integral(x)
        sage: P = parametric_plot([f, g], (x, 0.5 ,20))
        sage: show(P, frame=True, axes=False)

    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    REFERENCES:

    - http://en.wikipedia.org/wiki/Trigonometric_integral
    - mpmath documentation: `si`_

    .. _`si`: http://mpmath.googlecode.com/svn/trunk/doc/build/functions/expintegrals.html#si

    """
    def __init__(self):
        """
        See the docstring for ``Function_sin_integral``.

        EXAMPLES::

            sage: sin_integral(1)
            sin_integral(1)

        """
        BuiltinFunction.__init__(self, "sin_integral", nargs=1,
                                 latex_name=r'\operatorname{Si}',
                                 conversions=dict(maxima='expintegral_si'))

    def _eval_(self, z):
        """
        EXAMPLES::

            sage: z = var('z')
            sage: sin_integral(z)
            sin_integral(z)
            sage: sin_integral(3.0)
            1.84865252799947
            sage: sin_integral(0)
            0

        """
        if not isinstance(z, Expression) and is_inexact(z):
            return self._evalf_(z, parent(z))

        # special case: z = 0
        if isinstance(z, Expression):
            if z.is_trivial_zero():
                return z
        else:
            if not z:
                return z

        return None # leaves the expression unevaluated

    def _evalf_(self, z, parent=None):
        """
        EXAMPLES:

        The limit `\operatorname{Si}(z)` as `z \to \infty`  is `\pi/2`::

            sage: N(sin_integral(1e23) - pi/2)
            0.000000000000000

        At 200 bits of precision `\operatorname{Si}(10^{23})` agrees with `\pi/2` up to
        `10^{-24}`::

            sage: sin_integral(RealField(200)(1e23))
            1.5707963267948966192313288218697837425815368604836679189519
            sage: N(pi/2, prec=200)
            1.5707963267948966192313216916397514420985846996875529104875

        The exponential sine integral is analytic everywhere, even on the
        negative real axis::

            sage: sin_integral(-1.0)
            -0.946083070367183
            sage: sin_integral(-2.0)
            -1.60541297680269
            sage: sin_integral(-1e23)
            -1.57079632679490

        """
        import mpmath
        return mpmath_utils_call(mpmath.si, z, parent=parent)

    def _derivative_(self, z, diff_param=None):
        r"""
        The derivative of `\operatorname{Si}(z)` is `\sin(z)/z` if `z`
        is not zero. The derivative at `z = 0` is `1` (but this
        exception is not currently implemented).

        EXAMPLES::

            sage: x = var('x')
            sage: f = sin_integral(x)
            sage: f.diff(x)
            sin(x)/x

            sage: f = sin_integral(x^2)
            sage: f.diff(x)
            2*sin(x^2)/x

        """
        return sin(z)/z

Si = sin_integral = Function_sin_integral()


class Function_cos_integral(BuiltinFunction):
    r"""
    The trigonometric integral `\operatorname{Ci}(z)` defined by

    .. math::

        \operatorname{Ci}(z) = \gamma + \log(z) + \int_0^z \frac{\cos(t)-1}{t}\; dt,

    where `\gamma` is the Euler gamma constant (``euler_gamma`` in Sage),
    see [AS]_ 5.2.1.

    EXAMPLES:

    Numerical evaluation for real and complex arguments is handled using mpmath::

        sage: cos_integral(3.0)
        0.119629786008000

    The alias `Ci` can be used instead of `cos_integral`::

        sage: Ci(3.0)
        0.119629786008000

    Compare ``cos_integral(3.0)`` to the definition of the value using
    numerical integration::

        sage: N(euler_gamma + log(3.0) + integrate((cos(x)-1)/x, x, 0, 3.0) - cos_integral(3.0)) < 1e-14
        True

    Arbitrary precision and complex arguments are handled::

        sage: N(cos_integral(3), digits=30)
        0.119629786008000327626472281177
        sage: cos_integral(ComplexField(100)(3+I))
        0.078134230477495714401983633057 - 0.37814733904787920181190368789*I

    The limit `\operatorname{Ci}(z)` as `z \to \infty` is zero::

        sage: N(cos_integral(1e23))
        -3.24053937643003e-24

    Symbolic derivatives and integrals are handled by Sage and Maxima::

        sage: x = var('x')
        sage: f = cos_integral(x)
        sage: f.diff(x)
        cos(x)/x

        sage: f.integrate(x)
        x*cos_integral(x) - sin(x)

    The Nielsen spiral is the parametric plot of (Si(t), Ci(t))::

        sage: t=var('t')
        sage: f(t) = sin_integral(t)
        sage: g(t) = cos_integral(t)
        sage: P = parametric_plot([f, g], (t, 0.5 ,20))
        sage: show(P, frame=True, axes=False)

    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    REFERENCES:

    - http://en.wikipedia.org/wiki/Trigonometric_integral
    - mpmath documentation: `ci`_

    .. _`ci`: http://mpmath.googlecode.com/svn/trunk/doc/build/functions/expintegrals.html#ci

    """
    def __init__(self):
        """
        See the docstring for :class:`Function_cos_integral`.

        EXAMPLES::

            sage: cos_integral(1)
            cos_integral(1)

        """
        BuiltinFunction.__init__(self, "cos_integral", nargs=1,
                                 latex_name=r'\operatorname{Ci}',
                                 conversions=dict(maxima='expintegral_ci'))

    def _eval_(self, z):
        """
        EXAMPLES::

            sage: z = var('z')
            sage: cos_integral(z)
            cos_integral(z)
            sage: cos_integral(3.0)
            0.119629786008000
            sage: cos_integral(0)
            cos_integral(0)
            sage: N(cos_integral(0))
            -infinity

        """
        if not isinstance(z, Expression) and is_inexact(z):
            return self._evalf_(z, parent(z))

        return None # leaves the expression unevaluated

    def _evalf_(self, z, parent=None):
        """
        EXAMPLES::

            sage: N(cos_integral(1e23)) < 1e-20
            True
            sage: N(cos_integral(1e-10), digits=30)
            -22.4486352650389235918737540487
            sage: cos_integral(ComplexField(100)(I))
            0.83786694098020824089467857943 + 1.5707963267948966192313216916*I

        """
        import mpmath
        return mpmath_utils_call(mpmath.ci, z, parent=parent)

    def _derivative_(self, z, diff_param=None):
        r"""
        The derivative of `\operatorname{Ci}(z)` is `\cos(z)/z` if `z` is not zero.

        EXAMPLES::

            sage: x = var('x')
            sage: f = cos_integral(x)
            sage: f.diff(x)
            cos(x)/x

            sage: f = cos_integral(x^2)
            sage: f.diff(x)
            2*cos(x^2)/x

        """
        return cos(z)/z

Ci = cos_integral = Function_cos_integral()


class Function_sinh_integral(BuiltinFunction):
    r"""
    The trigonometric integral `\operatorname{Shi}(z)` defined by

    .. math::

        \operatorname{Shi}(z) = \int_0^z \frac{\sinh(t)}{t}\; dt,

    see [AS]_ 5.2.3.

    EXAMPLES:

    Numerical evaluation for real and complex arguments is handled using mpmath::

        sage: sinh_integral(3.0)
        4.97344047585981
        sage: sinh_integral(1.0)
        1.05725087537573
        sage: sinh_integral(-1.0)
        -1.05725087537573

    The alias `Shi` can be used instead of `sinh_integral`::

        sage: Shi(3.0)
        4.97344047585981

    Compare ``sinh_integral(3.0)`` to the definition of the value using
    numerical integration::

        sage: N(integrate((sinh(x))/x, x, 0, 3.0) - sinh_integral(3.0)) < 1e-14
        True

    Arbitrary precision and complex arguments are handled::

        sage: N(sinh_integral(3), digits=30)
        4.97344047585980679771041838252
        sage: sinh_integral(ComplexField(100)(3+I))
        3.9134623660329374406788354078 + 3.0427678212908839256360163759*I

    The limit `\operatorname{Shi}(z)` as `z \to \infty` is `\infty`::

        sage: N(sinh_integral(Infinity))
        +infinity

    Symbolic derivatives and integrals are handled by Sage and Maxima::

        sage: x = var('x')
        sage: f = sinh_integral(x)
        sage: f.diff(x)
        sinh(x)/x

        sage: f.integrate(x)
        x*sinh_integral(x) - cosh(x)

    Note that due to some problems with the way Maxima handles these
    expressions, definite integrals can sometimes give unexpected
    results (typically when using inexact endpoints) due to
    inconsistent branching::

        sage: integrate(sinh_integral(x), x, 0, 1/2)
        -cosh(1/2) + 1/2*sinh_integral(1/2) + 1
        sage: integrate(sinh_integral(x), x, 0, 1/2).n() # correct
        0.125872409703453
        sage: integrate(sinh_integral(x), x, 0, 0.5).n() # incorrect!
        0.125872409703453 + 1.57079632679490*I

    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    REFERENCES:

    - http://en.wikipedia.org/wiki/Trigonometric_integral
    - mpmath documentation: `shi`_

    .. _`shi`: http://mpmath.googlecode.com/svn/trunk/doc/build/functions/expintegrals.html#shi

    """
    def __init__(self):
        """
        See the docstring for ``Function_sinh_integral``.

        EXAMPLES::

            sage: sinh_integral(1)
            sinh_integral(1)

        """
        BuiltinFunction.__init__(self, "sinh_integral", nargs=1,
                                 latex_name=r'\operatorname{Shi}',
                                 conversions=dict(maxima='expintegral_shi'))

    def _eval_(self, z):
        """
        EXAMPLES::

            sage: z = var('z')
            sage: sinh_integral(z)
            sinh_integral(z)
            sage: sinh_integral(3.0)
            4.97344047585981
            sage: sinh_integral(0)
            0

        """
        if not isinstance(z, Expression) and is_inexact(z):
            return self._evalf_(z, parent(z))

        # special case: z = 0
        if isinstance(z, Expression):
            if z.is_trivial_zero():
                return z
        else:
            if not z:
                return z

        return None # leaves the expression unevaluated

    def _evalf_(self, z, parent=None):
        """
        EXAMPLES::

            sage: N(sinh_integral(1e-10), digits=30)
            1.00000000000000003643219731550e-10
            sage: sinh_integral(ComplexField(100)(I))
            0.94608307036718301494135331382*I

        """
        import mpmath
        return mpmath_utils_call(mpmath.shi, z, parent=parent)

    def _derivative_(self, z, diff_param=None):
        r"""
        The derivative of `\operatorname{Shi}(z)` is `\sinh(z)/z`.

        EXAMPLES::

            sage: x = var('x')
            sage: f = sinh_integral(x)
            sage: f.diff(x)
            sinh(x)/x

            sage: f = sinh_integral(ln(x))
            sage: f.diff(x)
            sinh(log(x))/(x*log(x))

        """
        return sinh(z)/z

Shi = sinh_integral = Function_sinh_integral()


class Function_cosh_integral(BuiltinFunction):
    r"""
    The trigonometric integral `\operatorname{Chi}(z)` defined by

    .. math::

        \operatorname{Chi}(z) = \gamma + \log(z) + \int_0^z \frac{\cosh(t)-1}{t}\; dt,

    see [AS]_ 5.2.4.

    EXAMPLES:

    Numerical evaluation for real and complex arguments is handled using mpmath::

        sage: cosh_integral(1.0)
        0.837866940980208

    The alias `Chi` can be used instead of `cosh_integral`::

        sage: Chi(1.0)
        0.837866940980208

    Here is an example from the mpmath documentation::

        sage: f(x) = cosh_integral(x)
        sage: find_root(f, 0.1, 1.0)
        0.5238225713894826

    Compare ``cosh_integral(3.0)`` to the definition of the value using
    numerical integration::

        sage: N(euler_gamma + log(3.0) + integrate((cosh(x)-1)/x, x, 0, 3.0) -
        ...     cosh_integral(3.0)) < 1e-14
        True

    Arbitrary precision and complex arguments are handled::

        sage: N(cosh_integral(3), digits=30)
        4.96039209476560976029791763669
        sage: cosh_integral(ComplexField(100)(3+I))
        3.9096723099686417127843516794 + 3.0547519627014217273323873274*I

    The limit of `\operatorname{Chi}(z)` as `z \to \infty` is `\infty`::

        sage: N(cosh_integral(Infinity))
        +infinity

    Symbolic derivatives and integrals are handled by Sage and Maxima::

        sage: x = var('x')
        sage: f = cosh_integral(x)
        sage: f.diff(x)
        cosh(x)/x

        sage: f.integrate(x)
        x*cosh_integral(x) - sinh(x)

    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    REFERENCES:

    - http://en.wikipedia.org/wiki/Trigonometric_integral
    - mpmath documentation: `chi`_

    .. _`chi`: http://mpmath.googlecode.com/svn/trunk/doc/build/functions/expintegrals.html#chi

    """
    def __init__(self):
        """
        See the docstring for ``Function_cosh_integral``.

        EXAMPLES::

            sage: cosh_integral(1)
            cosh_integral(1)

        """
        BuiltinFunction.__init__(self, "cosh_integral", nargs=1,
                                 latex_name=r'\operatorname{Chi}',
                                 conversions=dict(maxima='expintegral_chi'))

    def _eval_(self, z):
        """
        EXAMPLES::

            sage: z = var('z')
            sage: cosh_integral(z)
            cosh_integral(z)
            sage: cosh_integral(3.0)
            4.96039209476561

        """
        if not isinstance(z, Expression) and is_inexact(z):
            return self._evalf_(z, parent(z))

        return None

    def _evalf_(self, z, parent=None):
        """
        EXAMPLES::

            sage: N(cosh_integral(1e-10), digits=30)
            -22.4486352650389235918737540487
            sage: cosh_integral(ComplexField(100)(I))
            0.33740392290096813466264620389 + 1.5707963267948966192313216916*I

        """
        import mpmath
        return mpmath_utils_call(mpmath.chi, z, parent=parent)

    def _derivative_(self, z, diff_param=None):
        """
        The derivative of `\operatorname{Chi}(z)` is `\cosh(z)/z`.

        EXAMPLES::

            sage: x = var('x')
            sage: f = cosh_integral(x)
            sage: f.diff(x)
            cosh(x)/x

            sage: f = cosh_integral(ln(x))
            sage: f.diff(x)
            cosh(log(x))/(x*log(x))

        """
        return cosh(z)/z

Chi = cosh_integral = Function_cosh_integral()


###################################################################
## Code below here was moved from sage/functions/transcendental.py
## This occured as part of Trac #11143.
###################################################################
#
# This class has a name which is not specific enough
# see Function_exp_integral_e above, for example, which
# is the "generalized" exponential integral function. We
# are leaving the name the same for backwards compatibility
# purposes.
class Function_exp_integral(BuiltinFunction):
    r"""
    The generalized complex exponential integral Ei(z) defined by

    .. math::

        \operatorname{Ei}(x) = \int_{-\infty}^x \frac{e^t}{t}\; dt

    for x > 0 and for complex arguments by analytic continuation,
    see [AS]_ 5.1.2.

    EXAMPLES::

        sage: Ei(10)
        Ei(10)
        sage: Ei(I)
        Ei(I)
        sage: Ei(3+I)
        Ei(I + 3)
        sage: Ei(1.3)
        2.72139888023202

    The branch cut for this function is along the negative real axis::

        sage: Ei(-3 + 0.1*I)
        -0.0129379427181693 + 3.13993830250942*I
        sage: Ei(-3 - 0.1*I)
        -0.0129379427181693 - 3.13993830250942*I

    ALGORITHM: Uses mpmath.

    TESTS:

    Show that the evaluation and limit issue in :trac:`13271` is fixed::

        sage: var('Z')
        Z
        sage: (Ei(-Z)).limit(Z=oo)
        0
        sage: (Ei(-Z)).limit(Z=1000)
        Ei(-1000)
        sage: (Ei(-Z)).limit(Z=1000).n()
        -5.07089306023517e-438
    """
    def __init__(self):
        """
        Return the value of the complex exponential integral Ei(z) at a
        complex number z.

        EXAMPLES::

            sage: Ei(10)
            Ei(10)
            sage: Ei(I)
            Ei(I)
            sage: Ei(3+I)
            Ei(I + 3)
            sage: Ei(1.3)
            2.72139888023202

        The branch cut for this function is along the negative real axis::

            sage: Ei(-3 + 0.1*I)
            -0.0129379427181693 + 3.13993830250942*I
            sage: Ei(-3 - 0.1*I)
            -0.0129379427181693 - 3.13993830250942*I

        ALGORITHM: Uses mpmath.
        """
        BuiltinFunction.__init__(self, "Ei",
                                 conversions=dict(maxima='expintegral_ei'))

    def _eval_(self, x ):
        """
        EXAMPLES::

            sage: Ei(10)
            Ei(10)
            sage: Ei(I)
            Ei(I)
            sage: Ei(1.3)
            2.72139888023202
            sage: Ei(10r)
            Ei(10)
            sage: Ei(1.3r)
            2.7213988802320235
        """
        if not isinstance(x, Expression) and is_inexact(x):
            return self._evalf_(x, parent(x))
        return None

    def _evalf_(self, x, parent=None):
        """
        EXAMPLES::

            sage: Ei(10).n()
            2492.22897624188
            sage: Ei(20).n()
            2.56156526640566e7
            sage: Ei(I).n()
            0.337403922900968 + 2.51687939716208*I
            sage: Ei(3+I).n()
            7.82313467600158 + 6.09751978399231*I
        """
        import mpmath
        return mpmath_utils_call(mpmath.ei, x, parent=parent)

    def __call__(self, x, prec=None, coerce=True, hold=False ):
        """
        Note that the ``prec`` argument is deprecated. The precision for
        the result is deduced from the precision of the input. Convert
        the input to a higher precision explicitly if a result with higher
        precision is desired.

        EXAMPLES::

            sage: t = Ei(RealField(100)(2.5)); t
            7.0737658945786007119235519625
            sage: t.prec()
            100

            sage: Ei(1.1, prec=300)
            doctest:...: DeprecationWarning: The prec keyword argument is deprecated. Explicitly set the precision of the input, for example Ei(RealField(300)(1)), or use the prec argument to .n() for exact inputs, e.g., Ei(1).n(300), instead.
            See http://trac.sagemath.org/7748 for details.
            2.16737827956340306615064476647912607220394065907142504328679588538509331805598360907980986
        """
        if prec is not None:
            from sage.misc.superseded import deprecation
            deprecation(7748, "The prec keyword argument is deprecated. Explicitly set the precision of the input, for example Ei(RealField(300)(1)), or use the prec argument to .n() for exact inputs, e.g., Ei(1).n(300), instead.")
            import mpmath
            return mpmath_utils_call(mpmath.ei, x, prec=prec)

        return BuiltinFunction.__call__(self, x, coerce=coerce, hold=hold)

    def _derivative_(self, x, diff_param=None):
        """
        EXAMPLES::

            sage: Ei(x).diff(x)
            e^x/x
            sage: Ei(x).diff(x).subs(x=1)
            e
            sage: Ei(x^2).diff(x)
            2*e^(x^2)/x
            sage: f = function('f')
            sage: Ei(f(x)).diff(x)
            e^f(x)*D[0](f)(x)/f(x)
        """
        return exp(x)/x

Ei = exp_integral_ei = Function_exp_integral()


# moved here from sage/functions/transcendental.py
def exponential_integral_1(x, n=0):
    r"""
    Returns the exponential integral `E_1(x)`. If the optional
    argument `n` is given, computes list of the first
    `n` values of the exponential integral
    `E_1(x m)`.

    The exponential integral `E_1(x)` is

    .. math::

                      E_1(x) = \int_{x}^{\infty} e^{-t}/t dt

    INPUT:

    -  ``x`` - a positive real number

    -  ``n`` - (default: 0) a nonnegative integer; if
       nonzero, then return a list of values E_1(x\*m) for m =
       1,2,3,...,n. This is useful, e.g., when computing derivatives of
       L-functions.


    OUTPUT:

    -  ``float`` - if n is 0 (the default) or

    -  ``list`` - list of floats if n 0


    EXAMPLES::

        sage: exponential_integral_1(2)
        0.04890051070806112
        sage: exponential_integral_1(2,4)    # rel tol 1e-10
        [0.04890051070806112, 0.0037793524098489067, 0.00036008245216265873, 3.7665622843924751e-05]
        sage: exponential_integral_1(0)
        +Infinity

    IMPLEMENTATION: We use the PARI C-library functions eint1 and
    veceint1.

    REFERENCE:

    - See page 262, Prop 5.6.12, of Cohen's book "A Course in
      Computational Algebraic Number Theory".

    REMARKS: When called with the optional argument n, the PARI
    C-library is fast for values of n up to some bound, then very very
    slow. For example, if x=5, then the computation takes less than a
    second for n=800000, and takes "forever" for n=900000.
    """
    if isinstance(x, Expression):
        if x.is_trivial_zero():
            from sage.rings.infinity import Infinity
            return Infinity
        else:
            raise NotImplementedError("Use the symbolic exponential integral " +
                                      "function: exp_integral_e1.")
    elif not is_inexact(x): # x is exact and not an expression
        if not x: # test if exact x == 0 quickly
            from sage.rings.infinity import Infinity
            return Infinity

    # else x is in not an exact 0
    from sage.libs.pari.all import pari
    if n <= 0:
        return float(pari(x).eint1())
    else:
        return [float(z) for z in pari(x).eint1(n)]
