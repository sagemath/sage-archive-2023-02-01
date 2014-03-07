r"""
Special Functions

AUTHORS:

- David Joyner (2006-13-06): initial version

- David Joyner (2006-30-10): bug fixes to pari wrappers of Bessel
  functions, hypergeometric_U

- William Stein (2008-02): Impose some sanity checks.

- David Joyner (2008-04-23): addition of elliptic integrals

This module provides easy access to many of Maxima and PARI's
special functions.

Maxima's special functions package (which includes spherical
harmonic functions, spherical Bessel functions (of the 1st and 2nd
kind), and spherical Hankel functions (of the 1st and 2nd kind))
was written by Barton Willis of the University of Nebraska at
Kearney. It is released under the terms of the General Public
License (GPL).

Support for elliptic functions and integrals was written by Raymond
Toy. It is placed under the terms of the General Public License
(GPL) that governs the distribution of Maxima.

Next, we summarize some of the properties of the functions
implemented here.


-  Airy function The function `Ai(x)` and the related
   function `Bi(x)`, which is also called an Airy function,
   are solutions to the differential equation


   .. math::

         y'' - xy = 0,

   known as the Airy equation. They belong to the class of 'Bessel functions of
   fractional order'. The initial conditions
   `Ai(0) = (\Gamma(2/3)3^{2/3})^{-1}`,
   `Ai'(0) = -(\Gamma(1/3)3^{1/3})^{-1}` define
   `Ai(x)`. The initial conditions
   `Bi(0) = 3^{1/2}Ai(0)`, `Bi'(0) = -3^{1/2}Ai'(0)`
   define `Bi(x)`.

   They are named after the British astronomer George Biddell Airy.

-  Spherical harmonics: Laplace's equation in spherical coordinates
   is:

   .. math::

       {\frac{1}{r^2}}{\frac{\partial}{\partial r}}   \left(r^2 {\frac{\partial f}{\partial r}}\right) +   {\frac{1}{r^2}\sin\theta}{\frac{\partial}{\partial \theta}}   \left(\sin\theta {\frac{\partial f}{\partial \theta}}\right) +   {\frac{1}{r^2\sin^2\theta}}{\frac{\partial^2 f}{\partial \varphi^2}} = 0.


   Note that the spherical coordinates `\theta` and
   `\varphi` are defined here as follows: `\theta` is
   the colatitude or polar angle, ranging from
   `0\leq\theta\leq\pi` and `\varphi` the azimuth or
   longitude, ranging from `0\leq\varphi<2\pi`.

   The general solution which remains finite towards infinity is a
   linear combination of functions of the form

   .. math::

         r^{-1-\ell} \cos (m \varphi) P_\ell^m (\cos{\theta} )


   and

   .. math::

         r^{-1-\ell} \sin (m \varphi) P_\ell^m (\cos{\theta} )


   where `P_\ell^m` are the associated Legendre polynomials,
   and with integer parameters `\ell \ge 0` and `m`
   from `0` to `\ell`. Put in another way, the
   solutions with integer parameters `\ell \ge 0` and
   `- \ell\leq m\leq \ell`, can be written as linear
   combinations of:

   .. math::

         U_{\ell,m}(r,\theta , \varphi ) = r^{-1-\ell} Y_\ell^m( \theta , \varphi )


   where the functions `Y` are the spherical harmonic
   functions with parameters `\ell`, `m`, which can be
   written as:

   .. math::

         Y_\ell^m( \theta , \varphi )     = \sqrt{{\frac{(2\ell+1)}{4\pi}}{\frac{(\ell-m)!}{(\ell+m)!}}}       \cdot e^{i m \varphi } \cdot P_\ell^m ( \cos{\theta} ) .



   The spherical harmonics obey the normalisation condition


   .. math::

     \int_{\theta=0}^\pi\int_{\varphi=0}^{2\pi} Y_\ell^mY_{\ell'}^{m'*}\,d\Omega =\delta_{\ell\ell'}\delta_{mm'}\quad\quad d\Omega =\sin\theta\,d\varphi\,d\theta .



-  When solving for separable solutions of Laplace's equation in
   spherical coordinates, the radial equation has the form:

   .. math::

         x^2 \frac{d^2 y}{dx^2} + 2x \frac{dy}{dx} + [x^2 - n(n+1)]y = 0.


   The spherical Bessel functions `j_n` and `y_n`,
   are two linearly independent solutions to this equation. They are
   related to the ordinary Bessel functions `J_n` and
   `Y_n` by:

   .. math::

         j_n(x) = \sqrt{\frac{\pi}{2x}} J_{n+1/2}(x),



   .. math::

         y_n(x) = \sqrt{\frac{\pi}{2x}} Y_{n+1/2}(x)     = (-1)^{n+1} \sqrt{\frac{\pi}{2x}} J_{-n-1/2}(x).



-  For `x>0`, the confluent hypergeometric function
   `y = U(a,b,x)` is defined to be the solution to Kummer's
   differential equation


   .. math::

     xy'' + (b-x)y' - ay = 0,

   which satisfies `U(a,b,x) \sim x^{-a}`, as
   `x\rightarrow \infty`. (There is a linearly independent
   solution, called Kummer's function `M(a,b,x)`, which is not
   implemented.)

-  Jacobi elliptic functions can be thought of as generalizations
   of both ordinary and hyperbolic trig functions. There are twelve
   Jacobian elliptic functions. Each of the twelve corresponds to an
   arrow drawn from one corner of a rectangle to another.

   ::

                    n ------------------- d
                    |                     |
                    |                     |
                    |                     |
                    s ------------------- c

   Each of the corners of the rectangle are labeled, by convention, s,
   c, d and n. The rectangle is understood to be lying on the complex
   plane, so that s is at the origin, c is on the real axis, and n is
   on the imaginary axis. The twelve Jacobian elliptic functions are
   then pq(x), where p and q are one of the letters s,c,d,n.

   The Jacobian elliptic functions are then the unique
   doubly-periodic, meromorphic functions satisfying the following
   three properties:


   #. There is a simple zero at the corner p, and a simple pole at the
      corner q.

   #. The step from p to q is equal to half the period of the function
      pq(x); that is, the function pq(x) is periodic in the direction pq,
      with the period being twice the distance from p to q. Also, pq(x)
      is also periodic in the other two directions as well, with a period
      such that the distance from p to one of the other corners is a
      quarter period.

   #. If the function pq(x) is expanded in terms of x at one of the
      corners, the leading term in the expansion has a coefficient of 1.
      In other words, the leading term of the expansion of pq(x) at the
      corner p is x; the leading term of the expansion at the corner q is
      1/x, and the leading term of an expansion at the other two corners
      is 1.


   We can write

   .. math::

      pq(x)=\frac{pr(x)}{qr(x)}


   where `p`, `q`, and `r` are any of the
   letters `s`, `c`, `d`, `n`, with
   the understanding that `ss=cc=dd=nn=1`.

   Let

   .. math::

         u=\int_0^\phi \frac{d\theta} {\sqrt {1-m \sin^2 \theta}}


   Then the *Jacobi elliptic function* `sn(u)` is given by

   .. math::

        {sn}\; u = \sin \phi

   and `cn(u)` is given by

   .. math::

       {cn}\; u = \cos \phi

   and

   .. math::

     {dn}\; u = \sqrt {1-m\sin^2 \phi}.

   To emphasize the dependence on `m`, one can write
   `sn(u,m)` for example (and similarly for `cn` and
   `dn`). This is the notation used below.

   For a given `k` with `0 < k < 1` they therefore are
   solutions to the following nonlinear ordinary differential
   equations:


   -  `\mathrm{sn}\,(x;k)` solves the differential equations

      .. math::

          \frac{\mathrm{d}^2 y}{\mathrm{d}x^2} + (1+k^2) y - 2 k^2 y^3 = 0,

      and

      .. math::

         \left(\frac{\mathrm{d} y}{\mathrm{d}x}\right)^2 = (1-y^2) (1-k^2 y^2).

   -  `\mathrm{cn}\,(x;k)` solves the differential equations


      .. math::

         \frac{\mathrm{d}^2 y}{\mathrm{d}x^2} + (1-2k^2) y + 2 k^2 y^3 = 0,

      and `\left(\frac{\mathrm{d} y}{\mathrm{d}x}\right)^2 = (1-y^2) (1-k^2 + k^2 y^2)`.

   -  `\mathrm{dn}\,(x;k)` solves the differential equations

      .. math::

         \frac{\mathrm{d}^2 y}{\mathrm{d}x^2} - (2 - k^2) y + 2 y^3 = 0,

      and `\left(\frac{\mathrm{d} y}{\mathrm{d}x}\right)^2= y^2 (1 - k^2 - y^2)`.

      If `K(m)` denotes the complete elliptic integral of the
      first kind (denoted ``elliptic_kc``), the elliptic functions
      `sn (x,m)` and `cn (x,m)` have real periods
      `4K(m)`, whereas `dn (x,m)` has a period
      `2K(m)`. The limit `m\rightarrow 0` gives
      `K(0) = \pi/2` and trigonometric functions:
      `sn(x, 0) = \sin x`, `cn(x, 0) = \cos x`,
      `dn(x, 0) = 1`. The limit `m \rightarrow 1` gives
      `K(1) \rightarrow \infty` and hyperbolic functions:
      `sn(x, 1) = \tanh x`,
      `cn(x, 1) = \mbox{\rm sech} x`,
      `dn(x, 1) = \mbox{\rm sech} x`.

   -  The incomplete elliptic integrals (of the first kind, etc.) are:

      .. math::

         \begin{array}{c} \displaystyle\int_0^\phi \frac{1}{\sqrt{1 - m\sin(x)^2}}\, dx,\\ \displaystyle\int_0^\phi \sqrt{1 - m\sin(x)^2}\, dx,\\ \displaystyle\int_0^\phi \frac{\sqrt{1-mt^2}}{\sqrt(1 - t^2)}\, dx,\\ \displaystyle\int_0^\phi \frac{1}{\sqrt{1 - m\sin(x)^2\sqrt{1 - n\sin(x)^2}}}\, dx, \end{array}

      and the complete ones are obtained by taking `\phi =\pi/2`.


REFERENCES:

- Abramowitz and Stegun: Handbook of Mathematical Functions,
  http://www.math.sfu.ca/~cbm/aands/

- http://en.wikipedia.org/wiki/Airy_function

- http://en.wikipedia.org/wiki/Spherical_harmonics

- http://en.wikipedia.org/wiki/Helmholtz_equation

- http://en.wikipedia.org/wiki/Jacobi's_elliptic_functions

- A. Khare, U. Sukhatme, 'Cyclic Identities Involving
  Jacobi Elliptic Functions', Math ArXiv, math-ph/0201004

- Online Encyclopedia of Special Function
  http://algo.inria.fr/esf/index.html

TODO: Resolve weird bug in commented out code in hypergeometric_U
below.

AUTHORS:

- David Joyner and William Stein

Added 16-02-2008 (wdj): optional calls to scipy and replace all
'#random' by '...' (both at the request of William Stein)

.. warning::

   SciPy's versions are poorly documented and seem less
   accurate than the Maxima and PARI versions.
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Joyner <wdj@usna.edu>
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

from sage.plot.plot import plot
from sage.rings.real_mpfr import RealField
from sage.rings.complex_field import ComplexField
from sage.misc.sage_eval import sage_eval
from sage.rings.all import ZZ, RR, RDF
from sage.functions.other import real, imag, log_gamma
from sage.symbolic.function import BuiltinFunction
from sage.calculus.calculus import maxima

_done = False
def _init():
    """
    Internal function which checks if Maxima has loaded the
    "orthopoly" package.  All functions using this in this
    file should call this function first.

    TEST:

    The global starts ``False``::

        sage: sage.functions.special._done
        False

    Then after using one of these functions, it changes::

        sage: from sage.functions.special import airy_ai
        sage: airy_ai(1.0)
        0.135292416313
        sage: sage.functions.special._done
        True
    """
    global _done
    if _done:
        return
    maxima.eval('load("orthopoly");')
    maxima.eval('orthopoly_returns_intervals:false;')
    _done = True

def meval(x):
    """
    Returns ``x`` evaluated in Maxima, then returned to Sage.
    This is used to evaluate several of these special functions.

    TEST::

        sage: from sage.functions.special import airy_ai
        sage: airy_bi(1.0)
        1.20742359495
    """
    return maxima(x).sage()


class MaximaFunction(BuiltinFunction):
    """
    EXAMPLES::

        sage: from sage.functions.special import MaximaFunction
        sage: f = MaximaFunction("jacobi_sn")
        sage: f(1,1)
        tanh(1)
        sage: f(1/2,1/2).n()
        0.470750473655657
    """
    def __init__(self, name, nargs=2, conversions={}):
        """
        EXAMPLES::

            sage: from sage.functions.special import MaximaFunction
            sage: f = MaximaFunction("jacobi_sn")
            sage: f(1,1)
            tanh(1)
            sage: f(1/2,1/2).n()
            0.470750473655657
        """
        c = dict(maxima=name)
        c.update(conversions)
        BuiltinFunction.__init__(self, name=name, nargs=nargs,
                                   conversions=c)

    def _maxima_init_evaled_(self, *args):
        """
        Returns a string which represents this function evaluated at
        *args* in Maxima.

        EXAMPLES::

            sage: from sage.functions.special import MaximaFunction
            sage: f = MaximaFunction("jacobi_sn")
            sage: f._maxima_init_evaled_(1/2, 1/2)
            'jacobi_sn(1/2, 1/2)'

        TESTS:

        Check if complex numbers in the arguments are converted to maxima
        correctly :trac:`7557`::

            sage: t = jacobi('sn',1.2+2*I*elliptic_kc(1-.5),.5)
            sage: t._maxima_init_(maxima)
            '0.88771548861927...*%i'
            sage: t.n() # abs tol 1e-13
            0.887715488619280 - 1.79195288804672e-15*I

        """
        args_maxima = []
        for a in args:
            if isinstance(a, str):
                args_maxima.append(a)
            elif hasattr(a, '_maxima_init_'):
                args_maxima.append(a._maxima_init_())
            else:
                args_maxima.append(str(a))
        return "%s(%s)"%(self.name(), ', '.join(args_maxima))

    def _evalf_(self, *args, **kwds):
        """
        Returns a numerical approximation of this function using
        Maxima.  Currently, this is limited to 53 bits of precision.

        EXAMPLES::

            sage: from sage.functions.special import MaximaFunction
            sage: f = MaximaFunction("jacobi_sn")
            sage: f(1/2,1/2)
            jacobi_sn(1/2, 1/2)
            sage: f(1/2,1/2).n()
            0.470750473655657

        TESTS::

            sage: f(1/2,1/2).n(150)
            Traceback (most recent call last):
            ...
            NotImplementedError: jacobi_sn not implemented for precision > 53
        """
        parent = kwds['parent']
        if hasattr(parent, 'prec') and parent.prec() > 53:
            raise NotImplementedError, "%s not implemented for precision > 53"%self.name()
        _init()
        return parent(maxima("%s, numer"%self._maxima_init_evaled_(*args)))

    def _eval_(self, *args):
        """
        Returns a string which represents this function evaluated at
        *args* in Maxima.

        EXAMPLES::

            sage: from sage.functions.special import MaximaFunction
            sage: f = MaximaFunction("jacobi_sn")
            sage: f(1,1)
            tanh(1)

            sage: f._eval_(1,1)
            tanh(1)

        Here arccoth doesn't have 1 in its domain, so we just hold the expression:

            sage: elliptic_e(arccoth(1), x^2*e)
            elliptic_e(arccoth(1), x^2*e)
        """
        _init()
        try:
            s = maxima(self._maxima_init_evaled_(*args))
        except TypeError:
            return None
        if self.name() in repr(s):
            return None
        else:
            return s.sage()

from sage.misc.cachefunc import cached_function

@cached_function
def maxima_function(name):
    """
    Returns a function which is evaluated both symbolically and
    numerically via Maxima.  In particular, it returns an instance
    of :class:`MaximaFunction`.

    .. note::

       This function is cached so that duplicate copies of the same
       function are not created.

    EXAMPLES::

        sage: spherical_hankel2(2,i)
        -e
    """
    # The superclass of MaximaFunction, BuiltinFunction, assumes that there
    # will be only one symbolic function with the same name and class.
    # We create a new class for each Maxima function wrapped.
    class NewMaximaFunction(MaximaFunction):
        def __init__(self):
            """
            Constructs an object that wraps a Maxima function.

            TESTS::

                sage: spherical_hankel2(2,x)
                (-I*x^2 - 3*x + 3*I)*e^(-I*x)/x^3
            """
            MaximaFunction.__init__(self, name)

    return NewMaximaFunction()


def airy_ai(x):
   r"""
   The function `Ai(x)` and the related function `Bi(x)`,
   which is also called an *Airy function*, are
   solutions to the differential equation

   .. math::

      y'' - xy = 0,

   known as the *Airy equation*. The initial conditions
   `Ai(0) = (\Gamma(2/3)3^{2/3})^{-1}`,
   `Ai'(0) = -(\Gamma(1/3)3^{1/3})^{-1}` define `Ai(x)`.
   The initial conditions `Bi(0) = 3^{1/2}Ai(0)`,
   `Bi'(0) = -3^{1/2}Ai'(0)` define `Bi(x)`.

   They are named after the British astronomer George Biddell Airy.
   They belong to the class of "Bessel functions of fractional order".

   EXAMPLES::

       sage: airy_ai(1.0)        # last few digits are random
       0.135292416312881400
       sage: airy_bi(1.0)        # last few digits are random
       1.20742359495287099

   REFERENCE:

   - Abramowitz and Stegun: Handbook of Mathematical Functions,
     http://www.math.sfu.ca/~cbm/aands/

   - http://en.wikipedia.org/wiki/Airy_function
   """
   _init()
   return RDF(meval("airy_ai(%s)"%RDF(x)))

def airy_bi(x):
   r"""
   The function `Ai(x)` and the related function `Bi(x)`,
   which is also called an *Airy function*, are
   solutions to the differential equation

   .. math::

      y'' - xy = 0,

   known as the *Airy equation*. The initial conditions
   `Ai(0) = (\Gamma(2/3)3^{2/3})^{-1}`,
   `Ai'(0) = -(\Gamma(1/3)3^{1/3})^{-1}` define `Ai(x)`.
   The initial conditions `Bi(0) = 3^{1/2}Ai(0)`,
   `Bi'(0) = -3^{1/2}Ai'(0)` define `Bi(x)`.

   They are named after the British astronomer George Biddell Airy.
   They belong to the class of "Bessel functions of fractional order".

   EXAMPLES::

       sage: airy_ai(1)        # last few digits are random
       0.135292416312881400
       sage: airy_bi(1)        # last few digits are random
       1.20742359495287099

   REFERENCE:

   - Abramowitz and Stegun: Handbook of Mathematical Functions,
     http://www.math.sfu.ca/~cbm/aands/

   - http://en.wikipedia.org/wiki/Airy_function
   """
   _init()
   return RDF(meval("airy_bi(%s)"%RDF(x)))


def hypergeometric_U(alpha,beta,x,algorithm="pari",prec=53):
    r"""
    Default is a wrap of PARI's hyperu(alpha,beta,x) function.
    Optionally, algorithm = "scipy" can be used.

    The confluent hypergeometric function `y = U(a,b,x)` is
    defined to be the solution to Kummer's differential equation

    .. math::

             xy'' + (b-x)y' - ay = 0.

    This satisfies `U(a,b,x) \sim x^{-a}`, as
    `x\rightarrow \infty`, and is sometimes denoted
    ``x^{-a}2_F_0(a,1+a-b,-1/x)``. This is not the same as Kummer's
    `M`-hypergeometric function, denoted sometimes as
    ``_1F_1(alpha,beta,x)``, though it satisfies the same DE that
    `U` does.

    .. warning::

       In the literature, both are called "Kummer confluent
       hypergeometric" functions.

    EXAMPLES::

        sage: hypergeometric_U(1,1,1,"scipy")
        0.596347362323...
        sage: hypergeometric_U(1,1,1)
        0.59634736232319...
        sage: hypergeometric_U(1,1,1,"pari",70)
        0.59634736232319407434...
    """
    if algorithm=="scipy":
        if prec != 53:
            raise ValueError, "for the scipy algorithm the precision must be 53"
        import scipy.special
        ans = str(scipy.special.hyperu(float(alpha),float(beta),float(x)))
        ans = ans.replace("(","")
        ans = ans.replace(")","")
        ans = ans.replace("j","*I")
        return sage_eval(ans)
    elif algorithm=='pari':
        from sage.libs.pari.all import pari
        R = RealField(prec)
        return R(pari(R(alpha)).hyperu(R(beta), R(x), precision=prec))
    else:
        raise ValueError, "unknown algorithm '%s'"%algorithm

def spherical_bessel_J(n, var, algorithm="maxima"):
    r"""
    Returns the spherical Bessel function of the first kind for
    integers n >= 1.

    Reference: AS 10.1.8 page 437 and AS 10.1.15 page 439.

    EXAMPLES::

        sage: spherical_bessel_J(2,x)
        ((3/x^2 - 1)*sin(x) - 3*cos(x)/x)/x
        sage: spherical_bessel_J(1, 5.2, algorithm='scipy')
        -0.12277149950007...
        sage: spherical_bessel_J(1, 3, algorithm='scipy')
        0.345677499762355...
    """
    if algorithm=="scipy":
        from scipy.special.specfun import sphj
        return sphj(int(n), float(var))[1][-1]
    elif algorithm == 'maxima':
        _init()
        return meval("spherical_bessel_j(%s,%s)"%(ZZ(n),var))
    else:
        raise ValueError, "unknown algorithm '%s'"%algorithm

def spherical_bessel_Y(n,var, algorithm="maxima"):
    r"""
    Returns the spherical Bessel function of the second kind for
    integers n -1.

    Reference: AS 10.1.9 page 437 and AS 10.1.15 page 439.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: spherical_bessel_Y(2,x)
        -((3/x^2 - 1)*cos(x) + 3*sin(x)/x)/x
    """
    if algorithm=="scipy":
        import scipy.special
        ans = str(scipy.special.sph_yn(int(n),float(var)))
        ans = ans.replace("(","")
        ans = ans.replace(")","")
        ans = ans.replace("j","*I")
        return sage_eval(ans)
    elif algorithm == 'maxima':
        _init()
        return meval("spherical_bessel_y(%s,%s)"%(ZZ(n),var))
    else:
        raise ValueError, "unknown algorithm '%s'"%algorithm

def spherical_hankel1(n, var):
    r"""
    Returns the spherical Hankel function of the first kind for
    integers `n > -1`, written as a string. Reference: AS
    10.1.36 page 439.

    EXAMPLES::

        sage: spherical_hankel1(2, x)
        (I*x^2 - 3*x - 3*I)*e^(I*x)/x^3
    """
    return maxima_function("spherical_hankel1")(ZZ(n), var)

def spherical_hankel2(n,x):
    r"""
    Returns the spherical Hankel function of the second kind for
    integers `n > -1`, written as a string. Reference: AS 10.1.17 page
    439.

    EXAMPLES::

        sage: spherical_hankel2(2, x)
        (-I*x^2 - 3*x + 3*I)*e^(-I*x)/x^3

    Here I = sqrt(-1).
    """
    return maxima_function("spherical_hankel2")(ZZ(n), x)

def spherical_harmonic(m,n,x,y):
    r"""
    Returns the spherical Harmonic function of the second kind for
    integers `n > -1`, `|m|\leq n`. Reference:
    Merzbacher 9.64.

    EXAMPLES::

        sage: x,y = var('x,y')
        sage: spherical_harmonic(3,2,x,y)
        15/4*sqrt(7/30)*cos(x)*e^(2*I*y)*sin(x)^2/sqrt(pi)
        sage: spherical_harmonic(3,2,1,2)
        15/4*sqrt(7/30)*cos(1)*e^(4*I)*sin(1)^2/sqrt(pi)
    """
    _init()
    return meval("spherical_harmonic(%s,%s,%s,%s)"%(ZZ(m),ZZ(n),x,y))

####### elliptic functions and integrals

def elliptic_j(z):
   r"""
   Returns the elliptic modular `j`-function evaluated at `z`.

   INPUT:

   - ``z`` (complex) -- a complex number with positive imaginary part.

   OUTPUT:

   (complex) The value of `j(z)`.

   ALGORITHM:

   Calls the ``pari`` function ``ellj()``.

   AUTHOR:

   John Cremona

   EXAMPLES::

       sage: elliptic_j(CC(i))
       1728.00000000000
       sage: elliptic_j(sqrt(-2.0))
       8000.00000000000
       sage: z = ComplexField(100)(1,sqrt(11))/2
       sage: elliptic_j(z)
       -32768.000...
       sage: elliptic_j(z).real().round()
       -32768

   """
   CC = z.parent()
   from sage.rings.complex_field import is_ComplexField
   if not is_ComplexField(CC):
      CC = ComplexField()
      try:
         z = CC(z)
      except ValueError:
         raise ValueError, "elliptic_j only defined for complex arguments."
   from sage.libs.all import pari
   return CC(pari(z).ellj())

def jacobi(sym,x,m):
    r"""
    Here sym = "pq", where p,q in c,d,n,s. This returns the value of
    the Jacobi function pq(x,m), as described in the documentation for
    Sage's "special" module. There are a total of 12 functions
    described by this.

    EXAMPLES::

        sage: jacobi("sn",1,1)
        tanh(1)
        sage: jacobi("cd",1,1/2)
        jacobi_cd(1, 1/2)
        sage: RDF(jacobi("cd",1,1/2))
        0.724009721659
        sage: RDF(jacobi("cn",1,1/2)); RDF(jacobi("dn",1,1/2)); RDF(jacobi("cn",1,1/2)/jacobi("dn",1,1/2))
        0.595976567672
        0.823161001632
        0.724009721659
        sage: jsn = jacobi("sn",x,1)
        sage: P = plot(jsn,0,1)

    To view this, type P.show().
    """
    return maxima_function("jacobi_%s"%sym)(x,m)

def inverse_jacobi(sym,x,m):
    """
    Here sym = "pq", where p,q in c,d,n,s. This returns the value of
    the inverse Jacobi function `pq^{-1}(x,m)`. There are a
    total of 12 functions described by this.

    EXAMPLES::

        sage: jacobi("sn",1/2,1/2)
        jacobi_sn(1/2, 1/2)
        sage: float(jacobi("sn",1/2,1/2))
        0.4707504736556572
        sage: float(inverse_jacobi("sn",0.47,1/2))
        0.4990982313222196
        sage: float(inverse_jacobi("sn",0.4707504,0.5))
        0.4999999114665548
        sage: P = plot(inverse_jacobi('sn', x, 0.5), 0, 1, plot_points=20)

    Now to view this, just type show(P).
    """
    return maxima_function("inverse_jacobi_%s"%sym)(x,m)

#### elliptic integrals

class EllipticE(MaximaFunction):
    r"""
    This returns the value of the "incomplete elliptic integral of the
    second kind,"

    .. math::

        \int_0^\phi \sqrt{1 - m\sin(x)^2}\, dx,

    i.e., ``integrate(sqrt(1 - m*sin(x)^2), x, 0, phi)``.  Taking `\phi
    = \pi/2` gives ``elliptic_ec``.

    EXAMPLES::

        sage: z = var("z")
        sage: # this is still wrong: must be abs(sin(z)) + 2*round(z/pi)
        sage: elliptic_e(z, 1)
        2*round(z/pi) + sin(z)
        sage: elliptic_e(z, 0)
        z
        sage: elliptic_e(0.5, 0.1)
        0.498011394498832

        sage: loads(dumps(elliptic_e))
        elliptic_e
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: elliptic_e(0.5, 0.1)
            0.498011394498832
        """
        MaximaFunction.__init__(self, "elliptic_e")

elliptic_e = EllipticE()

class EllipticEC(MaximaFunction):
    """
    This returns the value of the "complete elliptic integral of the
    second kind,"

    .. math::

        \int_0^{\pi/2} \sqrt{1 - m\sin(x)^2}\, dx.

    EXAMPLES::

        sage: elliptic_ec(0.1)
        1.53075763689776
        sage: elliptic_ec(x).diff()
        1/2*(elliptic_ec(x) - elliptic_kc(x))/x

        sage: loads(dumps(elliptic_ec))
        elliptic_ec
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: elliptic_ec(0.1)
            1.53075763689776
        """
        MaximaFunction.__init__(self, "elliptic_ec", nargs=1)

    def _derivative_(self, *args, **kwds):
        """
        EXAMPLES::

            sage: elliptic_ec(x).diff()
            1/2*(elliptic_ec(x) - elliptic_kc(x))/x
        """
        diff_param = kwds['diff_param']
        assert diff_param == 0
        x = args[diff_param]
        return (elliptic_ec(x) - elliptic_kc(x))/(2*x)

elliptic_ec = EllipticEC()

class EllipticEU(MaximaFunction):
    r"""
    This returns the value of the "incomplete elliptic integral of the
    second kind,"

    .. math::

        \int_0^u \mathrm{dn}(x,m)^2\, dx = \int_0^\tau
        {\sqrt{1-m x^2}\over\sqrt{1-x^2}}\, dx.

    where `\tau = \mathrm{sn}(u, m)`.

    EXAMPLES::

        sage: elliptic_eu (0.5, 0.1)
        0.496054551286597
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: elliptic_eu (0.5, 0.1)
            0.496054551286597
        """
        MaximaFunction.__init__(self, "elliptic_eu")

elliptic_eu = EllipticEU()

class EllipticF(MaximaFunction):
    r"""
    This returns the value of the "incomplete elliptic integral of the
    first kind,"

    .. math::

        \int_0^\phi \frac{dx}{\sqrt{1 - m\sin(x)^2}},

    i.e., ``integrate(1/sqrt(1 - m*sin(x)^2), x, 0, phi)``.  Taking
    `\phi = \pi/2` gives ``elliptic_kc``.

    EXAMPLES::

        sage: z = var("z")
        sage: elliptic_f (z, 0)
        z
        sage: elliptic_f (z, 1)
        log(tan(1/4*pi + 1/2*z))
        sage: elliptic_f (0.2, 0.1)
        0.200132506747543
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: elliptic_f (0.2, 0.1)
            0.200132506747543
        """
        MaximaFunction.__init__(self, "elliptic_f")

elliptic_f = EllipticF()

class EllipticKC(MaximaFunction):
    r"""
    This returns the value of the "complete elliptic integral of the
    first kind,"

    .. math::

        \int_0^{\pi/2} \frac{dx}{\sqrt{1 - m\sin(x)^2}}.

    EXAMPLES::

        sage: elliptic_kc(0.5)
        1.85407467730137
        sage: elliptic_f(RR(pi/2), 0.5)
        1.85407467730137
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: elliptic_kc(0.5)
            1.85407467730137
            sage: elliptic_f(RR(pi/2), 0.5)
            1.85407467730137
        """
        MaximaFunction.__init__(self, "elliptic_kc", nargs=1)

elliptic_kc = EllipticKC()

class EllipticPi(MaximaFunction):
    r"""
    This returns the value of the "incomplete elliptic integral of the
    third kind,"

    .. math::

        \text{elliptic\_pi}(n, t, m) = \int_0^t \frac{dx}{(1 - n \sin(x)^2)
        \sqrt{1 - m \sin(x)^2}}.

    INPUT:

    - ``n`` -- a real number, called the "characteristic"

    - ``t`` -- a real number, called the "amplitude"

    - ``m`` -- a real number, called the "parameter"

    EXAMPLES::

        sage: N(elliptic_pi(1, pi/4, 1))
        1.14779357469632

    Compare the value computed by Maxima to the definition as a definite integral
    (using GSL)::

        sage: elliptic_pi(0.1, 0.2, 0.3)
        0.200665068220979
        sage: numerical_integral(1/(1-0.1*sin(x)^2)/sqrt(1-0.3*sin(x)^2), 0.0, 0.2)
        (0.2006650682209791, 2.227829789769088e-15)

    ALGORITHM:

    Numerical evaluation and symbolic manipulation are provided by `Maxima`_.

    REFERENCES:

    - Abramowitz and Stegun: Handbook of Mathematical Functions, section 17.7
      http://www.math.sfu.ca/~cbm/aands/
    - Elliptic Functions in `Maxima`_

    .. _`Maxima`: http://maxima.sourceforge.net/docs/manual/en/maxima_16.html#SEC91
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: elliptic_pi(0.1, 0.2, 0.3)
            0.200665068220979
        """
        MaximaFunction.__init__(self, "elliptic_pi", nargs=3)

elliptic_pi = EllipticPi()


def lngamma(t):
    r"""
    This method is deprecated, please use
    :meth:`~sage.functions.other.log_gamma` instead.

    See the :meth:`~sage.functions.other.log_gamma` function for '
    documentation and examples.

    EXAMPLES::

        sage: lngamma(RR(6))
        doctest:...: DeprecationWarning: The method lngamma() is deprecated. Use log_gamma() instead.
        See http://trac.sagemath.org/6992 for details.
        4.78749174278205
    """
    from sage.misc.superseded import deprecation
    deprecation(6992, "The method lngamma() is deprecated. Use log_gamma() instead.")
    return log_gamma(t)

def error_fcn(t):
    r"""
    The complementary error function
    `\frac{2}{\sqrt{\pi}}\int_t^\infty e^{-x^2} dx` (t belongs
    to RR).  This function is currently always
    evaluated immediately.

    EXAMPLES::

        sage: error_fcn(6)
        2.15197367124989e-17
        sage: error_fcn(RealField(100)(1/2))
        0.47950012218695346231725334611

    Note this is literally equal to `1 - erf(t)`::

        sage: 1 - error_fcn(0.5)
        0.520499877813047
        sage: erf(0.5)
        0.520499877813047
    """
    try:
        return t.erfc()
    except AttributeError:
        from sage.rings.real_mpfr import RR
        try:
            return RR(t).erfc()
        except Exception:
            raise NotImplementedError



