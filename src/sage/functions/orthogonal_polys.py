r"""
Orthogonal Polynomials

This module wraps some of the orthogonal/special functions in the
Maxima package "orthopoly". This package was written by Barton
Willis of the University of Nebraska at Kearney. It is released
under the terms of the General Public License (GPL). Send
Maxima-related bug reports and comments on this module to
willisb@unk.edu. In your report, please include Maxima and specfun
version information.


-  The Chebyshev polynomial of the first kind arises as a solution
   to the differential equation

   .. math::

         (1-x^2)\,y'' - x\,y' + n^2\,y = 0


   and those of the second kind as a solution to

   .. math::

         (1-x^2)\,y'' - 3x\,y' + n(n+2)\,y = 0.


   The Chebyshev polynomials of the first kind are defined by the
   recurrence relation

   .. math::

     T_0(x) = 1 \, T_1(x) = x \, T_{n+1}(x) = 2xT_n(x) - T_{n-1}(x). \,


   The Chebyshev polynomials of the second kind are defined by the
   recurrence relation

   .. math::

     U_0(x) = 1 \, U_1(x) = 2x \, U_{n+1}(x) = 2xU_n(x) - U_{n-1}(x). \,



   For integers `m,n`, they satisfy the orthogonality
   relations

   .. math::

     \int_{-1}^1 T_n(x)T_m(x)\,\frac{dx}{\sqrt{1-x^2}} =\left\{ \begin{matrix} 0 &: n\ne m~~~~~\\ \pi &: n=m=0\\ \pi/2 &: n=m\ne 0 \end{matrix} \right.


   and


   .. math::

     \int_{-1}^1 U_n(x)U_m(x)\sqrt{1-x^2}\,dx =\frac{\pi}{2}\delta_{m,n}.



   They are named after Pafnuty Chebyshev (alternative
   transliterations: Tchebyshef or Tschebyscheff).

-  The Hermite polynomials are defined either by

   .. math::

     H_n(x)=(-1)^n e^{x^2/2}\frac{d^n}{dx^n}e^{-x^2/2}


   (the "probabilists' Hermite polynomials"), or by


   .. math::

     H_n(x)=(-1)^n e^{x^2}\frac{d^n}{dx^n}e^{-x^2}


   (the "physicists' Hermite polynomials"). Sage (via Maxima)
   implements the latter flavor. These satisfy the orthogonality
   relation

   .. math::

     \int_{-\infty}^\infty H_n(x)H_m(x)\,e^{-x^2}\,dx ={n!2^n}{\sqrt{\pi}}\delta_{nm}



   They are named in honor of Charles Hermite.

-  Each *Legendre polynomial* `P_n(x)` is an `n`-th degree polynomial.
   It may be expressed using Rodrigues' formula:

   .. math::

      P_n(x) = (2^n n!)^{-1} {\frac{d^n}{dx^n} } \left[ (x^2 -1)^n \right].

   These are solutions to Legendre's differential equation:

   .. math::

      {\frac{d}{dx}} \left[ (1-x^2) {\frac{d}{dx}} P(x) \right] + n(n+1)P(x) = 0.

   and satisfy the orthogonality relation

   .. math::

      \int_{-1}^{1} P_m(x) P_n(x)\,dx = {\frac{2}{2n + 1}} \delta_{mn}

   The *Legendre function of the second kind* `Q_n(x)` is another
   (linearly independent) solution to the Legendre differential equation.
   It is not an "orthogonal polynomial" however.

   The associated Legendre functions of the first kind
   `P_\ell^m(x)` can be given in terms of the "usual"
   Legendre polynomials by

   .. math::

     \begin{array}{ll} P_\ell^m(x)    &=  (-1)^m(1-x^2)^{m/2}\frac{d^m}{dx^m}P_\ell(x) \\ &=  \frac{(-1)^m}{2^\ell \ell!} (1-x^2)^{m/2}\frac{d^{\ell+m}}{dx^{\ell+m}}(x^2-1)^\ell. \end{array}


   Assuming `0 \le m \le \ell`, they satisfy the orthogonality
   relation:

   .. math::

      \int_{-1}^{1} P_k ^{(m)} P_\ell ^{(m)} dx  = \frac{2 (\ell+m)!}{(2\ell+1)(\ell-m)!}\ \delta _{k,\ell},


   where `\delta _{k,\ell}` is the Kronecker delta.

   The associated Legendre functions of the second kind
   `Q_\ell^m(x)` can be given in terms of the "usual"
   Legendre polynomials by


   .. math::

     Q_\ell^m(x)   =  (-1)^m(1-x^2)^{m/2}\frac{d^m}{dx^m}Q_\ell(x).



   They are named after Adrien-Marie Legendre.

-  Laguerre polynomials may be defined by the Rodrigues formula

   .. math::

      L_n(x)=\frac{e^x}{n!}\frac{d^n}{dx^n}\left(e^{-x} x^n\right).


   They are solutions of Laguerre's equation:


   .. math::

      x\,y'' + (1 - x)\,y' + n\,y = 0\,

   and satisfy the orthogonality relation


   .. math::

      \int_0^\infty L_m(x) L_n(x) e^{-x}\,dx = \delta_{mn}.



   The generalized Laguerre polynomials may be defined by the
   Rodrigues formula:


   .. math::

       L_n^{(\alpha)}(x)   = {\frac{x^{-\alpha} e^x}{n!}}{\frac{d^n}{dx^n}} \left(e^{-x} x^{n+\alpha}\right) .


   (These are also sometimes called the associated Laguerre
   polynomials.) The simple Laguerre polynomials are recovered from
   the generalized polynomials by setting `\alpha =0`.

   They are named after Edmond Laguerre.

-  Jacobi polynomials are a class of orthogonal polynomials. They
   are obtained from hypergeometric series in cases where the series
   is in fact finite:

   .. math::

     P_n^{(\alpha,\beta)}(z) =\frac{(\alpha+1)_n}{n!} \,_2F_1\left(-n,1+\alpha+\beta+n;\alpha+1;\frac{1-z}{2}\right) ,


   where `()_n` is Pochhammer's symbol (for the rising
   factorial), (Abramowitz and Stegun p561.) and thus have the
   explicit expression


   .. math::

     P_n^{(\alpha,\beta)} (z) = \frac{\Gamma (\alpha+n+1)}{n!\Gamma (\alpha+\beta+n+1)} \sum_{m=0}^n {n\choose m} \frac{\Gamma (\alpha + \beta + n + m + 1)}{\Gamma (\alpha + m + 1)} \left(\frac{z-1}{2}\right)^m .



   They are named after Carl Jacobi.

-  Ultraspherical or Gegenbauer polynomials are given in terms of
   the Jacobi polynomials `P_n^{(\alpha,\beta)}(x)` with
   `\alpha=\beta=a-1/2` by


   .. math::

     C_n^{(a)}(x)= \frac{\Gamma(a+1/2)}{\Gamma(2a)}\frac{\Gamma(n+2a)}{\Gamma(n+a+1/2)} P_n^{(a-1/2,a-1/2)}(x).


   They satisfy the orthogonality relation

   .. math::

     \int_{-1}^1(1-x^2)^{a-1/2}C_m^{(a)}(x)C_n^{(a)}(x)\, dx =\delta_{mn}2^{1-2a}\pi \frac{\Gamma(n+2a)}{(n+a)\Gamma^2(a)\Gamma(n+1)} ,


   for `a>-1/2`. They are obtained from hypergeometric series
   in cases where the series is in fact finite:


   .. math::

     C_n^{(a)}(z) =\frac{(2a)^{\underline{n}}}{n!} \,_2F_1\left(-n,2a+n;a+\frac{1}{2};\frac{1-z}{2}\right)


   where `\underline{n}` is the falling factorial. (See
   Abramowitz and Stegun p561)

   They are named for Leopold Gegenbauer (1849-1903).


For completeness, the Pochhammer symbol, introduced by Leo August
Pochhammer, `(x)_n`, is used in the theory of special
functions to represent the "rising factorial" or "upper factorial"

.. math::

         (x)_n=x(x+1)(x+2)\cdots(x+n-1)=\frac{(x+n-1)!}{(x-1)!}.


On the other hand, the "falling factorial" or "lower factorial" is

.. math::

     x^{\underline{n}}=\frac{x!}{(x-n)!} ,


in the notation of Ronald L. Graham, Donald E. Knuth and Oren
Patashnik in their book Concrete Mathematics.

.. note::

   The first call of any of these will usually cost a bit extra
   (it loads "specfun", but I'm not sure if that is the real reason).
   The next call is usually faster but not always.

TODO: Implement associated Legendre polynomials and Zernike
polynomials. (Neither is in Maxima.)
http://en.wikipedia.org/wiki/Associated_Legendre_polynomials
http://en.wikipedia.org/wiki/Zernike_polynomials

REFERENCES:

-  Abramowitz and Stegun: Handbook of Mathematical Functions,
   http://www.math.sfu.ca/ cbm/aands/

-  http://en.wikipedia.org/wiki/Chebyshev_polynomials

-  http://en.wikipedia.org/wiki/Legendre_polynomials

-  http://en.wikipedia.org/wiki/Hermite_polynomials

-  http://mathworld.wolfram.com/GegenbauerPolynomial.html

-  http://en.wikipedia.org/wiki/Jacobi_polynomials

-  http://en.wikipedia.org/wiki/Laguerre_polynomia

-  http://en.wikipedia.org/wiki/Associated_Legendre_polynomials


AUTHORS:

- David Joyner (2006-06)
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

import copy
import sage.plot.plot
import sage.interfaces.all
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import RationalField
from sage.rings.real_mpfr import RealField
from sage.misc.sage_eval import sage_eval
from sage.rings.all import QQ, ZZ, CDF, RDF
import sage.rings.commutative_ring as commutative_ring
import sage.rings.ring as ring
from sage.calculus.calculus import maxima

_done = False
def _init():
    global _done
    if _done:
        return
    maxima.eval('load("orthopoly");')
    # TODO -- make it possible to use the intervals returned
    # instead of just discarding this info!
    maxima.eval('orthopoly_returns_intervals:false;')
    _done = True


def chebyshev_T(n,x):
    """
    Returns the Chebyshev function of the first kind for integers
    `n>-1`.

    REFERENCE:

    - AS 22.5.31 page 778 and AS 6.1.22 page 256.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: chebyshev_T(2,x)
        2*x^2 - 1
    """
    _init()
    return sage_eval(maxima.eval('chebyshev_t(%s,x)'%ZZ(n)), locals={'x':x})

def chebyshev_U(n,x):
    """
    Returns the Chebyshev function of the second kind for integers `n>-1`.

    REFERENCE:

    - AS, 22.8.3 page 783 and AS 6.1.22 page 256.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: chebyshev_U(2,x)
        4*x^2 - 1
    """
    _init()
    return sage_eval(maxima.eval('chebyshev_u(%s,x)'%ZZ(n)), locals={'x':x})

def gen_laguerre(n,a,x):
    """
    Returns the generalized Laguerre polynomial for integers `n > -1`.
    Typically, a = 1/2 or a = -1/2.

    REFERENCE:

    - table on page 789 in AS.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: gen_laguerre(2,1,x)
        1/2*x^2 - 3*x + 3
        sage: gen_laguerre(2,1/2,x)
        1/2*x^2 - 5/2*x + 15/8
        sage: gen_laguerre(2,-1/2,x)
        1/2*x^2 - 3/2*x + 3/8
        sage: gen_laguerre(2,0,x)
        1/2*x^2 - 2*x + 1
        sage: gen_laguerre(3,0,x)
        -1/6*x^3 + 3/2*x^2 - 3*x + 1
    """
    from sage.functions.all import sqrt
    _init()
    return sage_eval(maxima.eval('gen_laguerre(%s,%s,x)'%(ZZ(n),a)), locals={'x':x})

def gen_legendre_P(n,m,x):
    r"""
    Returns the generalized (or associated) Legendre function of the
    first kind for integers `n > -1, m > -1`.

    The awkward code for when m is odd and 1 results from the fact that
    Maxima is happy with, for example, `(1 - t^2)^3/2`, but
    Sage is not. For these cases the function is computed from the
    (m-1)-case using one of the recursions satisfied by the Legendre
    functions.

    REFERENCE:

    - Gradshteyn and Ryzhik 8.706 page 1000.

    EXAMPLES::

        sage: P.<t> = QQ[]
        sage: gen_legendre_P(2, 0, t)
        3/2*t^2 - 1/2
        sage: gen_legendre_P(2, 0, t) == legendre_P(2, t)
        True
        sage: gen_legendre_P(3, 1, t)
        -3/2*sqrt(-t^2 + 1)*(5*t^2 - 1)
        sage: gen_legendre_P(4, 3, t)
        (105*t^3 - 105*t)*sqrt(-t^2 + 1)
        sage: gen_legendre_P(7, 3, I).expand()
        -16695*sqrt(2)
        sage: gen_legendre_P(4, 1, 2.5)
        -583.562373654533*I
    """
    from sage.functions.all import sqrt
    _init()
    if m.mod(2).is_zero() or m.is_one():
        return sage_eval(maxima.eval('assoc_legendre_p(%s,%s,x)'%(ZZ(n),ZZ(m))), locals={'x':x})
    else:
        return sqrt(1-x**2)*(((n-m+1)*x*gen_legendre_P(n,m-1,x)-(n+m-1)*gen_legendre_P(n-1,m-1,x))/(1-x**2))

def gen_legendre_Q(n,m,x):
    """
    Returns the generalized (or associated) Legendre function of the
    second kind for integers `n>-1`, `m>-1`.

    Maxima restricts m = n. Hence the cases m n are computed using the
    same recursion used for gen_legendre_P(n,m,x) when m is odd and
    1.

    EXAMPLES::

        sage: P.<t> = QQ[]
        sage: gen_legendre_Q(2,0,t)
        3/4*t^2*log((-t - 1)/(t - 1)) - 3/2*t - 1/4*log((-t - 1)/(t - 1))
        sage: gen_legendre_Q(2,0,t) - legendre_Q(2, t)
        0
        sage: gen_legendre_Q(3,1,0.5)
        2.49185259170895
        sage: gen_legendre_Q(0, 1, x)
        -1/sqrt(-x^2 + 1)
        sage: gen_legendre_Q(2, 4, x).factor()
        48*x/((x - 1)^2*(x + 1)^2)
    """
    from sage.functions.all import sqrt
    if m <= n:
        _init()
        return sage_eval(maxima.eval('assoc_legendre_q(%s,%s,x)'%(ZZ(n),ZZ(m))), locals={'x':x})
    if m == n + 1 or n == 0:
        if m.mod(2).is_zero():
            denom = (1 - x**2)**(m/2)
        else:
            denom = sqrt(1 - x**2)*(1 - x**2)**((m-1)/2)
        if m == n + 1:
            return (-1)**m*(m-1).factorial()*2**n/denom
        else:
            return (-1)**m*(m-1).factorial()*((x+1)**m - (x-1)**m)/(2*denom)
    else:
        return ((n-m+1)*x*gen_legendre_Q(n,m-1,x)-(n+m-1)*gen_legendre_Q(n-1,m-1,x))/sqrt(1-x**2)

def hermite(n,x):
    """
    Returns the Hermite polynomial for integers `n > -1`.

    REFERENCE:

    - AS 22.5.40 and 22.5.41, page 779.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: hermite(2,x)
        4*x^2 - 2
        sage: hermite(3,x)
        8*x^3 - 12*x
        sage: hermite(3,2)
        40
        sage: S.<y> = PolynomialRing(RR)
        sage: hermite(3,y)
        8.00000000000000*y^3 - 12.0000000000000*y
        sage: R.<x,y> = QQ[]
        sage: hermite(3,y^2)
        8*y^6 - 12*y^2
        sage: w = var('w')
        sage: hermite(3,2*w)
        8*(8*w^2 - 3)*w
    """
    _init()
    return sage_eval(maxima.eval('hermite(%s,x)'%ZZ(n)), locals={'x':x})

def jacobi_P(n,a,b,x):
    r"""
    Returns the Jacobi polynomial `P_n^{(a,b)}(x)` for
    integers `n > -1` and a and b symbolic or `a > -1`
    and `b > -1`. The Jacobi polynomials are actually defined
    for all a and b. However, the Jacobi polynomial weight
    `(1-x)^a(1+x)^b` isn't integrable for `a \leq -1`
    or `b \leq -1`.

    REFERENCE:

    - table on page 789 in AS.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: jacobi_P(2,0,0,x)
        3/2*x^2 - 1/2
        sage: jacobi_P(2,1,2,1.2)        # random output of low order bits
        5.009999999999998
    """
    _init()
    return sage_eval(maxima.eval('jacobi_p(%s,%s,%s,x)'%(ZZ(n),a,b)), locals={'x':x})

def laguerre(n,x):
    """
    Returns the Laguerre polynomial for integers `n > -1`.

    REFERENCE:

    - AS 22.5.16, page 778 and AS page 789.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: laguerre(2,x)
        1/2*x^2 - 2*x + 1
        sage: laguerre(3,x)
        -1/6*x^3 + 3/2*x^2 - 3*x + 1
        sage: laguerre(2,2)
        -1
    """
    _init()
    return sage_eval(maxima.eval('laguerre(%s,x)'%ZZ(n)), locals={'x':x})

def legendre_P(n,x):
    """
    Returns the Legendre polynomial of the first kind for integers
    `n > -1`.

    REFERENCE:

    - AS 22.5.35 page 779.

    EXAMPLES::

        sage: P.<t> = QQ[]
        sage: legendre_P(2,t)
        3/2*t^2 - 1/2
        sage: legendre_P(3, 1.1)
        1.67750000000000
        sage: legendre_P(3, 1 + I)
        7/2*I - 13/2
        sage: legendre_P(3, MatrixSpace(ZZ, 2)([1, 2, -4, 7]))
        [-179  242]
        [-484  547]
        sage: legendre_P(3, GF(11)(5))
        8
    """
    _init()
    return sage_eval(maxima.eval('legendre_p(%s,x)'%ZZ(n)), locals={'x':x})

def legendre_Q(n,x):
    """
    Returns the Legendre function of the second kind for integers
    `n>-1`.

    Computed using Maxima.

    EXAMPLES::

        sage: P.<t> = QQ[]
        sage: legendre_Q(2, t)
        3/4*t^2*log((-t - 1)/(t - 1)) - 3/2*t - 1/4*log((-t - 1)/(t - 1))
        sage: legendre_Q(3, 0.5)
        -0.198654771479482
        sage: legendre_Q(4, 2)
        443/16*I*pi + 443/16*log(3) - 365/12
        sage: legendre_Q(4, 2.0)
        0.00116107583162324 + 86.9828465962674*I
    """
    _init()
    return sage_eval(maxima.eval('legendre_q(%s,x)'%ZZ(n)), locals={'x':x})

def ultraspherical(n,a,x):
    """
    Returns the ultraspherical (or Gegenbauer) polynomial for integers
    `n > -1`.

    Computed using Maxima.

    REFERENCE:

    - AS 22.5.27

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: ultraspherical(2,3/2,x)
        15/2*x^2 - 3/2
        sage: ultraspherical(2,1/2,x)
        3/2*x^2 - 1/2
        sage: ultraspherical(1,1,x)
        2*x
        sage: t = PolynomialRing(RationalField(),"t").gen()
        sage: gegenbauer(3,2,t)
        32*t^3 - 12*t
    """
    _init()
    return sage_eval(maxima.eval('ultraspherical(%s,%s,x)'%(ZZ(n),a)), locals={'x':x})

gegenbauer = ultraspherical
