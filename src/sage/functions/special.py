r"""
Special Functions

AUTHORS:
   -- David Joyner (2006-13-06), initial version
   -- David Joyner (2006-30-10), bug fixes to pari wrappers of Bessel
                                 functions, hypergeometric_U
   -- William Stein (2008-02): Impose some sanity checks.
   -- David Joyner (2008-04-23), addition of elliptic integrals

This module provides easy access to many of Maxima and PARI's
special functions.

Maxima's special functions package (which includes spherical harmonic
functions, spherical Bessel functions (of the 1st and 2nd kind), and
spherical Hankel functions (of the 1st and 2nd kind)) was written by
Barton Willis of the University of Nebraska at Kearney.  It is
released under the terms of the General Public License (GPL).

Support for elliptic functions and integrals was written by
Raymond Toy. It is placed under the terms of the General
Public License (GPL) that governs the distribution of Maxima.

The (usual) Bessel functions and Airy functions are part of the
standard Maxima package. Some Bessel functions also are implemented in
Pari. (Caution: The Pari versions are sometimes different than the
Maxima version.) For example, the K-Bessel function $K_\nu (z)$ can be
computed using either Maxima or Pari, depending on an optional
variable you pass to bessel_K.

Next, we summarize some of the properties of the functions implemented
here.

\begin{itemize}
\item
{\it Bessel functions}, first defined by the Swiss mathematician
Daniel Bernoulli and named after Friedrich Bessel, are canonical
solutions y(x) of {\it Bessel's differential equation}:

\[
    x^2 \frac{d^2 y}{dx^2} + x \frac{dy}{dx} + (x^2 - \alpha^2)y = 0,
\]
for an arbitrary real number $\alpha$ (the order).

\item
Another important formulation of the two linearly independent
solutions to Bessel's equation are the {\it Hankel functions}
$H_\alpha^{(1)}(x)$ and $H_\alpha^(2)(x)$, defined by:

\[
    H_\alpha^{(1)}(x) = J_\alpha(x) + i Y_\alpha(x)
\]
\[
    H_\alpha^{(2)}(x) = J_\alpha(x) - i Y_\alpha(x)
\]
where $i$ is the imaginary unit (and $J_*$ and $Y_*$ are the
usual J- and Y-Bessel functions). These linear combinations are
also known as {\it Bessel functions of the third kind}; they are two
linearly independent solutions of Bessel's differential equation.
They are named for Hermann Hankel.

\item
{\it Airy function} The function $Ai(x)$ and the related
function $Bi(x)$, which is also called an {\it Airy function}, are
solutions to the differential equation

\[
    y'' - xy = 0,
\]
known as the {\it Airy equation}. They belong to the class of
"Bessel functions of fractional order". The initial conditions
$Ai(0) = (\Gamma(2/3)3^{2/3})^{-1}$,
$Ai'(0) = -(\Gamma(1/3)3^{1/3})^{-1}$ define $Ai(x)$.
The initial conditions $Bi(0) = 3^{1/2}Ai(0)$,
$Bi'(0) = -3^{1/2}Ai'(0)$ define $Bi(x)$.

They are named after the British astronomer George Biddell Airy.

\item
Spherical harmonics:
Laplace's equation in spherical coordinates is:
\[
  {\frac{1}{r^2}}{\frac{\partial}{\partial r}}
  \left(r^2 {\frac{\partial f}{\partial r}}\right) +
  {\frac{1}{r^2}\sin\theta}{\frac{\partial}{\partial \theta}}
  \left(\sin\theta {\frac{\partial f}{\partial \theta}}\right) +
  {\frac{1}{r^2\sin^2\theta}}{\frac{\partial^2 f}{\partial \varphi^2}} = 0.
\]
Note that the spherical coordinates $\theta$
and $\varphi$ are defined here as follows:
$\theta$ is the colatitude or polar angle, ranging from
$0\leq\theta\leq\pi$ and $\varphi$ the azimuth or longitude,
ranging from $0\leq\varphi<2\pi$.

The general solution which remains finite towards infinity
is a linear combination of functions of the form
\[
    r^{-1-\ell} \cos (m \varphi) P_\ell^m (\cos{\theta} )
\]
and
\[
    r^{-1-\ell} \sin (m \varphi) P_\ell^m (\cos{\theta} )
\]
where $P_\ell^m$ are the associated {\it Legendre polynomials}, and
with integer parameters $\ell \ge 0$ and $m$ from $0$ to $\ell$.
Put in another way, the solutions with integer parameters
$\ell \ge 0$ and $- \ell\leq m\leq \ell$, can be written as
linear combinations of:
\[
    U_{\ell,m}(r,\theta , \varphi ) = r^{-1-\ell} Y_\ell^m( \theta , \varphi )
\]
where the functions $Y$ are the {\it spherical harmonic functions}
with parameters $\ell$, $m$, which can be written as:
\[
    Y_\ell^m( \theta , \varphi )
    = \sqrt{{\frac{(2\ell+1)}{4\pi}}{\frac{(\ell-m)!}{(\ell+m)!}}}
      \cdot e^{i m \varphi } \cdot P_\ell^m ( \cos{\theta} ) .
\]

The spherical harmonics obey the normalisation condition

\[
\int_{\theta=0}^\pi\int_{\varphi=0}^{2\pi}
Y_\ell^mY_{\ell'}^{m'*}\,d\Omega
=\delta_{\ell\ell'}\delta_{mm'}\quad\quad d\Omega
=\sin\theta\,d\varphi\,d\theta .
\]

\item
When solving for separable solutions of Laplace's equation in
spherical coordinates, the radial equation has the form:
\[
    x^2 \frac{d^2 y}{dx^2} + 2x \frac{dy}{dx} + [x^2 - n(n+1)]y = 0.
\]
The {\it spherical Bessel functions} $j_n$ and $y_n$,
are two linearly independent solutions to this equation.
They are related to the ordinary Bessel functions $J_n$ and $Y_n$ by:
\[
    j_n(x) = \sqrt{\frac{\pi}{2x}} J_{n+1/2}(x),
\]
\[
    y_n(x) = \sqrt{\frac{\pi}{2x}} Y_{n+1/2}(x)
    = (-1)^{n+1} \sqrt{\frac{\pi}{2x}} J_{-n-1/2}(x).
\]

\item
For $x>0$, the confluent hypergeometric function
$y = U(a,b,x)$ is defined to be the solution to Kummer's
differential equation

\[
xy'' + (b-x)y' - ay = 0,
\]
which satisfies $U(a,b,x) \sim x^{-a}$, as $x\rightarrow \infty$.
(There is a linearly independent solution, called Kummer's
function $M(a,b,x)$, which is not implemented.)


\item
Jacobi elliptic functions can be thought of as generalizations
of both ordinary and hyperbolic trig functions.
There are twelve Jacobian elliptic functions. Each of the twelve
corresponds to an arrow drawn from one corner of a rectangle to
another.
\begin{verbatim}
             n ------------------- d
             |                     |
             |                     |
             |                     |
             s ------------------- c
\end{verbatim}

Each of the corners of the rectangle are labeled, by convention,
s, c, d and n. The rectangle is understood to be lying on the complex
plane, so that s is at the origin, c is on the real axis,
and n is on the imaginary axis.
The twelve Jacobian elliptic functions are then pq(x), where p and q
are one of the letters s,c,d,n.

The {\it Jacobian elliptic functions} are then the unique
doubly-periodic, meromorphic functions satisfying the following
three properties:

\begin{enumerate}
\item
There is a simple zero at the corner p, and a simple pole
at the corner q.
\item
The step from p to q is equal to half the period of the function pq(x);
that is, the function pq(x) is periodic in the direction pq,
with the period being twice the distance from p to q.
Also, pq(x) is also periodic in the other two directions as well,
with a period such that the distance from p to one of the other corners
is a quarter period.
\item
If the function pq(x) is expanded in terms of x at one of the corners,
the leading term in the expansion has a coefficient of 1.
In other words, the leading term of the expansion of pq(x) at the
corner p is x; the leading term of the expansion at the corner
q is 1/x, and the leading term of an expansion at the other two
corners is 1.
\end{enumerate}

We can write
\[
pq(x)=\frac{pr(x)}{qr(x)}
\]
where $p$, $q$, and $r$ are any of the letters $s$, $c$, $d$, $n$,
with the understanding that $ss=cc=dd=nn=1$.

Let
\[
    u=\int_0^\phi \frac{d\theta} {\sqrt {1-m \sin^2 \theta}}
\]
Then the \emph{Jacobi elliptic function} $sn(u)$ is given by
\[
   {sn}\; u = \sin \phi
\]
and $cn(u)$ is given by

\[
  {cn}\; u = \cos \phi
\]
and
\[
{dn}\; u = \sqrt {1-m\sin^2 \phi}.
\]
To emphasize the dependence on $m$, one can write
$sn(u,m)$ for example (and similarly for $cn$ and $dn$).
This is the notation used below.

For a given $k$ with $0 < k < 1$ they therefore are solutions to
the following nonlinear ordinary differential equations:

\begin{itemize}
\item
$\mathrm{sn}\,(x;k)$ solves the differential equations
\[
\frac{\mathrm{d}^2 y}{\mathrm{d}x^2} + (1+k^2) y - 2 k^2 y^3 = 0,
\]
and $\left(\frac{\mathrm{d} y}{\mathrm{d}x}\right)^2 = (1-y^2) (1-k^2 y^2)$.

\item
$\mathrm{cn}\,(x;k)$ solves the differential equations

\[
\frac{\mathrm{d}^2 y}{\mathrm{d}x^2} + (1-2k^2) y + 2 k^2 y^3 = 0,
\]
and $\left(\frac{\mathrm{d} y}{\mathrm{d}x}\right)^2
= (1-y^2) (1-k^2 + k^2 y^2)$.

\item
$\mathrm{dn}\,(x;k)$ solves the differential equations

\[
\frac{\mathrm{d}^2 y}{\mathrm{d}x^2} - (2 - k^2) y + 2 y^3 = 0,
\]
and $\left(\frac{\mathrm{d} y}{\mathrm{d}x}\right)^2
= y^2 (1 - k^2 - y^2)$.


If $K(m)$ denotes the {\it complete elliptic integral of the first kind}
(denoted \verb+elliptic_kc+), the elliptic functions $sn (x,m)$ and
$cn (x,m)$ have real periods $4K(m)$, whereas $dn (x,m)$ has a period
$2K(m)$. The limit $m\rightarrow 0$ gives
$K(0) = \pi/2$ and trigonometric functions: $sn(x, 0) = \sin x$,
$cn(x, 0) = \cos x$, $dn(x, 0) = 1$. The limit $m \rightarrow 1$
gives $K(1) \rightarrow \infty$ and hyperbolic functions:
$sn(x, 1) = \tanh x$, $cn(x, 1) = \mbox{\rm sech} x$,
$dn(x, 1) = \mbox{\rm sech} x$.

\item
The {\it incomplete elliptic integrals} (of the first kind,
etc.) are:
\[
\begin{array}{c}
\displaystyle\int_0^\phi \frac{1}{\sqrt{1 - m\sin(x)^2}}\, dx,\\
\displaystyle\int_0^\phi \sqrt{1 - m\sin(x)^2}\, dx,\\
\displaystyle\int_0^\phi \frac{\sqrt{1-mt^2}}{\sqrt(1 - t^2)}\, dx,\\
\displaystyle\int_0^\phi \frac{1}{\sqrt{1 - m\sin(x)^2\sqrt{1 - n\sin(x)^2}}}\, dx,
\end{array}
\]
and the {\it complete} ones are obtained by taking $\phi =\pi/2$.


\end{itemize}
\end{itemize}

\begin{verbatim}
Methods implemented:
    * Bessel functions and Airy functions
    * spherical harmonic functions
    * spherical Bessel functions (of the 1st and 2nd kind)
    * spherical Hankel functions (of the 1st and 2nd kind)
    * Jacobi elliptic functions
    * complete/incomplete elliptic integrals
    * hyperbolic trig functions (for completeness, since
      they are special cases of elliptic functions)
    * Kummer confluent $U$-hypergeometric functions.

REFERENCE:
    * Abramowitz and Stegun: Handbook of Mathematical Functions,
      http://www.math.sfu.ca/~cbm/aands/
    * http://en.wikipedia.org/wiki/Bessel_function
    * http://en.wikipedia.org/wiki/Airy_function
    * http://en.wikipedia.org/wiki/Spherical_harmonics
    * http://en.wikipedia.org/wiki/Helmholtz_equation
    * http://en.wikipedia.org/wiki/Jacobi's_elliptic_functions
    * A. Khare, U. Sukhatme, "Cyclic Identities Involving
      Jacobi Elliptic Functions", Math ArXiv, math-ph/0201004
    * Online Encyclopedia of Special Function
      http://algo.inria.fr/esf/index.html
\end{verbatim}

TODO:
    Resolve weird bug in commented out code in hypergeometric_U below.

AUTHORS:
    David Joyner and William Stein

Added 16-02-2008 (wdj): optional calls to scipy and replace all
"\#random" by "..." (both at the request of William Stein)

WARNING:
   SciPy's versions are poorly documented and seem less
accurate than the Maxima and Pari versions.

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
from sage.plot.plot import plot
import sage.interfaces.all
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import RationalField
from sage.rings.real_mpfr import RealField
from sage.rings.complex_field import ComplexField
from sage.misc.sage_eval import sage_eval
from sage.rings.all import ZZ, QQ, RR
import sage.rings.commutative_ring as commutative_ring
import sage.rings.ring as ring
from sage.misc.functional import real, imag

from sage.interfaces.maxima import maxima

def meval(x):
    from sage.calculus.calculus import symbolic_expression_from_maxima_element
    return symbolic_expression_from_maxima_element(maxima(x))

from functions import *

_done = False
def _init():
    global _done
    if _done:
        return
    maxima.eval('load("orthopoly");')
    maxima.eval('orthopoly_returns_intervals:false;')
    _done = True

def bessel_I(nu,z,algorithm = "pari",prec=53):
    r"""
    Implements the "I-Bessel function", or
    "modified Bessel function, 1st kind", with
    index (or "order") nu and argument z.

    INPUT:
        nu -- a real (or complex, for pari) number
        z  -- a real (positive)
        algorithm - "pari" or "maxima" or "scipy"
        prec - real precision (for Pari only)

    DEFINITION:
    \begin{verbatim}
    Maxima:
                     inf
                    ====   - nu - 2 k  nu + 2 k
                    \     2          z
                     >    -------------------
                    /     k! Gamma(nu + k + 1)
                    ====
                    k = 0

    Pari:

                     inf
                    ====   - 2 k  2 k
                    \     2      z    Gamma(nu + 1)
                     >    -----------------------
                    /       k! Gamma(nu + k + 1)
                    ====
                    k = 0

    \end{verbatim}
    Sometimes \code{bessel_I(nu,z)} is denoted \code{I_nu(z)} in the
    literature.

    WARNING:
       In Maxima (the manual says) i0 is deprecated but \code{bessel_i(0,*)}
       is broken. (Was fixed in recent CVS patch though.)

    EXAMPLES:
        sage: bessel_I(1,1,"pari",500)
        0.565159103992485027207696027609863307328899621621092009480294489479255640964371134092664997766814410064677886055526302676857637684917179812041131208121
        sage: bessel_I(1,1)
        0.565159103992485
        sage: bessel_I(2,1.1,"maxima")
        0.16708949925104...
        sage: bessel_I(0,1.1,"maxima")
        1.32616018371265...
        sage: bessel_I(0,1,"maxima")
        1.26606587775200...
        sage: bessel_I(1,1,"scipy")
        0.565159103992...

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
        return sage_eval(ans)
    elif algorithm == "maxima":
        if prec != 53:
            raise ValueError, "for the maxima algorithm the precision must be 53"
        return sage_eval(maxima.eval("bessel_i(%s,%s)"%(float(nu),float(z))))
    else:
        raise ValueError, "unknown algorithm '%s'"%algorithm

def bessel_J(nu,z,algorithm="pari",prec=53):
    r"""
    Return value of the "J-Bessel function", or
    "Bessel function, 1st kind", with
    index (or "order") nu and argument z.

    \begin{verbatim}
    Defn:
    Maxima:
                     inf
                    ====          - nu - 2 k  nu + 2 k
                    \     (-1)^k 2           z
                     >    -------------------------
                    /        k! Gamma(nu + k + 1)
                    ====
                    k = 0

    Pari:

                     inf
                    ====          - 2k    2k
                    \     (-1)^k 2      z    Gamma(nu + 1)
                     >    ----------------------------
                    /         k! Gamma(nu + k + 1)
                    ====
                    k = 0
    \end{verbatim}


    Sometimes bessel_J(nu,z) is denoted J_nu(z) in the
    literature.

    WARNING:
        Inaccurate for small values of z.

    EXAMPLES:
        sage: bessel_J(2,1.1)
        0.136564153956658
        sage: bessel_J(0,1.1)
        0.719622018527511
        sage: bessel_J(0,1)
        0.765197686557967
        sage: bessel_J(0,0)
        1.00000000000000
        sage: bessel_J(0.1,0.1)
        0.777264368097005

    We check consistency of PARI and Maxima:
        sage: n(bessel_J(3,10,"maxima"))
        0.0583793793051...
        sage: n(bessel_J(3,10,"pari"))
        0.0583793793051868
        sage: bessel_J(3,10,"scipy")
        0.0583793793052...
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
        return sage_eval(ans)
    elif algorithm == "maxima":
        if prec != 53:
            raise ValueError, "for the maxima algorithm the precision must be 53"
        return meval("bessel_j(%s,%s)"%(nu, z))
    else:
        raise ValueError, "unknown algorithm '%s'"%algorithm

def bessel_K(nu,z,algorithm="pari",prec=53):
    r"""
    Implements the "K-Bessel function", or
    "modified Bessel function, 2nd kind", with
    index (or "order") nu and argument z. Defn:
    \begin{verbatim}
            pi*(bessel_I(-nu, z) - bessel_I(nu, z))
           ----------------------------------------
                        2*sin(pi*nu)
    \end{verbatim}

    if nu is not an integer and by taking a limit
    otherwise.

    Sometimes bessel_K(nu,z) is denoted K_nu(z) in the
    literature. In Pari, nu can be complex and
    z must be real and positive.

    EXAMPLES:
        sage: bessel_K(3,2,"scipy")
        0.64738539094...
        sage: bessel_K(3,2)
        0.64738539094...
        sage: bessel_K(1,1)
        0.60190723019...
        sage: bessel_K(1,1,"pari",10)
        0.60
        sage: bessel_K(1,1,"pari",100)
        0.60190723019723457473754000154
    """
    if algorithm=="scipy":
        if prec != 53:
            raise ValueError, "for the scipy algorithm the precision must be 53"
        import scipy.special
        ans = str(scipy.special.kv(float(nu),float(z)))
        ans = ans.replace("(","")
        ans = ans.replace(")","")
        ans = ans.replace("j","*I")
        return sage_eval(ans)
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
    else:
        raise ValueError, "unknown algorithm '%s'"%algorithm


def bessel_Y(nu,z,algorithm="maxima", prec=53):
    r"""
    Implements the "Y-Bessel function", or
    "Bessel function of the 2nd kind", with
    index (or "order") nu and argument z.

    NOTE: Currently only prec=53 is supported.

    Defn:
\begin{verbatim}
            cos(pi n)*bessel_J(nu, z) - bessel_J(-nu, z)
           -------------------------------------------------
                             sin(nu*pi)
\end{verbatim}
    if nu is not an integer and by taking a limit
    otherwise.

    Sometimes bessel_Y(n,z) is denoted Y_n(z) in the
    literature.

    This is computed using Maxima by default.

    EXAMPLES:
        sage: bessel_Y(2,1.1,"scipy")
        -1.4314714939...
        sage: bessel_Y(2,1.1)
        -1.4314714939590...
        sage: bessel_Y(3.001,2.1)
        -1.0299574976424...

    NOTE: Adding "0"+ inside sage_eval as a temporary bug work-around.
    """
    if prec != 53:
        raise ValueError, "for the scipy algorithm the precision must be 53"
    if algorithm=="scipy":
        import scipy.special
        ans = str(scipy.special.yv(float(nu),complex(real(z),imag(z))))
        ans = ans.replace("(","")
        ans = ans.replace(")","")
        ans = ans.replace("j","*I")
        return sage_eval(ans)
    elif algorithm == "maxima":
        return RR(maxima.eval("bessel_y(%s,%s)"%(float(nu),float(z))))
    else:
        raise ValueError, "unknown algorithm '%s'"%algorithm

class Bessel():
    """
    A class implementing the I, J, K, and Y Bessel functions.

    EXAMPLES:
        sage: g = Bessel(2); g
        J_{2}
        sage: print g
        J-Bessel function of order 2
        sage: g.plot(0,10)
    """
    def __init__(self, nu, typ = "J", algorithm = "pari", prec = 53):
        self._order = nu
        self._system = algorithm
        if not (typ in ['I', 'J', 'K', 'Y']):
            raise ValueError, "typ must be one of I, J, K, Y"
        self._type = typ
        prec = int(prec)
        if prec < 0:
            raise ValueError, "prec must be a positive integer"
        self._prec = int(prec)

    def __str__(self):
        return self.type()+"-Bessel function of order "+str(self.order())

    def __repr__(self):
        return self.type()+"_{"+str(self.order())+"}"

    def type(self):
        return self._type

    def prec(self):
        return self._prec

    def order(self):
        return self._order

    def system(self):
        return self._system

    def __call__(self,z):
        nu = self.order()
        t = self.type()
        s = self.system()
        p = self.prec()
        if t == "I":
            return bessel_I(nu,z,algorithm=s,prec=p)
        if t == "J":
            return bessel_J(nu,z,algorithm=s,prec=p)
        if t == "K":
            return bessel_K(nu,z,algorithm=s,prec=p)
        if t == "Y":
            return bessel_Y(nu,z,algorithm=s,prec=p)

    def plot(self,a,b):
        nu = self.order()
        s = self.system()
        t = self.type()
        if t == "I":
            f = lambda z: bessel_I(nu,z,s)
            P = plot(f,a,b)
        if t == "J":
            f = lambda z: bessel_J(nu,z,s)
            P = plot(f,a,b)
        if t == "K":
            f = lambda z: bessel_K(nu,z,s)
            P = plot(f,a,b)
        if t == "Y":
            f = lambda z: bessel_Y(nu,z,s)
            P = plot(f,a,b)
        return P

def hypergeometric_U(alpha,beta,x,algorithm="pari",prec=53):
    r"""
    Default is a wrap of Pari's hyperu(alpha,beta,x) function.
    Optionally, algorithm = "scipy" can be used.

    The confluent hypergeometric function $y = U(a,b,x)$ is defined
    to be the solution to Kummer's differential equation

    \[
    xy'' + (b-x)y' - ay = 0.
    \]
    This satisfies $U(a,b,x) \sim x^{-a}$, as $x\rightarrow \infty$,
    and is sometimes denoted \verb|x^{-a}2_F_0(a,1+a-b,-1/x)|.
    This is not the same as Kummer's $M$-hypergeometric
    function, denoted sometimes as \verb|_1F_1(alpha,beta,x)|, though
    it satisfies the same DE that $U$ does.

    WARNING:
        In the literature, both are called "Kummer confluent hypergeometric"
        functions.

    EXAMPLES:
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
    Returns the spherical Bessel function of the first kind
    for integers n > -1.

    Reference: A&S 10.1.8 page 437 and A&S 10.1.15 page 439.

    EXAMPLES:
        sage: spherical_bessel_J(2,x)
        (-(1 - 24/(8*x^2))*sin(x) - 3*cos(x)/x)/x
    """
    if algorithm=="scipy":
        import scipy.special
        ans = str(scipy.special.sph_jn(int(n),float(var)))
        ans = ans.replace("(","")
        ans = ans.replace(")","")
        ans = ans.replace("j","*I")
        return sage_eval(ans)
    elif algorithm == 'maxima':
        _init()
        return meval("spherical_bessel_j(%s,%s)"%(ZZ(n),var))
    else:
        raise ValueError, "unknown algorithm '%s'"%algorithm

def spherical_bessel_Y(n,var, algorithm="maxima"):
    r"""
    Returns the spherical Bessel function of the second kind
    for integers n > -1.

    Reference: A&S 10.1.9 page 437 and A&S 10.1.15 page 439.

    EXAMPLES:
        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: spherical_bessel_Y(2,x)
        -(3*sin(x)/x - (1 - 24/(8*x^2))*cos(x))/x
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

def spherical_hankel1(n,var):
    r"""
    Returns the spherical Hankel function of the first
    kind for integers $n > -1$, written as a string.
    Reference: A&S 10.1.36 page 439.

    EXAMPLES:
        sage: spherical_hankel1(2,'x')
        -3*I*(-x^2/3 - I*x + 1)*e^(I*x)/x^3
    """
    _init()
    return meval("spherical_hankel1(%s,%s)"%(ZZ(n),var))

def spherical_hankel2(n,x):
    r"""
    Returns the spherical Hankel function of the second
    kind for integers n > -1, written as a string.
    Reference: A&S 10.1.17 page 439.

    EXAMPLES:
        sage: spherical_hankel2(2,'x')
        '3*I*(-x^2/3+I*x+1)*%e^-(I*x)/x^3'

    Here I = sqrt(-1).

    """
    _init()
    y = str(x)
    return maxima.eval("spherical_hankel2(%s,%s)"%(n,y)).replace("%i","I")

def spherical_harmonic(m,n,x,y):
    r"""
    Returns the spherical Harmonic function of the second
    kind for integers $n > -1$, $|m|\leq n$.
    Reference: Merzbacher 9.64.

    EXAMPLES:
        sage: x,y = var('x,y')
        sage: spherical_harmonic(3,2,x,y)
        15*sqrt(7)*cos(x)*sin(x)^2*e^(2*I*y)/(4*sqrt(30)*sqrt(pi))
        sage: spherical_harmonic(3,2,1,2)
        15*sqrt(7)*e^(4*I)*cos(1)*sin(1)^2/(4*sqrt(30)*sqrt(pi))
    """
    _init()
    return meval("spherical_harmonic(%s,%s,%s,%s)"%(ZZ(m),ZZ(n),x,y))

####### elliptic functions and integrals

def jacobi(sym,x,m):
    r"""
    Here sym = "pq", where p,q in {c,d,n,s}. This returns the
    value of the Jacobi function pq(x,m), as described in the
    documentation for SAGE's "special" module. There are a
    total of 12 functions described by this.

    EXAMPLES:
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
    _init()
    return meval("jacobi_%s(%s,%s)"%(sym, x, m))

def inverse_jacobi(sym,x,m):
    """
    Here sym = "pq", where p,q in {c,d,n,s}. This returns the
    value of the inverse Jacobi function $pq^{-1}(x,m)$. There are a
    total of 12 functions described by this.

    EXAMPLES:
        sage: jacobi("sn",1/2,1/2)
        jacobi_sn(1/2, 1/2)
        sage: float(jacobi("sn",1/2,1/2))
        0.4707504736556572
        sage: float(inverse_jacobi("sn",0.47,1/2))
        0.4990982313222197
        sage: float(inverse_jacobi("sn",0.4707504,0.5))
        0.49999991146655459
        sage: P = plot(inverse_jacobi('sn', x, 0.5), 0, 1, plot_points=20)

    Now to view this, just type show(P).
    """
    _init()
    return meval("inverse_jacobi_%s(%s,%s)"%(sym, x,m))

#### elliptic integrals

def elliptic_e (phi, m):
    r"""
    This returns the value of the incomplete elliptic integral of the
    second kind,  $\int_0^\phi \sqrt(1 - m\sin(x)^2)\, dx$, ie,
    \code{integrate(sqrt(1 - m*sin(x)\^2), x, 0, phi)}. Taking
    $\phi = \pi/2$ gives \code{elliptic_ec}.

    EXAMPLES:
        sage: z = var("z")
        sage: elliptic_e (z, 1)
        sin(z)
        sage: elliptic_e (z, 0)
        z
        sage: elliptic_e (0.5, 0.1)
        0.498011394499

    """
    _init()
    return meval("elliptic_e(%s,%s)"%(phi,m))

def elliptic_ec (m):
    """
    This returns the value of the complete elliptic integral of the
    second kind,  $\int_0^{\pi/2} \sqrt(1 - m\sin(x)^2)\, dx$.

    EXAMPLES:
        sage: elliptic_ec (0.1)
        1.5307576369
        sage: elliptic_ec (x).diff()
        (elliptic_ec(x) - elliptic_kc(x))/(2*x)

    """
    _init()
    return meval("elliptic_ec(%s)"%m)


def elliptic_eu (u, m):
    r"""
    This returns the value of the incomplete elliptic integral of the
    second kind defined by $\int_0^u jacobi_dn(x,m)^2)\, dx$.

    EXAMPLES:
        sage: elliptic_eu (0.5, 0.1)
        0.496054551287

    """
    _init()
    return meval("elliptic_eu(%s,%s)"%(u,m))


def elliptic_f (phi, m):
    r"""
    This returns the value of the "incomplete elliptic integral
    of the first kind", $\int_0^\phi \frac{dx}{\sqrt{1 - m\sin(x)^2}}$,
    ie, \code{integrate(1/sqrt(1 - m*sin(x)\^2), x, 0, phi)}.
    Taking $\phi = \pi/2$ gives \code{elliptic_kc}.

    EXAMPLES:
        sage: z = var("z")
        sage: elliptic_f (z, 0)
        z
        sage: elliptic_f (z, 1)
        log(tan(z/2 + pi/4))
        sage: elliptic_f (0.2, 0.1)
        0.200132506748

    """
    _init()
    return meval("elliptic_f(%s,%s)"%(phi,m))

def elliptic_kc (m):
    r"""
    This returns the value of the "complete elliptic integral
    of the first kind", $\int_0^{\pi/2} \frac{dx}{\sqrt{1 - m\sin(x)^2}}$.

    EXAMPLES:
        sage: elliptic_kc (0.5)
        1.8540746773
        sage: elliptic_f (RR(pi/2), 0.5)
        1.8540746773

    """
    _init()
    return meval("elliptic_kc(%s)"%m)

def elliptic_pi (n, phi, m):
    r"""
    This returns the value of the "incomplete elliptic integral
    of the third kind",
    $$\int_0^\phi \frac{dx}{\sqrt{(1 - m\sin(x)^2)(1 - n\sin(x)^2)}}$$.

    EXAMPLES:
        sage: elliptic_pi(0.1, 0.2, 0.3)
        0.200665068221

    """
    _init()
    return meval("elliptic_pi(%s,%s,%s)"%(n,phi,m))

#### hyperboic trig functions (which are special cases
#### of Jacobi elliptic functions but faster to evaluate directly)

## def sinh(t):
##     try:
##         return t.sinh()
##     except AttributeError:
##         from sage.calculus.calculus import exp
##         return (exp(t)-exp(-t))/2

## def cosh(t):
##     try:
##         return t.cosh()
##     except AttributeError:
##         from sage.calculus.calculus import exp
##         return (exp(t)+exp(-t))/2

## def coth(t):
##     try:
##         return t.coth()
##     except AttributeError:
##         return 1/tanh(t)

## def sech(t):
##     try:
##         return t.sech()
##     except AttributeError:
##         return 1/cosh(t)

## def csch(t):
##     try:
##         return t.csch()
##     except AttributeError:
##         return 1/sinh(t)

## Now implemented using polylog in calculus.py:
## def dilog(t):
##    """
##    Te dilogarithm of t is the analytic continuation of the
##    power series $\sum_{n \geq 1} t^n/n^2$.
##    """
##    try:
##        return t.dilog()
##    except AttributeError:
##        raise NotImplementedError

def lngamma(t):
    """
    The principal branch of the logarithm of the Gamma function of t.
    """
    try:
        return t.lngamma()
    except AttributeError:
        raise NotImplementedError

def exp_int(t):
    r"""
    The exponential integral $\int_x^\infty e^{-x}/x dx$ (t belongs to RR).
    """
    try:
        return t.eint1()
    except AttributeError:
        raise NotImplementedError

def error_fcn(t):
    r"""
    The complementary error function
    $\frac{2}{\sqrt{\pi}}\int_t^\infty e^{-x^2} dx$ (t belongs to RR).
    """
    try:
        return t.erfc()
    except AttributeError:
        raise NotImplementedError



